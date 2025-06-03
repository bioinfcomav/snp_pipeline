from pathlib import Path
import os
import logging
import gzip
from collections import defaultdict
import tempfile
import json
from enum import Enum
import shutil
from functools import partial
from multiprocessing import Pool


from reads_pipeline.run_cmd import run_cmd
from reads_pipeline.paths import (
    get_project_dir,
    SAMTOOLS_BIN,
    GATK_PYTHON_BIN,
    get_tmp_dir,
    get_gatk_db_dir,
    BCFTOOLS_BIN,
    get_crams_dir,
    get_vcfs_per_sample_dir,
    TABIX_BIN,
)
from reads_pipeline.read_group import get_read_group_info, get_read_group_id_from_path

logger = logging.getLogger(__name__)


def create_faidx(genome_path: Path, project_dir: Path):
    project_dir = get_project_dir(project_dir)
    cmd = [SAMTOOLS_BIN, "faidx", str(genome_path)]
    run_cmd(cmd, project_dir=project_dir)
    fai_path = Path(str(genome_path) + ".fai")
    if not fai_path.exists():
        raise ValueError(f"fai genome file not craeted for {genome_path}")
    return {"fai_path": fai_path}


def get_chrom_lengths_from_fai(genome_fai_path):
    chrom_lengths = {}
    with open(genome_fai_path) as fhand:
        for line in fhand:
            chrom, length, *_ = line.strip().split("\t")
            chrom_lengths[chrom] = int(length)
    return chrom_lengths


def create_genome_reference(
    genome_fasta: Path, reference_out_dir: Path, project_dir: Path
):
    project_dir = get_project_dir(project_dir)
    genome_name = genome_fasta.name
    genome_path = reference_out_dir / genome_name
    os.symlink(genome_fasta, genome_path)

    fai_path = create_faidx(genome_path, project_dir)["fai_path"]

    gatk_dict_path = genome_path.with_suffix(".dict")

    cmd = ["uv", "run", str(GATK_PYTHON_BIN), "CreateSequenceDictionary"]
    cmd.append(f"R={genome_path}")
    cmd.append(f"O={gatk_dict_path}")
    run_cmd(cmd, project_dir=project_dir)
    return {
        "genome_path": genome_path,
        "gatk_dict_path": gatk_dict_path,
        "fai_path": fai_path,
    }


def do_sample_snv_calling_basic_germline(
    genome_fasta: Path,
    bams: list[Path],
    out_vcf: Path,
    project_dir,
    min_mapq: int = 10,
    mem_in_gb=4,
    allow_uncompressed_vcf=False,
):
    if not out_vcf.suffix == ".gz" and not allow_uncompressed_vcf:
        raise ValueError("Output VCF must have a .gz suffix")

    fasta_dict_path = genome_fasta.with_suffix(".dict")
    if not fasta_dict_path.exists():
        raise ValueError(
            f"Missing dict file for genome fasta, please create it with GATK: {fasta_dict_path}"
        )

    tmp_dir = get_tmp_dir(project_dir)
    tmp_dir.mkdir(exist_ok=True)
    cmd = [
        "uv",
        "run",
        str(GATK_PYTHON_BIN),
        "--java-options",
        f"-Xmx{mem_in_gb}g",
        "HaplotypeCaller",
        "--tmp-dir",
        str(tmp_dir),
    ]
    cmd.extend(["-R", str(genome_fasta)])
    for bam in bams:
        cmd.extend(["-I", str(bam)])
    cmd.extend(["-O", str(out_vcf)])
    cmd.extend(["--mapping-quality-threshold-for-genotyping", str(min_mapq)])
    cmd.extend(["-ERC", "BP_RESOLUTION"])
    run_cmd(cmd, project_dir=project_dir)


def get_crams(project_dir):
    base_crams_dir = get_crams_dir(project_dir)

    cram_paths = []
    for path in base_crams_dir.iterdir():
        if not path.is_dir():
            continue
        cram_paths.extend(path.glob("*.cram"))
    return cram_paths


def create_tabix_path_from_vcf(vcf_path):
    return Path(str(vcf_path) + ".tbi")


def _index_vcf(vcf_path, project_dir, num_threads=0, re_run=False):
    tabix_path = create_tabix_path_from_vcf(vcf_path)
    if tabix_path.exists():
        if re_run:
            os.remove(tabix_path)
        else:
            return tabix_path

    cmd = [TABIX_BIN]
    if num_threads:
        cmd.append(f"-@{num_threads}")
    cmd.append(str(vcf_path))
    run_cmd(cmd, project_dir=project_dir, verbose=True)

    if not tabix_path.exists():
        raise RuntimeError("tabix indexing run but tbi not found: {tabix_path}")
    return tabix_path


def _do_snv_calling_for_sample(
    sample_info,
    vcfs_per_sample_dir,
    genome_fasta,
    project_dir,
    min_mapq,
    tabix_num_threads,
    re_run,
):
    cram_paths = sample_info["cram_paths"]
    out_vcf = sample_info["out_vcf"]
    with tempfile.TemporaryDirectory(
        prefix="gatk_per_sample", dir=vcfs_per_sample_dir
    ) as tmp_dir:
        out_tmp_vcf = Path(tmp_dir) / out_vcf.name
        do_sample_snv_calling_basic_germline(
            genome_fasta=genome_fasta,
            bams=cram_paths,
            out_vcf=out_tmp_vcf,
            project_dir=project_dir,
            min_mapq=min_mapq,
        )
        vcf_tmp_index = _index_vcf(
            out_tmp_vcf, project_dir, num_threads=tabix_num_threads, re_run=re_run
        )
        vcf_index = create_tabix_path_from_vcf(out_vcf)
        shutil.move(out_tmp_vcf, out_vcf)
        shutil.move(vcf_tmp_index, vcf_index)


def do_snv_calling_per_sample(
    project_dir: Path,
    genome_fasta: Path,
    min_mapq: int = 10,
    verbose=False,
    re_run=False,
    num_snvs_in_parallel=1,
    tabix_num_threads=0,
):
    vcfs_per_sample_dir = get_vcfs_per_sample_dir(project_dir)
    vcfs_per_sample_dir.mkdir(parents=True, exist_ok=True)

    group_infos = get_read_group_info(project_dir)

    cram_paths = get_crams(project_dir)

    infos_per_sample = {}
    for cram_path in cram_paths:
        read_group_id = get_read_group_id_from_path(cram_path)
        if read_group_id not in group_infos:
            raise ValueError(f"Missing read group info for {read_group_id}")
        sample = group_infos[read_group_id]["sample"]
        if sample not in infos_per_sample:
            infos_per_sample[sample] = {"cram_paths": []}
        infos_per_sample[sample]["cram_paths"].append(cram_path)

    for sample in infos_per_sample.keys():
        out_vcf = vcfs_per_sample_dir / f"{sample}.g.vcf.gz"
        if out_vcf.exists():
            if re_run:
                os.remove(out_vcf)
            else:
                continue
        infos_per_sample[sample]["out_vcf"] = out_vcf

    if verbose:
        print(f"Num. samples: {len(infos_per_sample)}")
    sample_infos_todo = {
        sample: sample_info
        for sample, sample_info in infos_per_sample.items()
        if "out_vcf" in sample_info
    }
    if verbose:
        print(f"Num. samples/SNV callings to do: {len(sample_infos_todo)}")

    do_snv_calling_for_sample = partial(
        _do_snv_calling_for_sample,
        vcfs_per_sample_dir=vcfs_per_sample_dir,
        genome_fasta=genome_fasta,
        project_dir=project_dir,
        min_mapq=min_mapq,
        tabix_num_threads=tabix_num_threads,
        re_run=re_run,
    )

    if num_snvs_in_parallel > 1:
        with Pool(num_snvs_in_parallel) as pool:
            res = pool.map(do_snv_calling_for_sample, sample_infos_todo.values())
    else:
        res = map(do_snv_calling_for_sample, sample_infos_todo.values())
    list(res)

    return {"samples_done": sample_infos_todo}


def _get_sample_names_from_vcf(vcf):
    with gzip.open(vcf, "rt") as fhand:
        for line in fhand:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                return samples


class GATKDBFileMode(Enum):
    CREATE = "create"
    UPDATE = "update"


def create_db_with_independent_sample_snv_calls(
    vcfs: list[Path],
    project_dir: Path,
    genome_fai_path: Path,
    mode: GATKDBFileMode,
    batch_size: int = 50,
    reader_threads: int = 4,
):
    chrom_lengths = get_chrom_lengths_from_fai(genome_fai_path)
    genome_bed_path = genome_fai_path.with_suffix(".bed")
    with open(genome_bed_path, "wt") as fhand:
        for chrom, length in chrom_lengths.items():
            fhand.write(f"{chrom}\t0\t{length}\n")
        fhand.flush()

    if mode == GATKDBFileMode.UPDATE:
        samples_in_db = get_samples_in_gatk_db(project_dir)
    else:
        samples_in_db = []

    vcfs_per_sample = defaultdict(set)
    for vcf in vcfs:
        samples = _get_sample_names_from_vcf(vcf)
        if len(samples) > 1:
            raise ValueError(
                "VCF must contain only one sample VCF {vcf} has these samples: {samples}"
            )
        elif len(samples) == 0:
            raise ValueError(f"VCF must contain at least one sample {vcf}")

        if samples[0] in samples_in_db:
            raise ValueError(
                f"Sample {samples[0]} is already in the GATK DB, it should not be in another VCF"
            )

        vcfs_per_sample[samples[0]].add(vcf)

    for sample, vcfs in vcfs_per_sample.items():
        if len(vcfs) > 1:
            raise ValueError(
                "Every sample should be in just one VCF, but {sample} is in {vcfs}, otherwise the SNV calling accuracy would be degraded"
            )

    db_dir = get_gatk_db_dir(project_dir)
    if mode == GATKDBFileMode.CREATE:
        if db_dir.exists():
            raise ValueError(f"DB dir already exists, GATK would fail: {db_dir}")
    else:
        if not db_dir.exists():
            raise ValueError(f"DB dir does not exist: {db_dir}")

    genomic_cmd = (
        "--genomicsdb-update-workspace-path"
        if mode == GATKDBFileMode.UPDATE
        else "--genomicsdb-workspace-path"
    )

    tmp_dir = get_tmp_dir(project_dir)
    tmp_dir.mkdir(exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wt", prefix="samples_", suffix=".map"
    ) as samples_map:
        for sample, vcfs in vcfs_per_sample.items():
            samples_map.write(f"{sample}\t{list(vcfs)[0]}\n")
        samples_map.flush()
        cmd = [
            "uv",
            "run",
            str(GATK_PYTHON_BIN),
            "GenomicsDBImport",
            genomic_cmd,
            str(db_dir),
            "--batch-size",
            str(batch_size),
            "--sample-name-map",
            samples_map.name,
            "--reader-threads",
            str(reader_threads),
            "-bypass-feature-reader",
            "-L",
            str(genome_bed_path),
            # "--java-options",
            # "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true",
        ]
        run_cmd(cmd, project_dir=project_dir)


def get_samples_in_gatk_db(project_dir):
    db_dir = get_gatk_db_dir(project_dir)
    if not db_dir.exists():
        raise ValueError(f"DB dir does not exist: {db_dir}")
    json_path = db_dir / "callset.json"
    if not json_path.exists():
        raise ValueError(f"DB dir does not contain callset.json: {db_dir}")

    json_data = json.loads(json_path.open("rt").read())
    samples = [callset["sample_name"] for callset in json_data["callsets"]]
    return samples


def do_svn_joint_genotyping_for_all_samples_together(
    project_dir, genome_fasta: Path, out_vcf: Path, allow_uncompressed_vcf=False
):
    if not out_vcf.suffix == ".gz" and not allow_uncompressed_vcf:
        raise ValueError("Output VCF must have a .gz suffix")

    db_dir = get_gatk_db_dir(project_dir)
    cmd = [
        "uv",
        "run",
        str(GATK_PYTHON_BIN),
        "GenotypeGVCFs",
        "--reference",
        str(genome_fasta),
        "--variant",
        f"gendb://{db_dir}",
        "--output",
        str(out_vcf),
        "--add-output-vcf-command-line",
    ]
    run_cmd(cmd, project_dir=project_dir)


def filter_vcf_with_gatk(in_vcf, out_vcf, genome_fasta, filters, project_dir):
    cmd = [
        "uv",
        "run",
        str(GATK_PYTHON_BIN),
        "VariantFiltration",
        "--reference",
        str(genome_fasta),
        "--variant",
        str(in_vcf),
        "--output",
        str(out_vcf),
    ]
    for filter_name, filter_expression in filters.items():
        cmd.extend(
            (
                "--filter-expression",
                filter_expression,
                "--filter-name",
                filter_name,
            )
        )
    run_cmd(cmd, project_dir=project_dir)
