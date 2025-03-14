from pathlib import Path
import os
import logging
import gzip
from collections import defaultdict
import tempfile
import json
from enum import Enum

from reads_pipeline.run_cmd import run_cmd
from reads_pipeline.paths import (
    get_project_dir,
    SAMTOOLS_BIN,
    GATK_PYTHON_BIN,
    get_tmp_dir,
    get_gatk_db_dir,
)


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


def filter_and_merge_variants():
    """
# Filter Out Low-Quality Variants (Hard Filtering)

Since you cannot use VariantRecalibrator, you will apply hard filtering with VariantFiltration.
Filtering SNPs

gatk VariantFiltration \
    -R reference.fasta \
    -V raw_variants.vcf.gz \
    -O filtered_SNPs.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum"

Explanation of Filters

    Quality by Depth (QD < 2.0) → Low confidence variants.
    Fisher Strand Bias (FS > 60.0) → Highly biased variants.
    Mapping Quality (MQ < 40.0) → Low mapping confidence.
    Mapping Quality Rank Sum (MQRankSum < -12.5) → Bad read alignments.
    Read Position Bias (ReadPosRankSum < -8.0) → Potential alignment artifacts.

This step removes low-confidence SNPs based on their quality metrics.
Filtering Indels

Indels often have different quality metrics, so we use slightly different criteria:

gatk VariantFiltration \
    -R reference.fasta \
    -V filtered_SNPs.vcf.gz \
    -O filtered_variants.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum"

    FS cutoff is higher (FS > 200.0) because indels are more prone to strand bias.
    Read position bias cutoff is stricter (ReadPosRankSum < -20.0).

Now, you have a high-confidence filtered VCF.

# Merge SNPs and Indels into Complex Variants

To combine adjacent SNPs and indels into haplotypes and complex alleles, use bcftools norm:

bcftools norm -m +any -f reference.fasta filtered_variants.vcf.gz -O z -o final_variants.vcf.gz

    --multiallelics+both → join biallelic sites into multiallelic records, merges SNPs and indels into complex variants.
    --fasta-ref reference.fasta → Uses the reference genome to correct alleles.
    --check-ref e
    --collapse both
    --output FILE
    --output-type z → compressed VCF (.vcf.gz)
    --threads INT
    --write-index

# Validate the Final VCF

Validate the VCF to ensure correct formatting

bcftools validate final_variants.vcf.gz
"""
