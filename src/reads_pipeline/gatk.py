from pathlib import Path
import os
import logging

from reads_pipeline.run_cmd import run_cmd
from reads_pipeline.paths import (
    get_project_dir,
    SAMTOOLS_BIN,
    GATK_PYTHON_BIN,
    get_tmp_dir,
)


logger = logging.getLogger(__name__)


def create_genome_reference(
    genome_fasta: Path, reference_out_dir: Path, project_dir: Path
):
    project_dir = get_project_dir(project_dir)
    genome_name = genome_fasta.name
    genome_path = reference_out_dir / genome_name
    os.symlink(genome_fasta, genome_path)

    cmd = [SAMTOOLS_BIN, "faidx", str(genome_path)]
    run_cmd(cmd, project_dir=project_dir)

    gatk_dict_path = genome_path.with_suffix(".dict")

    cmd = ["uv", "run", str(GATK_PYTHON_BIN), "CreateSequenceDictionary"]
    cmd.append(f"R={genome_path}")
    cmd.append(f"O={gatk_dict_path}")
    run_cmd(cmd, project_dir=project_dir)
    return {"genome_path": genome_path, "gatk_dict_path": gatk_dict_path}


def do_sample_snv_calling_basic_germline(
    genome_fasta: Path,
    bam: Path,
    out_vcf: Path,
    project_dir,
    min_mapq: int = 10,
    mem_in_gb=4,
):
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
    cmd.extend(["-I", str(bam)])
    cmd.extend(["-O", str(out_vcf)])
    cmd.extend(["--mapping-quality-threshold-for-genotyping", str(min_mapq)])
    cmd.extend(["-ERC", "BP_RESOLUTION"])
    run_cmd(cmd, project_dir=project_dir)


def create_db_with_sample_snv_calls():
    pass


def do_global_svn_calling():
    pass
