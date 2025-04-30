import sys
import tempfile
from pathlib import Path

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import (
    do_svn_joint_genotyping_for_all_samples_together,
    filter_vcf_with_gatk,
)
from reads_pipeline.paths import get_joint_vcf, get_snv_dir


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)
    config_path = get_config_path(project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    if not config_path.exists():
        msg = f"The config file {config_path} does not exist"
        print(msg)
        sys.exit(3)

    config = PipelineConfig(project_dir=project_dir)

    genome_path = config["general"]["genome_path"]

    out_vcf = get_joint_vcf(project_dir)

    if not str(out_vcf).endswith(".vcf.gz"):
        raise RuntimeError("Joint filtered vcf should end with .vcf.gz")

    if out_vcf.exists():
        raise RuntimeError(f"Joint filtered VCF already exists: {out_vcf}")

    join_vcf_tmp = tempfile.NamedTemporaryFile(
        prefix="gatk.join.", suffix=".vcf.gz", dir=get_snv_dir(project_dir)
    )

    do_svn_joint_genotyping_for_all_samples_together(
        project_dir=project_dir,
        genome_fasta=genome_path,
        out_vcf=Path(join_vcf_tmp.name),
    )

    filter_vcf_with_gatk(
        in_vcf=Path(join_vcf_tmp.name),
        out_vcf=out_vcf,
        genome_fasta=genome_path,
        filters=config["gatk_filters"],
        project_dir=project_dir,
    )
    join_vcf_tmp.close()


if __name__ == "__main__":
    main()
