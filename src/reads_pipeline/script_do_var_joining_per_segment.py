import sys
import tempfile
from pathlib import Path
from functools import partial
import multiprocessing

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import (
    do_svn_joint_genotyping_for_all_samples_together,
)


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)
    config_path = get_config_path(project_dir)

    if not project_dir.exists():
        msg = f"The project directory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    if not config_path.exists():
        msg = f"The config file {config_path} does not exist"
        print(msg)
        sys.exit(3)

    config = PipelineConfig(project_dir=project_dir)

    genome_path = config["general"]["genome_path"]
    n_processes = config["snv_calling"]["gatk_segment_vcf_sample_joining_n_process"]
    filters = config["gatk_filters"]

    do_svn_joint_genotyping_for_all_samples_together(
        project_dir=project_dir,
        genome_fasta=genome_path,
        n_processes=n_processes,
        gatk_filters=filters,
    )


if __name__ == "__main__":
    main()
