import sys

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import do_snv_calling_per_sample


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

    do_snv_calling_per_sample(
        project_dir=project_dir,
        re_run=config["general"]["re_run"],
        verbose=config["general"]["verbose"],
        num_snvs_in_parallel=config["general"]["num_snvs_in_parallel"],
        min_mapq=config["gatk"]["per_sample_calling_min_mapq"],
        genome_fasta=config["general"]["genome_path"],
        tabix_num_threads=config["gatk"]["num_threads_tabix"],
    )


if __name__ == "__main__":
    main()
