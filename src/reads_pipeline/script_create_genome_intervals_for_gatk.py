import sys

from reads_pipeline.script_run_mapping import get_project_dir, get_config_path, get_args
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import (
    create_gatk_intervals_file_from_chromosomes,
    get_genome_fai_path,
)


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

    fai_path = get_genome_fai_path(config)

    create_gatk_intervals_file_from_chromosomes(
        project_dir=project_dir, genome_fai_path=fai_path
    )


if __name__ == "__main__":
    main()
