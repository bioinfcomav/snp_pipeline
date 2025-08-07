import sys
from pathlib import Path

from split_gvcf.split_genome_avoiding_vars import create_var_ranges

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.paths import get_vcfs_per_sample_dir, get_gvcf_ranges_working_dir


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

    vcfs = [
        path
        for path in get_vcfs_per_sample_dir(project_dir).iterdir()
        if str(path).endswith(".vcf.gz")
    ]
    _ranges = create_var_ranges(
        vcfs,
        get_gvcf_ranges_working_dir(project_dir),
        config["split_gvcf_vars"]["n_processes_gvcf_parsing"],
    )


if __name__ == "__main__":
    main()
