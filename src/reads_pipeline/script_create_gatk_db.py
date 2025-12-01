import sys

from reads_pipeline.script_run_mapping import get_project_dir, get_config_path, get_args
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import (
    create_db_with_independent_sample_snv_calls,
    GATKDBFileMode,
    get_samples_in_gatk_db,
)
from reads_pipeline.paths import get_per_sample_vcfs


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

    db_mode = config["gatk"]["db_mode"]
    if db_mode == "create":
        mode = GATKDBFileMode.CREATE
    elif db_mode == "update":
        mode = GATKDBFileMode.UPDATE
    else:
        raise ValueError(
            f"GATK db_mode should be either create or update but it is: {db_mode}"
        )

    vcfs_per_sample = create_db_with_independent_sample_snv_calls(
        vcfs=get_per_sample_vcfs(project_dir),
        project_dir=project_dir,
        mode=mode,
        batch_size=config["gatk"]["db_creation_batch_size"],
        reader_threads=config["gatk"]["db_creation_reader_threads"],
        n_gatk_db_interval_creations_in_parallel=config["gatk"][
            "n_gatk_db_interval_creations_in_parallel"
        ],
    )
    samples_in_db = get_samples_in_gatk_db(project_dir)
    print(f"Num samples added to db: {len(vcfs_per_sample)}")
    print(f"Num samples in db: {len(samples_in_db)}")


if __name__ == "__main__":
    main()
