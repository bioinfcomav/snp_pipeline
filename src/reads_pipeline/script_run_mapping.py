import os
import argparse
import tomllib
import sys

from reads_pipeline.paths import get_config_path, get_project_dir
from reads_pipeline.genomic_pipeline import PipelineConfig, run_pipeline


def get_args():
    parser = argparse.ArgumentParser(
        description="Cleans and aligns reads using fastp and minimap2."
    )

    # project_dir: Defaults to the current working directory
    parser.add_argument(
        "project_dir",
        nargs="?",  # Makes it optional
        default=os.getcwd(),
        help="Project directory (default: current working directory).",
    )

    return parser.parse_args()


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
        sys.exit(2)

    with config_path.open("rb") as fhand:
        config = PipelineConfig(tomllib.load(fhand))

    run_pipeline(project_dir, config)


if __name__ == "__main__":
    main()
