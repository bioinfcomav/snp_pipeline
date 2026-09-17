import os
import argparse
import sys

from reads_pipeline.paths import get_config_path, get_project_dir
from reads_pipeline.psp import generate_psps_for_samples
from reads_pipeline.pipeline_config import PipelineConfig


def get_args():
    parser = argparse.ArgumentParser(
        description="Generates the per sample pileup (psp) files using pop_var_caller."
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

    config = PipelineConfig(project_dir=project_dir)

    catalog_path = config["pop_var_caller"]["catalog_path"]
    regions_path = config["pop_var_caller"]["regions_path"]

    generate_psps_for_samples(
        project_dir=project_dir,
        genome_fasta=config["general"]["genome_path"],
        catalog_path=catalog_path if catalog_path else None,
        regions_path=regions_path if regions_path else None,
        min_copies=config["pop_var_caller"]["min_copies"],
        min_period=config["pop_var_caller"]["min_period"],
        max_period=config["pop_var_caller"]["max_period"],
        max_str_len=config["pop_var_caller"]["max_str_len"],
        min_purity=config["pop_var_caller"]["min_purity"],
        build_index_if_missing=config["pop_var_caller"]["build_index_if_missing"],
        re_run=config["general"]["re_run"],
        verbose=config["general"]["verbose"],
        num_psps_in_parallel=config["general"]["num_psps_in_parallel"],
    )


if __name__ == "__main__":
    main()
