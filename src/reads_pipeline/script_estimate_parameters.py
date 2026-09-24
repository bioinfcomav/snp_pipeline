import os
import argparse
import sys

from reads_pipeline.paths import get_config_path, get_project_dir
from reads_pipeline.parameter_estimation import estimate_parameters
from reads_pipeline.pipeline_config import PipelineConfig


def get_args():
    parser = argparse.ArgumentParser(
        description="Estimates the parameters used as priors in the SNP calling."
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
    inbreeding = config["pop_var_caller"]["inbreeding"]

    estimate_parameters(
        project_dir=project_dir,
        genome_fasta=config["general"]["genome_path"],
        catalog_path=catalog_path if catalog_path else None,
        ploidy=config["pop_var_caller"]["ploidy"],
        # Saying nothing is not the same as saying that the samples are not
        # inbred, an empty inbreeding means use the fitted coefficient
        inbreeding=float(inbreeding) if inbreeding != "" else None,
        declare_batches=config["pop_var_caller"]["declare_sequencing_batches"],
        num_threads=config["pop_var_caller"]["parameters_num_threads"],
        re_run=config["general"]["re_run"],
        verbose=config["general"]["verbose"],
    )


if __name__ == "__main__":
    main()
