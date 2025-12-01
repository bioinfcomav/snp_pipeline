import os
import argparse
import tomllib
import sys

from reads_pipeline.paths import get_config_path, get_project_dir
from reads_pipeline.fastp_minimap import run_fastp_minimap_for_fastqs
from reads_pipeline.pipeline_config import PipelineConfig


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

    config = PipelineConfig(project_dir=project_dir)

    run_fastp_minimap_for_fastqs(
        project_dir=project_dir,
        minimap_index=config["minimap"]["index_path"],
        genome_fasta=config["general"]["genome_path"],
        deduplicate=config["samtools"]["deduplicate"],
        min_read_len=config["fastp"]["min_read_len"],
        fastp_num_threads=config["fastp"]["num_threads"],
        fastp_trim_front1=config["fastp"]["trim_front1"],
        fastp_trim_tail1=config["fastp"]["trim_tail1"],
        fastp_trim_front2=config["fastp"]["trim_front2"],
        fastp_trim_tail2=config["fastp"]["trim_tail2"],
        minimap_num_threads=config["minimap"]["num_threads"],
        sort_num_threads=config["samtools"]["sort_num_threads"],
        duplicates_num_threads=config["samtools"]["duplicates_num_threads"],
        calmd_num_threads=config["samtools"]["calmd_num_threads"],
        samtools_stats_num_threads=config["samtools"]["samtools_stats_num_threads"],
        re_run=config["general"]["re_run"],
        verbose=config["general"]["verbose"],
        num_mappings_in_parallel=config["general"]["num_mappings_in_parallel"],
        cmd1=config["mapping_command_hooks"]["cmd1"],
        force_cram_version=config["samtools"]["force_cram_version"],
    )


if __name__ == "__main__":
    main()
