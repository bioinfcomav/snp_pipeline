from pathlib import Path
from collections import defaultdict
import sys
import os
import argparse

from reads_pipeline.paths import (
    get_project_dir,
    get_raw_reads_parent_dir,
    get_reads_stats_fastp_parent_dir,
    get_paired_and_unpaired_read_files_in_dir,
    get_read_group_info_xls,
)
from reads_pipeline.fastp_minimap import _get_read_group_id_from_fastq_pair_paths
from reads_pipeline.read_group import get_read_group_info


def print_fastq_pair(fastq_pair):
    print(",".join(map(str, fastq_pair)))


def check_read_groups_to_map(project_dir: Path):
    read_groups_info = get_read_group_info(project_dir)

    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)
    stats_parent_dir = get_reads_stats_fastp_parent_dir(project_dir)
    stats_parent_dir.mkdir(exist_ok=True, parents=True)
    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]

    fastq_pairs_by_read_group = defaultdict(list)
    for raw_reads_dir in raw_reads_dirs:
        for fastq_pair in get_paired_and_unpaired_read_files_in_dir(raw_reads_dir):
            read_group_id = _get_read_group_id_from_fastq_pair_paths(fastq_pair)
            fastq_pairs_by_read_group[read_group_id].append(fastq_pair)

    print(f"Number of read groups: {len(fastq_pairs_by_read_group)}")

    for read_group, fastq_pairs in fastq_pairs_by_read_group.items():
        if len(fastq_pairs) > 1:
            print(f"ERROR: more than one fastq pairs for read group: {read_group}")
            for fastq_pair in fastq_pairs:
                print_fastq_pair(fastq_pair)
        elif not fastq_pairs:
            print_fastq_pair(f"ERROR: no fastq pair for read group: {read_group}")
        else:
            fastq_pair = fastq_pairs[0]
            if len(fastq_pair) > 2:
                print(
                    f"ERROR: fastq pair has more than two fastq files for read group: {read_group}"
                )
                print_fastq_pair(fastq_pair)

    read_group_info_excel = get_read_group_info_xls(project_dir)
    for read_group in fastq_pairs_by_read_group.keys():
        if read_group not in read_groups_info:
            print(
                f"ERROR: read group {read_group} not defined in read group excel file: {read_group_info_excel}"
            )


def get_args():
    parser = argparse.ArgumentParser(
        description="Check the read groups and fastq pairs"
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

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    check_read_groups_to_map(project_dir)
