import sys
import argparse
import os

import pandas

from reads_pipeline.paths import (
    get_project_dir,
    get_reads_stats_fastp_excel_report_path,
    get_crams_stats_dir,
    get_crams_stats_excel_report_path,
)
from reads_pipeline.fastp_minimap import get_fastq_pairs_to_process, collect_cram_stats
from reads_pipeline.fastp import collect_fastp_stats
from reads_pipeline.read_group import get_read_group_info


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


def check_num_clean_reads_vs_mapped_reads(fastp_stats, cram_stats):
    col1 = fastp_stats["num_clean_reads"]
    col2 = cram_stats["total sequences"]
    if not all(col1 == col2):
        combined = pandas.concat(
            [col1, col2],
            axis=1,
            keys=["s1", "s2"],
        )
        diff_any = combined[combined["s1"] != combined["s2"]]
        msg = "Not matching rows when comparing fastp clean reads and cram total sequences results:\n"
        msg += str(diff_any)
        raise RuntimeError(msg)


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    read_groups_info = get_read_group_info(project_dir)
    pairs_to_process = get_fastq_pairs_to_process(project_dir, read_groups_info)
    print(f"Num. pairs to process: {len(pairs_to_process)}")

    fastp_stats = collect_fastp_stats(project_dir=project_dir)
    fastp_stats.to_excel(
        get_reads_stats_fastp_excel_report_path(project_dir), index=False
    )
    print(f"Num. fastp stats: {fastp_stats.shape[0]}")

    cram_stats = collect_cram_stats(project_dir=project_dir)
    stats_dir = get_crams_stats_dir(project_dir)
    stats_dir.mkdir(exist_ok=True)
    stats_path = get_crams_stats_excel_report_path(project_dir)
    cram_stats.to_excel(stats_path, index=False)
    print(f"Num. cram stats: {cram_stats.shape[0]}")

    index = [
        (fastp_stats_row["dir"], fastp_stats_row["file_name"].split(".")[0])
        for _, fastp_stats_row in fastp_stats.iterrows()
    ]
    fastp_stats.index = index
    index = [
        (stats_row["dir"], stats_row["file_name"].split(".")[0])
        for _, stats_row in cram_stats.iterrows()
    ]
    cram_stats.index = index
    check_num_clean_reads_vs_mapped_reads(fastp_stats, cram_stats)
