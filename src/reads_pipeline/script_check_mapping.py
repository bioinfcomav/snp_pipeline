import argparse
import os
import sys
from pathlib import Path

import pandas

from reads_pipeline.paths import (
    get_paired_and_unpaired_read_files_in_dir,
    get_project_dir,
    get_raw_reads_parent_dir,
    is_gzip_file,
    get_raw_reads_fastqc_stats_parent_dir,
    FASTQC_XLS_STATS_FNAME,
)
from reads_pipeline.run_cmd import run_bash_script


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


COUNT_SEQS_SCRIPT = """
set -e -o pipefail
{cat_bin} {fastq} | wc -l
"""


def count_seqs_in_fastq(fastq: Path, project_dir: Path):
    cat_bin = "zcat" if is_gzip_file(fastq) else "cat"
    script = COUNT_SEQS_SCRIPT.format(cat_bin=cat_bin, fastq=fastq)
    process = run_bash_script(script, project_dir=project_dir)["process"]
    num_lines = int(process.stdout.decode().strip())
    if num_lines % 4 != 0:
        msg = f"Number of lines in fastq file {fastq} is not a multiple of 4"
        raise RuntimeError(msg)
    return num_lines // 4


def count_seqs_in_raw_reads_dir(project_dir):
    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)
    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]
    num_raw_reads_per_pair = {}
    for raw_reads_dir in raw_reads_dirs:
        bioproject_dir_name = raw_reads_dir.name
        for pair in get_paired_and_unpaired_read_files_in_dir(raw_reads_dir):
            num_reads_in_pair = [
                count_seqs_in_fastq(file, project_dir) for file in pair
            ]
            if not all(
                num_reads == num_reads_in_pair[0] for num_reads in num_reads_in_pair
            ):
                msg = f"Number of reads in pair {pair} do not match"
                raise RuntimeError(msg)
            num_reads_in_pair = num_reads_in_pair[0]
            pair = sorted(pair)
            pair_name = pair[0].name.split(".")[0]
            pair_key = (bioproject_dir_name, pair_name)
            num_raw_reads_per_pair[pair_key] = num_reads_in_pair
    return num_raw_reads_per_pair


def check_fastqc_stats(project_dir, num_raw_reads_per_pair):
    fastqc_report_xls = (
        get_raw_reads_fastqc_stats_parent_dir(project_dir) / FASTQC_XLS_STATS_FNAME
    )
    if not fastqc_report_xls.exists():
        msg = f"The fastqc report {fastqc_report_xls} does not exist"
        print(msg)
        sys.exit(3)
    fastqc_stats = pandas.read_excel(fastqc_report_xls)
    keys_in_fastqc_report = set()
    for _, fastq_file_stats in fastqc_stats.iterrows():
        pair_name = fastq_file_stats["file_name"].split(".")[0]
        bioproject = fastq_file_stats["dir"]
        if (
            fastq_file_stats["num_seqs"]
            != num_raw_reads_per_pair[(bioproject, pair_name)]
        ):
            msg = f"Number of reads in fastq file {fastq_file_stats['file_name']} according to fastqc does not match the number of reads counted in the raw fastq files"
            raise RuntimeError(msg)
        keys_in_fastqc_report.add((bioproject, pair_name))

    if keys_in_fastqc_report != set(num_raw_reads_per_pair.keys()):
        print("keys in raw reads dirs:", sorted(set(num_raw_reads_per_pair.keys())))
        print("keys in fastqc report:", sorted(keys_in_fastqc_report))
        msg = "fastq files in fastqc report do not match the fastqc in the raw reads directories"
        raise RuntimeError(msg)


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    num_raw_reads_per_pair = count_seqs_in_raw_reads_dir(project_dir)

    check_fastqc_stats(project_dir, num_raw_reads_per_pair)


if __name__ == "__main__":
    main()
