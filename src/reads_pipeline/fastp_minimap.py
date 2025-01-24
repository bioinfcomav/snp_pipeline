from pathlib import Path
import logging

from .paths import (
    get_project_dir,
    get_log_path,
    get_raw_reads_parent_dir,
    get_reads_stats_fastp_parent_dir,
    get_paired_and_unpaired_read_files_in_dir,
    get_crams_dir,
    FASTP_BIN,
    MINIMAP2_BIN,
    SAMTOOLS_BIN,
)
from .run_cmd import run_bash_script

logger = logging.getLogger(__name__)

PIPE_TEMPLATE = """
#!/usr/bin/env bash
# Exit immediately if any command has a non-zero exit status.
set -e
# Enable pipefail to make the entire pipeline fail if any part of the pipe fails
#  If any command in a pipeline fails, that return code will be used as t   he return code of the whole pipeline
set -o pipefail

# fastp
{fastp_bin} {fastp_in1} {fastp_in2} --stdout -h {fastp_html_report_path} -j {fastp_json_report_path} --length_required {min_read_len} --overrepresentation_analysis --thread {fastp_num_threads} | \

# minimap2
{minimap2_bin} -t {minimap_num_threads} -a -x sr {minimap_index} - | \

# This adds mate cigar (MC) and mate score tags (ms) which will be used later by samtools markdup proper
{samtools_bin} fixmate -u -m - - | \

# When estimating the total number of concurrent threads to allocate,
# consider that the sort step is a crunch point that separates the steps before it from the step afterwards.
# The mapping, fixmate and partial sort to temporary file steps will operate in parallel.
# Once complete, the sort merge (from temporary files) and markdup steps will then run in parallel.
# Sort is highly parallel so the -@8 option here enables to use of 8 additional CPU threads.
# It can also be sped up by providing it with more memory, but note the memory option (-m) is per-thread.
# The -l 1 indicates level 1 compression again. We could also specify -O bam,level=1 as used above.
{samtools_bin} sort -u -@{sort_num_threads} -T {tmp_dir} - | \
# The main core of marking duplicates may now be ran on the position-sorted file
{samtools_bin} markdup -@{duplicates_num_threads} --reference {genome_fasta} - {cram_path}

# create index
samtools index {cram_path}

# stats
samtools stats {cram_path} > {cram_stats_path}
"""


def _run_fastp_minimap_for_pair(
    pair: tuple[Path],
    project_dir: Path,
    stats_dir: Path,
    crams_dir: Path,
    min_read_len: int,
    fastp_num_threads: int,
    minimap_index: Path,
    minimap_num_threads: int,
    sort_num_threads: int,
    tmp_dir: Path,
    duplicates_num_threads: int,
    genome_fasta: Path,
):
    if len(pair) == 2:
        fastp_in1 = f"--in1 {pair[0]}"
        fastp_in2 = f"--in2 {pair[1]}"
    elif len(pair) == 1:
        fastp_in1 = f"--in1 {pair[0]}"
        fastp_in2 = ""

    html_report_path = stats_dir / pair[0].with_suffix(".html").name
    json_report_path = stats_dir / pair[0].with_suffix(".json").name

    cram_path = crams_dir / pair[0].with_suffix(".cram").name

    cram_stats_path = crams_dir / pair[0].with_suffix(".cram.stats").name

    logging.info(
        "Running fastp-minimap pipeline for files: " + " ".join(map(str, pair))
    )

    script = PIPE_TEMPLATE.format(
        fastp_bin=FASTP_BIN,
        fastp_in1=fastp_in1,
        fastp_in2=fastp_in2,
        fastp_html_report_path=html_report_path,
        fastp_json_report_path=json_report_path,
        min_read_len=min_read_len,
        fastp_num_threads=fastp_num_threads,
        minimap2_bin=MINIMAP2_BIN,
        minimap_index=minimap_index,
        minimap_num_threads=minimap_num_threads,
        samtools_bin=SAMTOOLS_BIN,
        sort_num_threads=sort_num_threads,
        tmp_dir=tmp_dir,
        duplicates_num_threads=duplicates_num_threads,
        genome_fasta=genome_fasta,
        cram_path=cram_path,
        cram_stats_path=cram_stats_path,
    )
    run_bash_script(script, project_dir=project_dir)


def run_fastp_minimap(
    project_dir: Path,
    minimap_index: Path,
    genome_fasta: Path,
    min_read_len=30,
    fastp_num_threads=3,
    minimap_num_threads=3,
    sort_num_threads=8,
    duplicates_num_threads=8,
):
    project_dir = get_project_dir(project_dir)
    tmp_dir = project_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        force=True,
    )
    logging.debug("hola")
    print(get_log_path(project_dir))

    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)
    stats_parent_dir = get_reads_stats_fastp_parent_dir(project_dir)
    stats_parent_dir.mkdir(exist_ok=True, parents=True)
    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]
    crams_parent_dir = get_crams_dir(project_dir)
    crams_parent_dir.mkdir(exist_ok=True, parents=True)

    for raw_reads_dir in raw_reads_dirs:
        dir_name = raw_reads_dir.name
        stats_dir = stats_parent_dir / dir_name
        stats_dir.mkdir(exist_ok=True)
        crams_dir = crams_parent_dir / dir_name
        crams_dir.mkdir(exist_ok=True)

        for pair in get_paired_and_unpaired_read_files_in_dir(raw_reads_dir):
            _run_fastp_minimap_for_pair(
                pair,
                project_dir=project_dir,
                stats_dir=stats_dir,
                crams_dir=crams_dir,
                min_read_len=min_read_len,
                fastp_num_threads=fastp_num_threads,
                minimap_index=minimap_index,
                minimap_num_threads=minimap_num_threads,
                sort_num_threads=sort_num_threads,
                tmp_dir=tmp_dir,
                duplicates_num_threads=duplicates_num_threads,
                genome_fasta=genome_fasta,
            )
    input("hola")
