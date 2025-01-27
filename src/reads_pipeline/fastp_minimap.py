from pathlib import Path
import logging
import os

import numpy
import pandas
import matplotlib.pyplot as plt


from .paths import (
    get_project_dir,
    get_log_path,
    get_raw_reads_parent_dir,
    get_reads_stats_fastp_parent_dir,
    get_paired_and_unpaired_read_files_in_dir,
    get_crams_dir,
    get_crams_stats_dir,
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

    cram_path = crams_dir / (pair[0].name.removesuffix(".fastq.gz") + ".cram")

    cram_stats_path = crams_dir / (
        pair[0].name.removesuffix(".fastq.gz") + ".cram.stats"
    )

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


def _parse_cram_stats(cram_stats_path):
    with cram_stats_path.open("rt") as fhand:
        lines = list(fhand)

    file_line = lines[2]
    if not file_line.startswith("# The command line was"):
        raise RuntimeError(
            f"Error reading the command file line in cram stats: {cram_stats_path}"
        )
    file_name = file_line.split(os.sep)[-1]
    result = {"file_name": file_name}

    summary = dict(
        [
            line.strip().split("\t", 1)[1].replace(":", "").split("\t")[:2]
            for line in lines
            if line.startswith("SN")
        ]
    )
    result["total sequences"] = int(summary["raw total sequences"])
    result["filtered sequences"] = int(summary["filtered sequences"])
    result["sequences"] = int(summary["sequences"])
    result["1st fragments"] = int(summary["1st fragments"])
    result["reads mapped"] = int(summary["reads mapped"])
    result[r"% reads mapped"] = result["reads mapped"] / result["sequences"] * 100.0
    result["reads mapped and paired"] = int(summary["reads mapped and paired"])
    result[r"% reads mapped and paired"] = (
        result["reads mapped and paired"] / result["sequences"] * 100.0
    )
    result["reads unmapped"] = int(summary["reads unmapped"])
    result[r"% reads unmapped"] = result["reads unmapped"] / result["sequences"] * 100.0
    result["reads duplicated"] = int(summary["reads duplicated"])
    result[r"% reads duplicated"] = (
        result["reads duplicated"] / result["sequences"] * 100.0
    )

    result["reads properly paired"] = int(summary["reads properly paired"])
    result["reads paired"] = int(summary["reads paired"])
    result[r"% reads properly paired"] = (
        result["reads properly paired"] / result["reads paired"] * 100.0
    )
    result["inward oriented pairs"] = int(summary["inward oriented pairs"])
    result["outward oriented pairs"] = int(summary["outward oriented pairs"])
    result["insert size average"] = float(summary["insert size average"])

    result["reads average length"] = float(summary["average length"])
    result["reads average quality"] = float(summary["average quality"])
    return result


def collect_cram_stats(project_dir):
    project_dir = get_project_dir(project_dir)
    crams_parent_dir = get_crams_dir(project_dir)
    crams_dirs = [path for path in crams_parent_dir.iterdir() if path.is_dir()]

    results = []
    for crams_dir in crams_dirs:
        for cram_stats_path in filter(
            lambda x: x.suffixes[-2:] == [".cram", ".stats"], crams_dir.iterdir()
        ):
            results.append({"dir": crams_dir.name} | _parse_cram_stats(cram_stats_path))
    results = pandas.DataFrame(results)

    stats_dir = get_crams_stats_dir(project_dir)
    stats_dir.mkdir(exist_ok=True)
    stats_path = stats_dir / "cram_stats.xlsx"
    results.to_excel(stats_path, index=False)


def _get_cram_stat_paths(project_dir):
    project_dir = get_project_dir(project_dir)
    crams_parent_dir = get_crams_dir(project_dir)
    crams_dirs = [path for path in crams_parent_dir.iterdir() if path.is_dir()]

    for crams_dir in crams_dirs:
        for cram_stats_path in filter(
            lambda x: x.suffixes[-2:] == [".cram", ".stats"], crams_dir.iterdir()
        ):
            yield cram_stats_path


def _parse_cram_mapq(cram_stats_path):
    with cram_stats_path.open("rt") as fhand:
        mapq_num_reads = dict(
            [
                map(int, line.removeprefix("MAPQ\t").strip().split("\t"))
                for line in fhand
                if line.startswith("MAPQ\t")
            ]
        )
        mapqs = numpy.arange(0, max(60, max(mapq_num_reads.keys()) + 1), dtype=int)
        num_reads = [mapq_num_reads.get(mapq, 0) for mapq in mapqs]
    return mapqs, num_reads


def plot_mapq_distributions(project_dir):
    parent_stats_dir = get_crams_stats_dir(project_dir)
    parent_stats_dir.mkdir(exist_ok=True)

    for cram_stats_path in _get_cram_stat_paths(project_dir):
        dir_name = cram_stats_path.parent.name
        stats_dir = parent_stats_dir / dir_name
        stats_dir.mkdir(exist_ok=True)
        fname = (
            "mapq_distrib."
            + str(cram_stats_path.name).removesuffix(".cram.stats")
            + ".svg"
        )

        mapqs, num_reads = _parse_cram_mapq(cram_stats_path)

        fig, axes = plt.subplots()
        axes.bar(mapqs, num_reads)
        axes.set_xlim((0, max(axes.get_xlim()[1], 60)))
        axes.set_xlabel("MAPQ")
        axes.set_ylabel("Num. reads (non-duplicated)")

        plot_path = stats_dir / fname
        fig.savefig(plot_path)
        plt.close(fig)
