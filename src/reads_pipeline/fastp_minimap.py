from pathlib import Path
import logging
import os
import shutil
import tempfile

import numpy
import pandas
import matplotlib.pyplot as plt


from .paths import (
    get_project_dir,
    get_read_group_info_xls,
    get_log_path,
    get_raw_reads_parent_dir,
    get_reads_stats_fastp_parent_dir,
    get_paired_and_unpaired_read_files_in_dir,
    get_crams_dir,
    get_crams_stats_dir,
    FASTP_BIN,
    MINIMAP2_BIN,
    SAMTOOLS_BIN,
    TRIM_QUALS_BIN,
    FILE_BIN,
    MD5BIN,
    remove_file,
)
from .run_cmd import run_bash_script, run_cmd
from .read_group import get_read_group_info, create_minimap_rg_str
from .utils_file_system import move_files_and_dirs

logger = logging.getLogger(__name__)

PIPE_TEMPLATE = """
#!/usr/bin/env bash
# Exit immediately if any command has a non-zero exit status.
set -e
# Enable pipefail to make the entire pipeline fail if any part of the pipe fails
#  If any command in a pipeline fails, that return code will be used as t   he return code of the whole pipeline
set -o pipefail

export REF_PATH={ref_path_dir}

# fastp
{fastp_bin} {fastp_in1} {fastp_in2} --stdout -h {fastp_html_report_path} -j {fastp_json_report_path} --length_required {min_read_len} --overrepresentation_analysis {fastp_gobal_trim} --thread {fastp_num_threads} | \\

# minimap2
{minimap2_bin} -R {rg_str} -t {minimap_num_threads} -a -x sr {minimap_index} - | \\

# This adds mate cigar (MC) and mate score tags (ms) which will be used later by samtools markdup proper
{samtools_bin} fixmate -u -m - - | \\

{deduplicate_and_sort_lines}

# calmd
# Calmd works best on position-sorted input files, as with these it can stream through the reference sequence
# and so doesn't have to store much reference data at any one time
# -A when used jointly with -r this option overwrites the original base quality
{samtools_bin} calmd -Ar -@{calmd_num_threads} - {genome_fasta} | \\

# trim_quals it reduces the qualities from the read edges
# In the past there has been problems running trim_quals before calmd, so run it after calmd
{trim_quals_line}

samtools view --reference {genome_fasta} -o {cram_path} - 

# create index
samtools index {cram_path}

# stats
samtools stats {cram_path} > {cram_stats_path}
"""

SORT_AND_DEDUPLIATE_LINES = """# When estimating the total number of concurrent threads to allocate,
# consider that the sort step is a crunch point that separates the steps before it from the step afterwards.
# The mapping, fixmate and partial sort to temporary file steps will operate in parallel.
# Once complete, the sort merge (from temporary files) and markdup steps will then run in parallel.
# Sort is highly parallel so the -@8 option here enables to use of 8 additional CPU threads.
# It can also be sped up by providing it with more memory, but note the memory option (-m) is per-thread.
# The -l 1 indicates level 1 compression again. We could also specify -O bam,level=1 as used above.
{samtools_bin} sort -u -@{sort_num_threads} -T {tmp_dir} - | \\

# The main core of marking duplicates may now be ran on the position-sorted file
{samtools_bin} markdup -@{duplicates_num_threads} --reference {genome_fasta} - - | \\"""

SORT_LINES = """# When estimating the total number of concurrent threads to allocate,
# consider that the sort step is a crunch point that separates the steps before it from the step afterwards.
# The mapping, fixmate and partial sort to temporary file steps will operate in parallel.
# Once complete, the sort merge (from temporary files) and markdup steps will then run in parallel.
# Sort is highly parallel so the -@8 option here enables to use of 8 additional CPU threads.
# It can also be sped up by providing it with more memory, but note the memory option (-m) is per-thread.
# The -l 1 indicates level 1 compression again. We could also specify -O bam,level=1 as used above.
{samtools_bin} sort -u -@{sort_num_threads} -T {tmp_dir} --reference {genome_fasta} - | \\"""

TRIM_QUALS_LINE = (
    "trim_quals --num-bases {num_bases} --qual-reduction {qual_reduction} - - | \\"
)


def _run_fastp_minimap_for_pair(
    pair: tuple[Path],
    project_dir: Path,
    stats_dir: Path,
    crams_dir: Path,
    min_read_len: int,
    fastp_num_threads: int,
    fastp_trim_front1: int,
    fastp_trim_tail1: int,
    fastp_trim_front2: int,
    fastp_trim_tail2: int,
    minimap_index: Path,
    minimap_num_threads: int,
    sort_num_threads: int,
    calmd_num_threads: int,
    duplicates_num_threads: int,
    genome_fasta: Path,
    deduplicate: bool,
    re_run: bool,
    read_groups_info: dict,
    trim_quals_num_bases: int,
    trim_quals_qual_reduction: int,
    genome_md5: str,
    verbose: bool,
):
    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        force=True,
    )
    read_id = pair[0].name.split(".")[0]
    excel_file = get_read_group_info_xls(project_dir)
    if read_id not in read_groups_info:
        msg = f"The read group {read_id} does not have read group info in the excel file: {excel_file}"
        logging.error(msg)
        raise ValueError(msg)
    rg_str = create_minimap_rg_str(
        read_id, read_groups_info[read_id], project_dir, excel_file
    )

    if len(pair) == 2:
        fastp_in1 = f"--in1 {pair[0]}"
        fastp_in2 = f"--in2 {pair[1]}"
    elif len(pair) == 1:
        fastp_in1 = f"--in1 {pair[0]}"
        fastp_in2 = ""
    fastp_gobal_trim = ""
    if fastp_trim_front1:
        fastp_gobal_trim += f"--trim_front1={fastp_trim_front1} "
    if fastp_trim_front2:
        fastp_gobal_trim += f"--trim_front2={fastp_trim_front2} "
    if fastp_trim_tail1:
        fastp_gobal_trim += f"--trim_tail1={fastp_trim_tail1} "
    if fastp_trim_tail2:
        fastp_gobal_trim += f"--trim_tail2={fastp_trim_tail2} "

    with tempfile.TemporaryDirectory(
        prefix="cram_tmp_dir_", dir=crams_dir
    ) as crams_tmp_dir:
        crams_tmp_dir_path = Path(crams_tmp_dir)

        ref_path_dir = crams_tmp_dir_path / "ref_path"
        ref_path_dir.mkdir(exist_ok=True)
        genome_in_ref_cache = ref_path_dir / genome_md5
        if not genome_in_ref_cache.exists():
            os.symlink(genome_fasta, genome_in_ref_cache)

        tmp_dir = tempfile.TemporaryDirectory(prefix="tmp_dir_", dir=crams_tmp_dir)

        html_report_path = stats_dir / pair[0].with_suffix(".html").name
        json_report_path = stats_dir / pair[0].with_suffix(".json").name
        html_report_tmp_path = crams_tmp_dir_path / pair[0].with_suffix(".html").name
        json_report_tmp_path = crams_tmp_dir_path / pair[0].with_suffix(".json").name

        cram_path = crams_dir / (pair[0].name.removesuffix(".fastq.gz") + ".cram")
        cram_tmp_path = crams_tmp_dir_path / (
            pair[0].name.removesuffix(".fastq.gz") + ".cram"
        )

        cram_stats_path = crams_dir / (
            pair[0].name.removesuffix(".fastq.gz") + ".cram.stats"
        )
        cram_stats_tmp_path = crams_tmp_dir_path / (
            pair[0].name.removesuffix(".fastq.gz") + ".cram.stats"
        )

        if re_run:
            if verbose:
                print(
                    f"Removing cram and cram stats files to rerun analysis: {cram_path}, {cram_stats_path}"
                )
            remove_file(cram_path, not_exist_ok=True)
            remove_file(cram_stats_path, not_exist_ok=True)
        else:
            if cram_path.exists() and cram_stats_path.exists():
                print(f"Skipping analysis for existing file: {cram_path}")
                return
            else:
                remove_file(cram_path, not_exist_ok=True)
                remove_file(cram_stats_path, not_exist_ok=True)
                logging.warning(
                    f"One mapped file or stat was missing the existing file was removed and the analysis was re done: {cram_path} or {cram_stats_path}"
                )

        logging.info(
            "Running fastp-minimap pipeline for files: " + " ".join(map(str, pair))
        )
        if verbose:
            print(
                "Running fastp-minimap pipeline for files: " + " ".join(map(str, pair))
            )

        if deduplicate:
            deduplicate_line = SORT_AND_DEDUPLIATE_LINES.format(
                samtools_bin=SAMTOOLS_BIN,
                duplicates_num_threads=duplicates_num_threads,
                genome_fasta=genome_fasta,
                cram_path=cram_path,
                tmp_dir=tmp_dir.name,
                sort_num_threads=sort_num_threads,
            )
        else:
            deduplicate_line = SORT_LINES.format(
                samtools_bin=SAMTOOLS_BIN,
                genome_fasta=genome_fasta,
                cram_path=cram_path,
                tmp_dir=tmp_dir.name,
                sort_num_threads=sort_num_threads,
            )

        if trim_quals_num_bases > 0:
            trim_quals_line = TRIM_QUALS_LINE.format(
                num_bases=trim_quals_num_bases, qual_reduction=trim_quals_qual_reduction
            )
        else:
            trim_quals_line = ""

        script = PIPE_TEMPLATE.format(
            fastp_bin=FASTP_BIN,
            fastp_in1=fastp_in1,
            fastp_in2=fastp_in2,
            fastp_html_report_path=html_report_tmp_path,
            fastp_json_report_path=json_report_tmp_path,
            min_read_len=min_read_len,
            fastp_num_threads=fastp_num_threads,
            fastp_gobal_trim=fastp_gobal_trim,
            minimap2_bin=MINIMAP2_BIN,
            minimap_index=minimap_index,
            minimap_num_threads=minimap_num_threads,
            rg_str=rg_str,
            samtools_bin=SAMTOOLS_BIN,
            trim_quals_bin=TRIM_QUALS_BIN,
            sort_num_threads=sort_num_threads,
            tmp_dir=tmp_dir,
            genome_fasta=genome_fasta,
            cram_path=cram_tmp_path,
            cram_stats_path=cram_stats_tmp_path,
            deduplicate_and_sort_lines=deduplicate_line,
            calmd_num_threads=calmd_num_threads,
            trim_quals_line=trim_quals_line,
            ref_path_dir=ref_path_dir,
        )
        try:
            run_bash_script(script, project_dir=project_dir)
        except RuntimeError:
            remove_file(cram_path, not_exist_ok=True)
            remove_file(cram_stats_path, not_exist_ok=True)
            raise

        tmp_dir.cleanup()
        tmp_dir_path = Path(tmp_dir.name)
        if tmp_dir_path.exists():
            tmp_dir_path.rmdir()

        shutil.move(html_report_tmp_path, html_report_path)
        shutil.move(json_report_tmp_path, json_report_path)
        move_files_and_dirs(crams_tmp_dir, crams_dir)
    return {"cram_path": cram_path}


def _get_text_file_md5(genome_fasta, project_dir, uncompress_if_gzipped=False):
    if uncompress_if_gzipped:
        cmd = [FILE_BIN, "--mime-type", str(genome_fasta)]
        process = run_cmd(cmd, project_dir=project_dir)["process"]
        mimetype = process.stdout.decode().strip()
        is_gzipped = "gzip" in mimetype
    else:
        is_gzipped = False

    if is_gzipped:
        script_content = f"zcat {genome_fasta} | {MD5BIN}"
    else:
        script_content = f"{MD5BIN} {genome_fasta}"

    process = run_bash_script(script_content, project_dir=project_dir)["process"]
    genome_md5 = process.stdout.decode().split()[0]
    return genome_md5


def run_fastp_minimap(
    project_dir: Path,
    minimap_index: Path,
    genome_fasta: Path,
    deduplicate: bool,
    min_read_len=30,
    fastp_num_threads=3,
    fastp_trim_front1=0,
    fastp_trim_tail1=0,
    fastp_trim_front2=0,
    fastp_trim_tail2=0,
    minimap_num_threads=3,
    sort_num_threads=8,
    duplicates_num_threads=8,
    calmd_num_threads=2,
    trim_quals_num_bases=4,
    trim_quals_qual_reduction=20,
    re_run=False,
    verbose=True,
):
    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome fasta file not found: {genome_fasta}")
    genome_fasta = genome_fasta.resolve()
    if not minimap_index.exists():
        raise FileNotFoundError(f"Minimap index file not found: {minimap_index}")
    minimap_index = minimap_index.resolve()

    project_dir = get_project_dir(project_dir)

    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        force=True,
    )

    read_groups_info = get_read_group_info(project_dir)

    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)
    stats_parent_dir = get_reads_stats_fastp_parent_dir(project_dir)
    stats_parent_dir.mkdir(exist_ok=True, parents=True)
    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]
    crams_parent_dir = get_crams_dir(project_dir)
    crams_parent_dir.mkdir(exist_ok=True, parents=True)

    genome_md5 = _get_text_file_md5(
        genome_fasta, project_dir, uncompress_if_gzipped=True
    )

    cram_paths = []
    for raw_reads_dir in raw_reads_dirs:
        dir_name = raw_reads_dir.name
        stats_dir = stats_parent_dir / dir_name
        stats_dir.mkdir(exist_ok=True)
        crams_dir = crams_parent_dir / dir_name
        crams_dir.mkdir(exist_ok=True)

        for idx, pair in enumerate(
            get_paired_and_unpaired_read_files_in_dir(raw_reads_dir)
        ):
            if verbose:
                pair_str = ", ".join(map(str, pair))
                print(f"Cleaning and mapping read pair num {idx}: {pair_str}")
            res = _run_fastp_minimap_for_pair(
                pair,
                project_dir=project_dir,
                stats_dir=stats_dir,
                crams_dir=crams_dir,
                min_read_len=min_read_len,
                fastp_num_threads=fastp_num_threads,
                fastp_trim_front1=fastp_trim_front1,
                fastp_trim_tail1=fastp_trim_tail1,
                fastp_trim_front2=fastp_trim_front2,
                fastp_trim_tail2=fastp_trim_tail2,
                minimap_index=minimap_index,
                minimap_num_threads=minimap_num_threads,
                sort_num_threads=sort_num_threads,
                calmd_num_threads=calmd_num_threads,
                duplicates_num_threads=duplicates_num_threads,
                genome_fasta=genome_fasta,
                deduplicate=deduplicate,
                re_run=re_run,
                read_groups_info=read_groups_info,
                trim_quals_num_bases=trim_quals_num_bases,
                trim_quals_qual_reduction=trim_quals_qual_reduction,
                genome_md5=genome_md5,
                verbose=verbose,
            )
            if res:
                cram_paths.append(res["cram_path"])
    return {"cram_paths": cram_paths}


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

    mapqs, num_reads = _parse_cram_mapq(cram_stats_path)
    for mapq in (10, 20, 30):
        idx10 = mapqs >= 10
        result[f"%reads above {mapq} mapq"] = (
            numpy.sum(num_reads[idx10]) / numpy.sum(num_reads)
        ) * 100.0
    return result


def collect_cram_stats(project_dir):
    project_dir = get_project_dir(project_dir)

    results = []
    for cram_stats_path in _get_cram_stat_paths(project_dir):
        crams_dir_name = cram_stats_path.parent.name
        results.append({"dir": crams_dir_name} | _parse_cram_stats(cram_stats_path))
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
    return mapqs, numpy.array(num_reads)


def _process_cov(cov_range: str):
    if cov_range[0] != "[" or cov_range[-1] != "]":
        raise ValueError(
            f"The coverage range is malformed, expected: [int-int]: {cov_range}"
        )
    range = tuple(map(int, cov_range[1:-1].split("-")))
    if range[0] == range[1]:
        range = range[0]
    return range


def _parse_cram_cov(cram_stats_path):
    with cram_stats_path.open("rt") as fhand:
        cov_num_readss = [
            line.strip().removeprefix("COV\t").split("\t")
            for line in fhand
            if line.startswith("COV")
        ]
        cov_num_readss = [
            (_process_cov(cov_num_reads[0]), int(cov_num_reads[2]))
            for cov_num_reads in cov_num_readss
        ]
        covs, num_reads = zip(*cov_num_readss)
        all_ranges_are_int = all([isinstance(cov, int) for cov in covs])

        if all_ranges_are_int:
            covs_num_readss = dict(zip(covs, num_reads))
            new_covs = list(range(min(covs), max(covs)))
            new_num_reads = []
            for cov in new_covs:
                new_num_reads.append(covs_num_readss.get(cov, 0))
            covs, num_reads = new_covs, new_num_reads

        return covs, num_reads


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


def plot_coverage_distributions(project_dir):
    parent_stats_dir = get_crams_stats_dir(project_dir)
    parent_stats_dir.mkdir(exist_ok=True)

    for cram_stats_path in _get_cram_stat_paths(project_dir):
        dir_name = cram_stats_path.parent.name
        stats_dir = parent_stats_dir / dir_name
        stats_dir.mkdir(exist_ok=True)
        fname = (
            "cov_distrib."
            + str(cram_stats_path.name).removesuffix(".cram.stats")
            + ".svg"
        )
        covs, num_reads = _parse_cram_cov(cram_stats_path)
        if not all([isinstance(cov, int) for cov in covs]):
            covs = list(map(str, covs))
        fig, axes = plt.subplots()
        axes.bar(covs, num_reads)
        axes.set_xlabel("Coverage")
        axes.set_ylabel("Num. reads (non-duplicated)")

        plot_path = stats_dir / fname
        fig.savefig(plot_path)
        plt.close(fig)
