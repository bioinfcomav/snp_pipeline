from pathlib import Path
import json

import pandas

from .paths import (
    get_paired_and_unpaired_read_files_in_dir,
    get_raw_reads_parent_dir,
    get_clean_reads_parent_dir,
    get_reads_stats_fastp_parent_dir,
    FASTP_BIN,
    remove_file,
)
from .run_cmd import run_cmd


def _run_fastp_for_pair(
    pair: tuple[Path],
    clean_reads_dir: Path,
    stats_dir: Path,
    min_len: int,
    deduplicate: bool,
    threads: int,
    project_dir: Path,
    re_run: bool,
):
    clean_path1 = clean_reads_dir / pair[0].name

    if re_run:
        remove_file(clean_path1, not_exist_ok=True)

    cmd = [FASTP_BIN]
    cmd.extend(["-i", str(pair[0]), "-o", str(clean_path1)])

    if len(pair) == 2:
        clean_path2 = clean_reads_dir / pair[1].name
        if re_run:
            remove_file(clean_path2, not_exist_ok=True)
        cmd.extend(["-I", str(pair[1]), "-O", str(clean_path2)])

    html_report_path = stats_dir / pair[1].with_suffix(".html").name
    json_report_path = stats_dir / pair[1].with_suffix(".json").name
    if re_run:
        remove_file(html_report_path, not_exist_ok=True)
        remove_file(json_report_path, not_exist_ok=True)

    cmd.extend(["-h", str(html_report_path), "-j", str(json_report_path)])

    cmd.extend(["--length_required", str(min_len)])

    cmd.extend(["--cut_front", "--cut_tail"])

    cmd.append("--overrepresentation_analysis")

    if deduplicate:
        cmd.append("--dedup")
    else:
        cmd.append("--dedup")

    cmd.extend(["--thread", str(threads)])
    cmd.append("--dont_overwrite")

    run_cmd(cmd, project_dir=project_dir)


def run_fastp(project_dir, min_len=30, deduplicate=True, threads=3, re_run=False):
    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)
    clean_reads_parent_dir = get_clean_reads_parent_dir(project_dir)
    clean_reads_parent_dir.mkdir(exist_ok=True)
    stats_parent_dir = get_reads_stats_fastp_parent_dir(project_dir)
    stats_parent_dir.mkdir(exist_ok=True, parents=True)

    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]
    for raw_reads_dir in raw_reads_dirs:
        dir_name = raw_reads_dir.name
        clean_reads_dir = clean_reads_parent_dir / dir_name
        clean_reads_dir.mkdir(exist_ok=True)
        stats_dir = stats_parent_dir / dir_name
        stats_dir.mkdir(exist_ok=True)
        for pair in get_paired_and_unpaired_read_files_in_dir(raw_reads_dir):
            _run_fastp_for_pair(
                pair,
                clean_reads_dir,
                stats_dir,
                min_len=min_len,
                deduplicate=deduplicate,
                threads=threads,
                project_dir=project_dir,
                re_run=re_run,
            )


def _parse_fastp_json(path):
    with path.open("rt") as fhand:
        report = json.load(fhand)

    result = {}
    summary = report["summary"]["before_filtering"]
    result["num_raw_reads"] = summary["total_reads"]
    result["num_raw_bases"] = summary["total_bases"]
    result["raw_q20"] = summary["q20_rate"]
    result["raw_q30"] = summary["q30_rate"]
    result["raw_read1_mean_length"] = summary["read1_mean_length"]
    result["raw_read2_mean_length"] = summary["read2_mean_length"]

    if "after_filtering" in report["summary"]:
        summary = report["summary"]["after_filtering"]
        result["num_clean_reads"] = summary["total_reads"]
        result["num_clean_bases"] = summary["total_bases"]
        result["clean_q20"] = summary["q20_rate"]
        result["clean_q30"] = summary["q30_rate"]
        result["clean_read1_mean_length"] = summary["read1_mean_length"]
        result["clean_read2_mean_length"] = summary["read2_mean_length"]

    if "filtering_results" in report:
        summary = report["filtering_results"]
        result["num_low_qual_reads"] = summary["filtering_results"]["low_quality_reads"]
        result["num_too_many_N_reads"] = summary["filtering_results"][
            "too_many_N_reads"
        ]
        result["too_short_reads"] = summary["filtering_results"]["too_short_reads"]
        result["too_long_reads"] = summary["filtering_results"]["too_long_reads"]
    if "duplication" in report:
        result["duplication_rate"] = report["duplication"]["rate"]
    if "insert_size" in report:
        result["insert_size"] = report["insert_size"]["peak"]
    if "adapter_cutting" in report:
        result["adapter_trimmed_reads"] = report["adapter_cutting"][
            "adapter_trimmed_reads"
        ]

    if "read1_before_filtering" in report:
        result["raw_read1_q30"] = (
            report["read1_before_filtering"]["q30_bases"]
            / report["read1_before_filtering"]["total_bases"]
        )
    if "read2_before_filtering" in report:
        result["raw_read2_q30"] = (
            report["read2_before_filtering"]["q30_bases"]
            / report["read2_before_filtering"]["total_bases"]
        )
    if "read1_after_filtering" in report:
        result["clean_read1_q30"] = (
            report["read1_after_filtering"]["q30_bases"]
            / report["read1_after_filtering"]["total_bases"]
        )
    if "read2_after_filtering" in report:
        result["clean_read2_q30"] = (
            report["read2_after_filtering"]["q30_bases"]
            / report["read2_after_filtering"]["total_bases"]
        )

    return result


def collect_fastp_stats(project_dir):
    stats_parent_dir = get_reads_stats_fastp_parent_dir(project_dir)
    stats_dirs = [path for path in stats_parent_dir.iterdir() if path.is_dir()]
    reports = []
    for stats_dir in stats_dirs:
        json_reports = [path for path in stats_dir.iterdir() if path.suffix == ".json"]
        for path in json_reports:
            read_report = {}
            read_report["dir"] = stats_dir.name
            read_report["file_name"] = path.name
            read_report |= _parse_fastp_json(path)
            reports.append(read_report)
    return pandas.DataFrame(reports)
