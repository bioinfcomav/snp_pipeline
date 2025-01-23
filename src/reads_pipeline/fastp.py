from pathlib import Path


from .paths import (
    get_paired_and_unpaired_read_files_in_dir,
    get_raw_reads_parent_dir,
    get_clean_reads_parent_dir,
    get_clean_reads_stats_parent_dir,
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
    stats_parent_dir = get_clean_reads_stats_parent_dir(project_dir)

    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]
    for raw_reads_dir in raw_reads_dirs:
        dir_name = raw_reads_dir.name
        clean_reads_dir = clean_reads_parent_dir / dir_name
        clean_reads_dir.mkdir(exist_ok=True)
        stats_dir = stats_parent_dir / dir_name
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
