from pathlib import Path
import os
import zipfile
import statistics
import tempfile

import pandas

from .paths import (
    get_reads_stats_fastqc_parent_dir,
    get_read_files_in_dir,
    FASTQC_BIN,
    FASTQ_EXT,
    FASTQC_XLS_STATS_FNAME,
    get_reads_stats_parent_dir,
    get_raw_reads_parent_dir,
    get_clean_reads_parent_dir,
    get_raw_reads_fastqc_stats_parent_dir,
    get_clean_reads_fastqc_stats_parent_dir,
)
from .run_cmd import run_cmd
from .utils_file_system import move_files_and_dirs


def run_fastqc_for_file(
    reads_path,
    out_stats_dir,
    project_dir,
    re_run,
    threads: int,
    verbose: bool,
    dry_run: bool,
    n_files_to_process: int,
    n_file_processing: int,
):
    with tempfile.TemporaryDirectory(
        dir=out_stats_dir, prefix="fastqc_tmp_dir_"
    ) as fastq_result_tmp_dir:
        fastq_result_tmp_dir_path = Path(fastq_result_tmp_dir)
        expected_zip_file = out_stats_dir / reads_path.name.replace(
            FASTQ_EXT, "_fastqc.zip"
        )
        expected_html_file = out_stats_dir / reads_path.name.replace(
            FASTQ_EXT, "_fastqc.html"
        )
        if re_run:
            if expected_zip_file.exists():
                if verbose:
                    print(f"Removing previous fastqc analisys for {reads_path}")
                os.remove(expected_zip_file)
            if expected_zip_file.exists():
                os.remove(expected_html_file)
        else:
            if expected_zip_file.exists():
                if verbose:
                    print(f"Skipping fastqc analisys for {reads_path}")
                return {"run": False, "should_have_run": False}

        if verbose:
            print(
                f"Running fastqc {n_file_processing} out of {n_files_to_process} for {reads_path}"
            )

        cmd = [FASTQC_BIN, "-o", str(fastq_result_tmp_dir_path)]
        if threads > 1:
            cmd.extend(["--threads", str(threads)])

        cmd.append(str(reads_path))

        if verbose:
            print("cmd: ", " ".join(cmd))
        if not dry_run:
            run_cmd(cmd, project_dir)

        move_files_and_dirs(fastq_result_tmp_dir, out_stats_dir)
        return {"run": True, "should_have_run": True}


def _run_fastqc(
    project_dir: Path | str | None = None,
    re_run=False,
    threads=1,
    verbose=False,
    read_types=("raw", "clean"),
    dry_run=False,
    n_files_to_process=None,
):
    get_reads_stats_parent_dir(project_dir).mkdir(exist_ok=True)
    fastq_stats_dir = get_reads_stats_fastqc_parent_dir(project_dir)
    fastq_stats_dir.mkdir(exist_ok=True)

    n_files_processed = 1
    for read_type in read_types:
        if read_type == "raw":
            reads_parent_dir = get_raw_reads_parent_dir(project_dir)
            stats_dir = get_raw_reads_fastqc_stats_parent_dir(project_dir)
        elif read_type == "clean":
            reads_parent_dir = get_clean_reads_parent_dir(project_dir)
            if not reads_parent_dir.exists():
                continue
            stats_dir = get_clean_reads_fastqc_stats_parent_dir(project_dir)
        reads_dirs = [path for path in reads_parent_dir.iterdir() if path.is_dir()]

        stats_dir.mkdir(exist_ok=True)
        for reads_dir in reads_dirs:
            dir_name = reads_dir.name
            out_stats_dir = stats_dir / dir_name
            out_stats_dir.mkdir(exist_ok=True)
            for reads_path in get_read_files_in_dir(reads_dir):
                res = run_fastqc_for_file(
                    reads_path,
                    out_stats_dir,
                    project_dir,
                    re_run=re_run,
                    threads=threads,
                    verbose=verbose,
                    dry_run=dry_run,
                    n_files_to_process=n_files_to_process,
                    n_file_processing=n_files_processed,
                )
                n_files_processed += int(res["should_have_run"])
    return {"n_files_processed": n_files_processed - 1}


def run_fastqc(
    project_dir: Path | str | None = None,
    re_run=False,
    threads=1,
    verbose=False,
    read_types=("raw", "clean"),
):
    if verbose:
        res = _run_fastqc(
            project_dir=project_dir,
            re_run=re_run,
            threads=threads,
            verbose=False,
            read_types=read_types,
            dry_run=True,
        )
        n_files_to_process = res["n_files_processed"]
    else:
        n_files_to_process = None

    return _run_fastqc(
        project_dir=project_dir,
        re_run=re_run,
        threads=threads,
        verbose=verbose,
        read_types=read_types,
        dry_run=False,
        n_files_to_process=n_files_to_process,
    )


def _parse_fastqc_data(file_content: str, bioproject_dir: str):
    modules = file_content.split(">>END_MODULE")
    # for module in modules:
    #    print(module)
    basic_stats = dict([line.split("\t") for line in modules[0].splitlines()])
    result = {}
    result["dir"] = bioproject_dir
    result["file_name"] = basic_stats["Filename"]
    result["num_seqs"] = basic_stats["Total Sequences"]
    result[r"%GC"] = basic_stats[r"%GC"]

    mean_qual_first_9_nucleotides = statistics.mean(
        [float(line.split("\t")[1]) for line in modules[1].splitlines()[3:12]]
    )
    result["mean_qual_first_9_nucleotides"] = mean_qual_first_9_nucleotides
    return result


def _parse_fastqc_zip_file(zip_path: Path, bioproject_dir: str):
    with zipfile.ZipFile(zip_path, "r") as zip_file:
        zip_name = zip_path.name.replace(".zip", "")
        with zip_file.open(f"{zip_name}/fastqc_data.txt") as file:
            file_content = file.read().decode()
            return _parse_fastqc_data(file_content, bioproject_dir)


def collect_fastqc_stats(project_dir):
    result = {}
    for read_type in ["raw", "clean"]:
        if read_type == "raw":
            stats_dir = get_raw_reads_fastqc_stats_parent_dir(project_dir)
        elif read_type == "clean":
            stats_dir = get_clean_reads_fastqc_stats_parent_dir(project_dir)
        if not stats_dir.exists():
            continue
        fastq_stats_dirs = [path for path in stats_dir.iterdir() if path.is_dir()]

        stats = []
        for fastq_stats_dir in fastq_stats_dirs:
            bioproject_dir = fastq_stats_dir.name
            fastqc_zip_paths = [
                path
                for path in fastq_stats_dir.iterdir()
                if str(path).endswith("_fastqc.zip")
            ]
            for zip_path in fastqc_zip_paths:
                stats.append(_parse_fastqc_zip_file(zip_path, bioproject_dir))
        result[read_type] = pandas.DataFrame(stats)
        excel_path = stats_dir / FASTQC_XLS_STATS_FNAME
        result[read_type].to_excel(excel_path, index=False)
    return result
