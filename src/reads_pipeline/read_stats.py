from pathlib import Path
from subprocess import run
import logging
import os

from .paths import (
    get_reads_stats_fastqc_parent_dir,
    get_read_files_in_dir,
    get_log_path,
    FASTQC_BIN,
    FASTQ_EXT,
    get_reads_stats_parent_dir,
    get_raw_reads_parent_dir,
    get_clean_reads_parent_dir,
    get_raw_reads_stats_parent_dir,
    get_clean_reads_stats_parent_dir,
)

logger = logging.getLogger(__name__)


def run_cmd(cmd, project_dir):
    logging.basicConfig(filename=get_log_path(project_dir), level=logging.INFO)
    logging.info("Running cmd: " + " ".join(cmd))

    process = run(cmd, check=False, capture_output=True)
    if process.returncode:
        msg = f"There was a problem running: {cmd[0]}\n"
        msg += "stderr:\n" + process.stderr.decode()
        msg += "stdout:\n" + process.stdout.decode()
        logging.error(msg)
        raise RuntimeError("There was an error running the command: " + " ".join(cmd))

    return {"process": process}


def run_fastqc_for_file(reads_path, out_stats_dir, project_dir, re_run, threads: int):
    expected_zip_file = out_stats_dir / reads_path.name.replace(
        FASTQ_EXT, "_fastqc.zip"
    )
    expected_html_file = out_stats_dir / reads_path.name.replace(
        FASTQ_EXT, "_fastqc.html"
    )
    if re_run:
        if expected_zip_file.exists():
            os.remove(expected_zip_file)
        if expected_zip_file.exists():
            os.remove(expected_html_file)
    else:
        if expected_zip_file.exists():
            return

    cmd = [FASTQC_BIN, "-o", str(out_stats_dir)]
    if threads > 1:
        cmd.extend(["--threads", str(threads)])

    cmd.append(str(reads_path))

    run_cmd(cmd, project_dir)


def run_fastqc(project_dir: Path | str | None = None, re_run=False, threads=1):
    get_reads_stats_parent_dir(project_dir).mkdir(exist_ok=True)
    fastq_stats_dir = get_reads_stats_fastqc_parent_dir(project_dir)
    fastq_stats_dir.mkdir(exist_ok=True)

    for read_type in ["raw", "clean"]:
        if read_type == "raw":
            reads_parent_dir = get_raw_reads_parent_dir(project_dir)
            stats_dir = get_raw_reads_stats_parent_dir(project_dir)
        elif read_type == "clean":
            reads_parent_dir = get_clean_reads_parent_dir(project_dir)
            if not reads_parent_dir.exists():
                continue
            stats_dir = get_clean_reads_stats_parent_dir(project_dir)
        reads_dirs = [path for path in reads_parent_dir.iterdir() if path.is_dir()]

        stats_dir.mkdir(exist_ok=True)
        for reads_dir in reads_dirs:
            dir_name = reads_dir.name
            out_stats_dir = stats_dir / dir_name
            out_stats_dir.mkdir(exist_ok=True)
            for reads_path in get_read_files_in_dir(reads_dir):
                run_fastqc_for_file(
                    reads_path,
                    out_stats_dir,
                    project_dir,
                    re_run=re_run,
                    threads=threads,
                )


# num_reads
# %q30
# mean qual first 10 and last 10 nucleotides
# %gc
# composition acgt first 10 bases


def run_fastp():
    pass
