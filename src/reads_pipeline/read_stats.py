from pathlib import Path
from subprocess import run
import logging

from .paths import (
    get_raw_reads_parent_dir,
    get_reads_stats_fastqc_parent_dir,
    get_read_files_in_dir,
    get_log_path,
    FASTQC_BIN,
    get_reads_stats_parent_dir,
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


def run_fastqc_for_file(reads_path, out_stats_dir, project_dir):
    cmd = [FASTQC_BIN, "-o", str(out_stats_dir), str(reads_path)]
    run_cmd(cmd, project_dir)


def run_fastqc(project_dir: Path | str | None = None, re_run=False):
    raw_reads_parent_dir = get_raw_reads_parent_dir(project_dir)

    raw_reads_dirs = [path for path in raw_reads_parent_dir.iterdir() if path.is_dir()]

    get_reads_stats_parent_dir(project_dir).mkdir(exist_ok=True)
    fastq_stats_dir = get_reads_stats_fastqc_parent_dir(project_dir)
    fastq_stats_dir.mkdir(exist_ok=True)

    raw_stats_dir = fastq_stats_dir / "raw"
    raw_stats_dir.mkdir(exist_ok=True)
    for raw_reads_dir in raw_reads_dirs:
        dir_name = raw_reads_dir.name
        out_stats_dir = raw_stats_dir / dir_name
        out_stats_dir.mkdir(exist_ok=True)
        for reads_path in get_read_files_in_dir(raw_reads_dir):
            run_fastqc_for_file(reads_path, out_stats_dir, project_dir)


def run_fastp():
    pass
