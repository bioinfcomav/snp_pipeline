from pathlib import Path
import os

FASTQC_BIN = "fastqc"
FASTQ_EXT = ".fastq.gz"


def get_project_dir(project_dir: None | str | Path) -> Path:
    if project_dir is not None:
        return Path(project_dir)
    else:
        return Path(os.getcwd())


def get_reads_dir(project_dir, check_exists=True) -> Path:
    project_dir = get_project_dir(project_dir=project_dir)
    path = project_dir / "reads"
    if check_exists and not path.exists():
        raise ValueError(f"The reads parent dir does not exist: {path}")
    return path


def get_raw_reads_parent_dir(project_dir, check_exists=True) -> Path:
    reads_dir = get_reads_dir(project_dir, check_exists=check_exists)
    path = reads_dir / "raw"
    if check_exists and not path.exists():
        raise ValueError(f"The raw reads parent dir does not exist: {path}")
    return path


def get_clean_reads_parent_dir(project_dir) -> Path:
    reads_dir = get_reads_dir(project_dir)
    path = reads_dir / "clean"
    return path


def get_reads_stats_parent_dir(project_dir) -> Path:
    reads_dir = get_reads_dir(project_dir=project_dir)
    return reads_dir / "stats"


def get_reads_stats_fastqc_parent_dir(project_dir) -> Path:
    return get_reads_stats_parent_dir(project_dir) / "fastqc"


def get_raw_reads_stats_parent_dir(project_dir) -> Path:
    return get_reads_stats_fastqc_parent_dir(project_dir) / "raw"


def get_clean_reads_stats_parent_dir(project_dir) -> Path:
    return get_reads_stats_fastqc_parent_dir(project_dir) / "clean"


def is_read_file(path: Path):
    if path.is_dir():
        return False

    if str(path).endswith(FASTQ_EXT):
        return True
    return False


def get_read_files_in_dir(dir_path: Path) -> list[Path]:
    return [path for path in dir_path.iterdir() if is_read_file(path)]


def get_log_path(project_dir) -> Path:
    log_path = get_project_dir(project_dir) / "log.txt"
    return log_path
