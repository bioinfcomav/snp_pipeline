from pathlib import Path
import os
from collections import defaultdict
import re

FASTQC_BIN = "fastqc"
FASTQ_EXT = ".fastq.gz"
FASTP_BIN = "fastp"
MINIMAP2_BIN = "minimap2"
SAMTOOLS_BIN = "samtools"
TABIX_BIN = "tabix"
BCFTOOLS_BIN = "bcftools"
TRIM_QUALS_BIN = "trim_quals"
SEQ_STATS_BIN = "seq_stats"
GATK_PYTHON_BIN = Path("/opt/gatk/gatk")
MD5BIN = "md5sum"
FILE_BIN = "file"
FASTQC_XLS_STATS_FNAME = "fastqc_stats.xls"


def get_project_dir(project_dir: None | str | Path) -> Path:
    if project_dir is not None:
        return Path(project_dir)
    else:
        return Path(os.getcwd())


def get_config_path(project_dir) -> Path:
    project_dir = get_project_dir(project_dir)
    return project_dir / "pipeline.toml"


def get_read_group_info_xls(project_dir) -> Path:
    project_dir = get_project_dir(project_dir)
    return project_dir / "reads.xlsx"


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


def get_reads_stats_fastp_parent_dir(project_dir) -> Path:
    return get_reads_stats_parent_dir(project_dir) / "fastp"


def get_reads_stats_fastp_excel_report_path(project_dir) -> Path:
    return get_reads_stats_fastp_parent_dir(project_dir) / "fastp_stats.xlsx"


def get_raw_reads_fastqc_stats_parent_dir(project_dir) -> Path:
    return get_reads_stats_fastqc_parent_dir(project_dir) / "raw"


def get_clean_reads_fastqc_stats_parent_dir(project_dir) -> Path:
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


def _parse_read_file_name(path):
    if not is_read_file(path):
        raise ValueError(f"The given path does not correspond to a read file: {path}")

    file_base_name = path.name.removesuffix(FASTQ_EXT)
    items = file_base_name.split(".")

    if items[-1].startswith("p") and items[-1][1:].isdigit():
        pair_number = int(items[-1][1:])
        base_name = ".".join(items[:-1]) + ".p"
    else:
        pair_number = None
        base_name = ".".join(items)
    return {"base_name": base_name, "pair_number": pair_number}


def get_paired_and_unpaired_read_files_in_dir(dir_path: Path):
    read_files = get_read_files_in_dir(dir_path=dir_path)

    paired_reads = defaultdict(list)
    for read_file in read_files:
        res = _parse_read_file_name(read_file)
        paired_reads[res["base_name"]].append((read_file, res["pair_number"]))

    for pair in paired_reads.values():
        if len(pair) == 1:
            path, pair_number = pair[0]
            if pair_number is None:
                yield (path,)
            else:
                raise ValueError(
                    f"read file seems paired but has no corresponding pair: {path}"
                )
        elif len(pair) == 2:
            pair = sorted(pair, key=lambda path_pairnum: path_pairnum[1])
            yield (pair[0][0], pair[1][0])
        else:
            files = ",".join([str(path_pairnum[0]) for path_pairnum in pair])
            raise ValueError(f"More than two paired read files: {files}")


def get_crams_dir(project_dir) -> Path:
    project_dir = get_project_dir(project_dir=project_dir)
    path = project_dir / "crams"
    return path


def get_crams_stats_excel_report_path(project_dir) -> Path:
    return get_crams_dir(project_dir) / "cram_stats.xlsx"


def get_snv_dir(project_dir) -> Path:
    project_dir = get_project_dir(project_dir=project_dir)
    path = project_dir / "snv_calling"
    return path


def get_gatk_db_dir(project_dir) -> Path:
    return get_snv_dir(project_dir) / "gatk_db"


def get_vcfs_per_sample_dir(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    path = snv_dir / "vcfs_per_sample"
    return path


def get_per_sample_vcfs(project_dir) -> list[Path]:
    vcf_dir = get_vcfs_per_sample_dir(project_dir)
    return [path for path in vcf_dir.iterdir() if str(path).endswith(".vcf.gz")]


def get_joint_vcf(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    return snv_dir / "joint_gatk.vcf.gz"


def get_joint_vcfs_per_segment_dir(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    snv_dir.mkdir(exist_ok=True)
    joint_dir = snv_dir / "joint_vcfs_per_segment"
    joint_dir.mkdir(exist_ok=True)
    return joint_dir


def get_joint_vcfs(project_dir) -> list[Path]:
    vcfs = []
    for path in get_joint_vcfs_per_segment_dir(project_dir).iterdir():
        if str(path).endswith(".vcf.gz"):
            vcfs.append(path)
    return vcfs


def get_joint_gatk_segments_bed(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    snv_dir.mkdir(exist_ok=True)
    return snv_dir / "segments_for_gatk_joint_var_calling.bed"


def get_gatk_intervals_bed(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    snv_dir.mkdir(exist_ok=True)
    return snv_dir / "intervals_for_gatk_db.bed"


def get_joint_var_calling_intervals_bed(project_dir) -> Path:
    snv_dir = get_snv_dir(project_dir)
    snv_dir.mkdir(exist_ok=True)
    return snv_dir / "intervals_for_var_calling.bed"


def get_gatk_interval_db_dir(project_dir, chrom, start, end):
    base_dir = get_gatk_db_dir(project_dir)
    return base_dir / f"{chrom}:{start}-{end}"


def get_gatk_interval_db_dirs(project_dir):
    genome_segment_pattern = re.compile(
        r"""
    ^\.?                # optional leading dot
    (?P<chrom>[^:]+)    # chromosome (everything until :)
    :                   # separator
    (?P<start>\d+)      # start position (digits)
    -                   # separator
    (?P<end>\d+)        # end position (digits)
    \.?$                # optional trailing dot
    """,
        re.VERBOSE,
    )

    base_dir = get_gatk_db_dir(project_dir)
    db_dirs = []
    for dir_ in base_dir.iterdir():
        if not dir_.is_dir():
            continue
        match = genome_segment_pattern.match(dir_.name)
        if not match:
            continue
        match = match.groupdict()
        db_dirs.append(
            {
                "path": dir_,
                "chrom": match["chrom"],
                "start": int(match["start"]),
                "end": int(match["end"]),
            }
        )

    return db_dirs


def get_crams_stats_dir(project_dir) -> Path:
    return get_crams_dir(project_dir) / "stats"


def get_tmp_dir(project_dir) -> Path:
    return get_project_dir(project_dir) / "tmp"


def remove_file(path, not_exist_ok=False):
    if not not_exist_ok and not path.exists():
        raise ValueError(f"Path can't be removed because it doesn't exist: {path}")
    if path.exists():
        os.remove(path)


def is_gzip_file(path: Path) -> bool:
    with open(path, "rb") as f:
        magic = f.read(2)
    return magic == b"\x1f\x8b"


def get_cache_dir(project_dir) -> Path:
    project_dir = get_project_dir(project_dir)
    cache_dir = project_dir / "cache"
    cache_dir.mkdir(exist_ok=True)
    return cache_dir


def get_gvcf_ranges_working_dir(project_dir) -> Path:
    dir_ = get_cache_dir(project_dir) / "gvcf_ranges"
    dir_.mkdir(exist_ok=True)
    return dir_
