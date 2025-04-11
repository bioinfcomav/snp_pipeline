import sys

from reads_pipeline.script_check_mapping import get_args
from reads_pipeline.paths import (
    get_project_dir,
    get_reads_stats_fastp_excel_report_path,
    get_crams_stats_dir,
    get_crams_stats_excel_report_path,
)
from reads_pipeline.fastp_minimap import get_fastq_pairs_to_process, collect_cram_stats
from reads_pipeline.fastp import collect_fastp_stats


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    pairs_to_process = get_fastq_pairs_to_process(project_dir)
    print(f"Num. pairs to process: {len(pairs_to_process)}")

    fastp_stats = collect_fastp_stats(project_dir=project_dir)
    fastp_stats.to_excel(
        get_reads_stats_fastp_excel_report_path(project_dir), index=False
    )
    print(f"Num. fastp stats: {fastp_stats.shape[0]}")

    cram_stats = collect_cram_stats(project_dir=project_dir)
    stats_dir = get_crams_stats_dir(project_dir)
    stats_dir.mkdir(exist_ok=True)
    stats_path = get_crams_stats_excel_report_path(project_dir)
    cram_stats.to_excel(stats_path, index=False)
    print(f"Num. cram stats: {cram_stats.shape[0]}")
