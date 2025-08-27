import sys

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.paths import (
    get_joint_var_calling_intervals_bed,
)
from reads_pipeline.gatk import get_chrom_lengths_from_fai, get_genome_fai_path


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)
    config_path = get_config_path(project_dir)

    if not project_dir.exists():
        msg = f"The project directory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    if not config_path.exists():
        msg = f"The config file {config_path} does not exist"
        print(msg)
        sys.exit(3)

    config = PipelineConfig(project_dir=project_dir)
    genome_fai_path = get_genome_fai_path(config)
    chrom_sizes = get_chrom_lengths_from_fai(genome_fai_path)
    join_segments_bed = get_joint_var_calling_intervals_bed(project_dir)

    snv_join_calling_segment_sizes = config["snv_calling"][
        "snv_join_calling_segment_sizes"
    ]

    with join_segments_bed.open("wt") as fhand:
        for chrom, chrom_length in chrom_sizes.items():
            for start in range(1, chrom_length + 1, snv_join_calling_segment_sizes):
                end = min(start + snv_join_calling_segment_sizes - 1, chrom_length)
                fhand.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()
