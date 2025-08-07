import sys

from split_gvcf.split_genome_avoiding_vars import create_var_ranges
from split_gvcf.split_genome import split_in_empty_loci

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.paths import (
    get_vcfs_per_sample_dir,
    get_gvcf_ranges_working_dir,
    get_joint_gatk_segments_bed,
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
    gatk_vcf_sample_joining_region_size = config["snv_calling"][
        "gatk_vcf_sample_joining_region_size"
    ]
    genome_fai_path = get_genome_fai_path(config)
    chrom_sizes = get_chrom_lengths_from_fai(genome_fai_path)

    vcfs = [
        path
        for path in get_vcfs_per_sample_dir(project_dir).iterdir()
        if str(path).endswith(".vcf.gz")
    ]
    ranges = create_var_ranges(
        vcfs,
        get_gvcf_ranges_working_dir(project_dir),
        config["split_gvcf_vars"]["n_processes_gvcf_parsing"],
    )

    segments = split_in_empty_loci(
        ranges=ranges,
        chrom_sizes=chrom_sizes,
        desired_region_size=gatk_vcf_sample_joining_region_size,
    )

    segments_bed_path = get_joint_gatk_segments_bed(project_dir)
    with segments_bed_path.open("wt") as fhand:
        for chrom, start, end in segments:
            fhand.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()
