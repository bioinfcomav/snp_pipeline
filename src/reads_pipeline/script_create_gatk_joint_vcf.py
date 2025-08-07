import sys
import tempfile
from pathlib import Path

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.gatk import (
    do_svn_joint_genotyping_for_all_samples_together,
    filter_vcf_with_gatk,
)
from reads_pipeline.paths import (
    get_joint_vcf,
    get_snv_dir,
    get_joint_vcfs_per_segment_dir,
    get_joint_gatk_segments_bed,
)


def _read_segments_from_bed(project_dir):
    segment_bed = get_joint_gatk_segments_bed(project_dir)
    with segment_bed.open("rt") as fhand:
        for line in fhand:
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            yield chrom, start, end


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)
    config_path = get_config_path(project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    if not config_path.exists():
        msg = f"The config file {config_path} does not exist"
        print(msg)
        sys.exit(3)

    config = PipelineConfig(project_dir=project_dir)

    genome_path = config["general"]["genome_path"]

    out_vcf_per_segment_dir = get_joint_vcfs_per_segment_dir(project_dir)
    segments = _read_segments_from_bed(project_dir)

    segments = list(segments)
    bigger_start = max([segment[1] for segment in segments])
    num_digits = len(str(bigger_start))
    vcf_name_formatting = (
        f"{{chrom}}-{{start:0{num_digits}}}-{{end:0{num_digits}}}.vcf.gz"
    )

    for chrom, start, end in segments:
        out_vcf = out_vcf_per_segment_dir / vcf_name_formatting.format(
            chrom=chrom, start=start, end=end
        )

        if out_vcf.exists():
            print(f"VCF already created for segment: {out_vcf}")
            continue

        join_vcf_tmp = tempfile.NamedTemporaryFile(
            prefix="gatk.join.", suffix=".vcf.gz", dir=get_snv_dir(project_dir)
        )

        do_svn_joint_genotyping_for_all_samples_together(
            project_dir=project_dir,
            genome_fasta=genome_path,
            out_vcf=Path(join_vcf_tmp.name),
        )

        filter_vcf_with_gatk(
            in_vcf=Path(join_vcf_tmp.name),
            out_vcf=out_vcf,
            genome_fasta=genome_path,
            filters=config["gatk_filters"],
            project_dir=project_dir,
        )
        join_vcf_tmp.close()


if __name__ == "__main__":
    main()
