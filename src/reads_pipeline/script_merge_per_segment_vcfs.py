import sys
from pathlib import Path

from reads_pipeline.script_run_mapping import get_args, get_project_dir, get_config_path
from reads_pipeline.pipeline_config import PipelineConfig
from reads_pipeline.paths import get_joint_vcfs, BCFTOOLS_BIN, get_joint_vcf
from reads_pipeline.run_cmd import run_bash_script

# bcftools concat -o - SL4.0ch12\:28100001-28200000.joint.vcf.gz SL4.0ch12\:28200001-28300000.joint.vcf.gz
# SL4.0ch12\:40000001-40100000.joint.vcf.gz |
# bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<10'  |
# bcftools filter --include "F_MISSING<0.8 & QUAL>20 & INFO/MAF>=0.01"  | wc -l
PIPE_TEMPLATE = """
{bcf_bin} concat -o - {vcfs} | \\
{gt_filter_str}
{snv_filter_str}
gzip > {out_vcf}
"""
JOINT_VCF_EXT = ".joint.vcf.gz"


def _get_chrom_start_end_from_joint_vcf(vcf: Path):
    if vcf.name.endswith(JOINT_VCF_EXT):
        chrom_poss = vcf.name[: -len(JOINT_VCF_EXT)]
    else:
        raise RuntimeError(
            f"Unexpected name for joint VCF, extension is not {JOINT_VCF_EXT}"
        )
    chrom, poss = chrom_poss.split(":")
    start, end = poss.split("-")
    start = int(start)
    end = int(end)
    return chrom, start, end


def merge_per_segment_vcfs(
    project_dir,
    gt_min_depth=None,
    gt_min_qual=None,
    snv_max_missing_rate=None,
    snv_min_qual=None,
    snv_min_maf=None,
):
    vcfs = get_joint_vcfs(project_dir)
    vcfs = sorted(vcfs, key=lambda vcf: _get_chrom_start_end_from_joint_vcf(vcf)[:2])

    vcfs = " ".join(map(str, vcfs))

    if gt_min_depth and gt_min_qual:
        gt_filter_str = f'{BCFTOOLS_BIN} filter -S . -e "FMT/DP<{gt_min_depth} | FMT/GQ<{gt_min_qual}" | \\'
    elif gt_min_depth and not gt_min_qual:
        gt_filter_str = f'{BCFTOOLS_BIN} filter -S . -e "FMT/DP<{gt_min_depth}" | \\'
    elif not gt_min_depth and gt_min_qual:
        gt_filter_str = f'{BCFTOOLS_BIN} filter -S . -e "FMT/GQ<{gt_min_qual}" | \\'
    else:
        gt_filter_str = ""

    snv_filters = []
    if snv_max_missing_rate:
        snv_filters.append(f"F_MISSING<{snv_max_missing_rate}")
    if snv_min_qual:
        snv_filters.append(f"QUAL>{snv_min_qual}")
    if snv_min_maf:
        snv_filters.append(f"INFO/MAF>{snv_min_maf}")
    if snv_filters:
        snv_filter_str = (
            f'{BCFTOOLS_BIN} filter --include "' + "&".join(snv_filters) + '" | \\'
        )
    else:
        snv_filter_str = ""

    pipe = PIPE_TEMPLATE.format(
        bcf_bin=BCFTOOLS_BIN,
        vcfs=vcfs,
        gt_filter_str=gt_filter_str,
        snv_filter_str=snv_filter_str,
        out_vcf=str(get_joint_vcf(project_dir)),
    )

    run_bash_script(pipe, project_dir=project_dir)


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
    merge_per_segment_vcfs(
        project_dir,
        gt_min_depth=config["merge_vcf_segments"]["gt_min_depth"],
        gt_min_qual=config["merge_vcf_segments"]["gt_min_qual"],
        snv_max_missing_rate=config["merge_vcf_segments"]["snv_max_missing_rate"],
        snv_min_qual=config["merge_vcf_segments"]["snv_min_qual"],
        snv_min_maf=config["merge_vcf_segments"]["snv_min_maf"],
    )


if __name__ == "__main__":
    main()
