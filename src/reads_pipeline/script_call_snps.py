import os
import argparse
import sys

from reads_pipeline.paths import get_config_path, get_project_dir
from reads_pipeline.snp_calling import call_snps
from reads_pipeline.pipeline_config import PipelineConfig


def get_args():
    parser = argparse.ArgumentParser(
        description="Calls the SNPs, indels and repeat tracts of the cohort from the psp files."
    )

    # project_dir: Defaults to the current working directory
    parser.add_argument(
        "project_dir",
        nargs="?",  # Makes it optional
        default=os.getcwd(),
        help="Project directory (default: current working directory).",
    )

    return parser.parse_args()


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
        sys.exit(2)

    config = PipelineConfig(project_dir=project_dir)

    catalog_path = config["pop_var_caller"]["catalog_path"]
    paralog_fdr = config["pop_var_caller"]["paralog_fdr"]

    call_snps(
        project_dir=project_dir,
        genome_fasta=config["general"]["genome_path"],
        catalog_path=catalog_path if catalog_path else None,
        use_default_parameters=config["pop_var_caller"]["use_default_parameters"],
        ploidy=config["pop_var_caller"]["ploidy"],
        num_threads=config["pop_var_caller"]["calling_num_threads"],
        # An empty paralog_fdr means that the pop_var_caller default is used, a
        # zero turns the hidden duplication filter off
        paralog_fdr=float(paralog_fdr) if paralog_fdr != "" else None,
        paralog_filter_tag=config["pop_var_caller"]["paralog_filter_tag"],
        max_cohort_locus_span=config["pop_var_caller"]["max_cohort_locus_span"],
        max_candidate_alleles=config["pop_var_caller"]["max_candidate_alleles"],
        cohort_locus_builder_regions_len=config["pop_var_caller"][
            "cohort_locus_builder_regions_len"
        ],
        # The repeat criteria have to be the ones the psps were walked with
        min_copies=config["pop_var_caller"]["min_copies"],
        min_period=config["pop_var_caller"]["min_period"],
        max_period=config["pop_var_caller"]["max_period"],
        max_str_len=config["pop_var_caller"]["max_str_len"],
        min_purity=config["pop_var_caller"]["min_purity"],
        index_vcf=config["pop_var_caller"]["index_vcf"],
        re_run=config["general"]["re_run"],
        verbose=config["general"]["verbose"],
    )


if __name__ == "__main__":
    main()
