"""The SNP calling itself.

pop_var_caller calls the cohort from the psp files, the stored evidence of
every sample, and writes one multi sample VCF with the SNPs, the indels and the
repeat tracts.

The numbers the calling scores with are the ones fitted by the parameter
estimation step.  Calling with the defaults compiled into pop_var_caller is a
different claim, it assumes no base quality calibration, no contamination and
no inbreeding, so it has to be asked for.
"""

from pathlib import Path
import logging
import shutil
import tempfile

from .paths import (
    get_project_dir,
    get_log_path,
    get_psps_dir,
    get_snp_calling_dir,
    get_parameters_path,
    get_vcf_path,
    POP_VAR_CALLER_BIN,
    TABIX_BIN,
)
from .run_cmd import run_cmd, setup_logging
from .psp import get_psp_paths, get_default_repeat_catalog_path

logger = logging.getLogger(__name__)

BGZF_VCF_SUFFIXES = (".vcf.gz", ".vcf.bgz")
# pop_var_caller writes the parameters it called with beside the vcf, the name
# of the vcf without its format suffixes plus this one
PARAMETERS_BESIDE_VCF_SUFFIX = ".parameters.toml"


def is_bgzf_vcf(vcf_path: Path) -> bool:
    return str(vcf_path).endswith(BGZF_VCF_SUFFIXES)


def get_vcf_index_path(vcf_path: Path) -> Path:
    return Path(str(vcf_path) + ".tbi")


def _index_vcf(vcf_path: Path, project_dir: Path, verbose: bool):
    cmd = [TABIX_BIN, "-p", "vcf", str(vcf_path)]
    run_cmd(cmd, project_dir=project_dir, verbose=verbose)

    index_path = get_vcf_index_path(vcf_path)
    if not index_path.exists():
        msg = f"tabix ran, but it did not create the vcf index: {index_path}"
        logging.error(msg)
        raise RuntimeError(msg)
    return index_path


def call_snps(
    project_dir: Path,
    genome_fasta: Path,
    catalog_path: None | Path = None,
    parameters_path: None | Path = None,
    use_default_parameters: bool = False,
    ploidy: int = 2,
    num_threads: int = 0,
    paralog_fdr: None | float = None,
    paralog_filter_tag: bool = False,
    max_cohort_locus_span: int = 0,
    max_candidate_alleles: int = 0,
    cohort_locus_builder_regions_len: int = 0,
    min_copies: str = "",
    min_period: int = 0,
    max_period: int = 0,
    max_str_len: int = 0,
    min_purity: float = 0.0,
    index_vcf: bool = True,
    re_run: bool = False,
    verbose: bool = True,
):
    project_dir = get_project_dir(project_dir)

    setup_logging(project_dir)

    genome_fasta = Path(genome_fasta)
    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome fasta file not found: {genome_fasta}")
    genome_fasta = genome_fasta.resolve()

    if catalog_path is None:
        catalog_path = get_default_repeat_catalog_path(genome_fasta)
    else:
        catalog_path = Path(catalog_path)
    if not catalog_path.exists():
        raise FileNotFoundError(
            f"The tandem repeat catalog was not found: {catalog_path}\n"
            f"Create it with: {POP_VAR_CALLER_BIN} repeat-catalog --reference {genome_fasta}"
        )
    catalog_path = catalog_path.resolve()

    # Calling with the compiled in defaults assumes no base quality
    # calibration, no contamination and no inbreeding, so it is not something
    # to do just because the parameters have not been fitted yet
    if use_default_parameters:
        parameters_path = None
    else:
        if parameters_path is None:
            parameters_path = get_parameters_path(project_dir)
        else:
            parameters_path = Path(parameters_path)
        if not parameters_path.exists():
            raise FileNotFoundError(
                f"The parameters file was not found: {parameters_path}\n"
                "Fit it with the estimate_parameters script, or ask for the pop_var_caller "
                "defaults, which assume no base quality calibration, no contamination and no inbreeding"
            )
        parameters_path = parameters_path.resolve()

    vcf_path = get_vcf_path(project_dir)
    if vcf_path.exists() and not re_run:
        if verbose:
            print(f"Skipping analysis for existing vcf file: {vcf_path}")
        return {"vcf_path": vcf_path, "should_have_run": False}

    psps_dir = get_psps_dir(project_dir)
    psp_paths = get_psp_paths(project_dir)
    if not psp_paths:
        raise RuntimeError(
            f"There is no psp file to call in {psps_dir}, generate them first"
        )
    if verbose:
        print(f"Num. samples to call: {len(psp_paths)}")

    snp_calling_dir = get_snp_calling_dir(project_dir)
    snp_calling_dir.mkdir(parents=True, exist_ok=True)

    # The cohort is called in a temporary dir inside the snp calling dir and the
    # vcf is moved, once it is complete, to its final path.  Both paths are in
    # the same file system, so the move is a rename, and an interrupted run
    # leaves either the vcf of the previous run or no vcf at all, but never a
    # partial one.  The temporary files that the calling writes beside its
    # output, the vcf that is being written and the spill of the hidden
    # duplication filter, are left in the temporary dir, which is removed.
    with tempfile.TemporaryDirectory(
        prefix="vcf_tmp_dir_", dir=snp_calling_dir
    ) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        tmp_vcf_path = tmp_dir_path / vcf_path.name

        cmd = [POP_VAR_CALLER_BIN, "call-from-psps"]
        cmd.extend(["--reference", str(genome_fasta)])
        cmd.extend(["--catalog", str(catalog_path)])
        cmd.extend(["--psp", str(psps_dir)])
        cmd.extend(["--output", str(tmp_vcf_path)])
        if parameters_path is None:
            cmd.append("--defaults")
        else:
            cmd.extend(["--parameters", str(parameters_path)])
        cmd.extend(["--ploidy", str(ploidy)])
        if paralog_fdr is not None:
            cmd.extend(["--paralog-fdr", str(paralog_fdr)])
        if paralog_filter_tag:
            cmd.append("--paralog-filter-tag")
        if max_cohort_locus_span:
            cmd.extend(["--max-cohort-locus-span", str(max_cohort_locus_span)])
        if max_candidate_alleles:
            cmd.extend(["--max-candidate-alleles", str(max_candidate_alleles)])
        if cohort_locus_builder_regions_len:
            cmd.extend(
                [
                    "--cohort-locus-builder-regions-len",
                    str(cohort_locus_builder_regions_len),
                ]
            )
        # The repeat criteria have to be the ones the psps were walked with, or
        # the cohort is refused, so they are taken from the same config keys
        if min_copies:
            cmd.extend(["--min-copies", str(min_copies)])
        if min_period:
            cmd.extend(["--min-period", str(min_period)])
        if max_period:
            cmd.extend(["--max-period", str(max_period)])
        if max_str_len:
            cmd.extend(["--max-str-len", str(max_str_len)])
        if min_purity:
            cmd.extend(["--min-purity", str(min_purity)])
        # Unlike the other steps, call-from-psps has its own threads option,
        # which builds the rayon global pool, so no environment is required
        # here.  The vcf written does not depend on how many threads are used
        if num_threads:
            cmd.extend(["--threads", str(num_threads)])

        logging.info("Calling the cohort snps")
        run_cmd(cmd, project_dir=project_dir, verbose=verbose, stream_output=True)

        if not tmp_vcf_path.exists():
            msg = f"pop_var_caller call-from-psps ran, but it did not create the expected vcf file: {tmp_vcf_path}"
            logging.error(msg)
            raise RuntimeError(msg)

        if index_vcf and is_bgzf_vcf(tmp_vcf_path):
            _index_vcf(tmp_vcf_path, project_dir=project_dir, verbose=verbose)

        # Every run writes the parameters that it called with beside its vcf,
        # whatever the numbers came from, so that file is moved as well.  The
        # temporary files that the calling may have left behind are not.
        moved_paths = []
        for path in sorted(tmp_dir_path.iterdir()):
            if path.name.endswith(".tmp"):
                continue
            dest_path = snp_calling_dir / path.name
            shutil.move(path, dest_path)
            moved_paths.append(dest_path)

    if verbose:
        print(f"VCF written in: {vcf_path}")
    return {
        "vcf_path": vcf_path,
        "created_paths": moved_paths,
        "should_have_run": True,
    }
