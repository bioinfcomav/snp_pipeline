"""Generation of the per sample pileup (psp) files used by pop_var_caller.

The first step of the SNP calling is to walk the alignments of every sample and
to store what its reads showed at each position in a psp file.  pop_var_caller
walks the samples given to one invocation one after the other, so the
parallelization is done here, by running one invocation per sample.
"""

from pathlib import Path
import logging
import shutil
import tempfile
from functools import partial
from multiprocessing import Pool

from .paths import (
    get_project_dir,
    get_log_path,
    get_crams_dir,
    get_psps_dir,
    get_psp_path,
    CRAM_EXT,
    PSP_EXT,
    POP_VAR_CALLER_BIN,
)
from .run_cmd import run_cmd, setup_logging
from .read_group import get_samples_in_cram

logger = logging.getLogger(__name__)

REPEAT_CATALOG_SUFFIX = ".repeats.parquet"
# pop_var_caller parallelizes with rayon.  generate-psps and estimate-parameters
# take no --threads option, so the only way of narrowing them is the environment
# variable that sets the size of the rayon global thread pool
RAYON_NUM_THREADS_ENV_VAR = "RAYON_NUM_THREADS"


def get_rayon_env(num_threads: int) -> None | dict:
    "The environment that bounds the pop_var_caller threads, zero means every core"
    if not num_threads:
        return None
    return {RAYON_NUM_THREADS_ENV_VAR: str(num_threads)}


def get_default_repeat_catalog_path(genome_fasta: Path) -> Path:
    "The path in which pop_var_caller repeat-catalog writes the catalog"
    return Path(str(genome_fasta) + REPEAT_CATALOG_SUFFIX)


def get_cram_paths(project_dir) -> list[Path]:
    crams_parent_dir = get_crams_dir(project_dir)
    if not crams_parent_dir.exists():
        raise RuntimeError(
            f"The crams dir does not exist, map the reads first: {crams_parent_dir}"
        )

    cram_paths = []
    for dir_ in sorted(crams_parent_dir.iterdir()):
        if not dir_.is_dir():
            continue
        cram_paths.extend(sorted(dir_.glob("*" + CRAM_EXT)))
    return cram_paths


def _check_sample_name_can_be_a_file_name(sample: str):
    # pop_var_caller names every psp after the @RG SM tag, so a sample named
    # like a path would write outside the psps dir
    if not sample or sample in (".", "..") or sample != Path(sample).name:
        raise ValueError(
            f"The sample name can not be used as a file name: {sample}, "
            "the psp files are named <sample>.psp"
        )


def get_samples_to_process(project_dir) -> list[dict]:
    """The crams to walk for every sample.

    All the crams that share a sample are walked together into one psp,
    because they are the alignments of that one sample.

    The sample of every cram is read from the SM tag of its own @RG header
    lines, which is what the mapping step wrote there, so the read group excel
    file is not required for this step.
    """
    cram_paths = get_cram_paths(project_dir)

    crams_per_sample = {}
    crams_with_no_sample = []
    crams_with_several_samples = []
    for cram_path in cram_paths:
        samples = get_samples_in_cram(cram_path, project_dir)
        if not samples:
            crams_with_no_sample.append(cram_path)
            continue
        if len(samples) > 1:
            crams_with_several_samples.append((cram_path, sorted(samples)))
            continue
        sample = samples.pop()
        _check_sample_name_can_be_a_file_name(sample)
        crams_per_sample.setdefault(sample, []).append(cram_path)

    if crams_with_no_sample or crams_with_several_samples:
        msg = ""
        if crams_with_no_sample:
            msg += (
                "The sample of a cram is read from the SM tag of its @RG header "
                f"lines, and {len(crams_with_no_sample)} cram(s) have no SM tag:\n"
            )
            msg += "".join(f"  {path}\n" for path in crams_with_no_sample)
            msg += (
                "Map those reads again, or add the sample to the cram header with: "
                "samtools addreplacerg\n"
            )
        if crams_with_several_samples:
            msg += (
                f"{len(crams_with_several_samples)} cram(s) hold several samples, "
                "so it is not known which psp they belong to:\n"
            )
            msg += "".join(
                f"  {path}: {', '.join(samples)}\n"
                for path, samples in crams_with_several_samples
            )
        logging.error(msg)
        raise RuntimeError(msg)

    samples_to_process = []
    for idx, sample in enumerate(sorted(crams_per_sample.keys()), start=1):
        samples_to_process.append(
            {
                "idx": idx,
                "sample": sample,
                "cram_paths": crams_per_sample[sample],
            }
        )
    return samples_to_process


def _generate_psp_for_sample(
    sample_info: dict,
    project_dir: Path,
    psps_dir: Path,
    genome_fasta: Path,
    catalog_path: Path,
    regions_path: None | Path,
    min_copies: str,
    min_period: int,
    max_period: int,
    max_str_len: int,
    min_purity: float,
    build_index_if_missing: bool,
    num_threads: int,
    re_run: bool,
    verbose: bool,
    num_analyses_to_do: int,
):
    sample = sample_info["sample"]
    cram_paths = sample_info["cram_paths"]
    psp_path = get_psp_path(project_dir, sample)

    setup_logging(project_dir)

    if psp_path.exists() and not re_run:
        if verbose:
            print(f"Skipping analysis for existing psp file: {psp_path}")
        return {"psp_path": psp_path, "should_have_run": False}

    if verbose:
        crams_str = ", ".join(map(str, cram_paths))
        print(
            f"Generating psp for sample with idx: {sample_info['idx']}, total to process: {num_analyses_to_do} : {sample} ({crams_str})"
        )
    logging.info(
        f"Generating psp for sample {sample} from crams: "
        + " ".join(map(str, cram_paths))
    )

    # The psp is written in a temporary dir inside the psps dir and it is moved,
    # once it is complete, to its final path.  Both paths are in the same file
    # system, so the move is a rename, and an interrupted run leaves either the
    # psp of the previous run or no psp at all, but never a partial one.
    with tempfile.TemporaryDirectory(prefix="psp_tmp_dir_", dir=psps_dir) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)

        cmd = [POP_VAR_CALLER_BIN, "generate-psps"]
        cmd.extend(["--reference", str(genome_fasta)])
        cmd.extend(["--catalog", str(catalog_path)])
        for cram_path in cram_paths:
            cmd.extend(["--alignment", str(cram_path)])
        cmd.extend(["--output-dir", str(tmp_dir_path)])
        if regions_path is not None:
            cmd.extend(["--regions", str(regions_path)])
        if build_index_if_missing:
            cmd.append("--build-index-if-missing")
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

        run_cmd(
            cmd,
            project_dir=project_dir,
            verbose=verbose,
            env=get_rayon_env(num_threads),
            stream_output=True,
        )

        tmp_psp_path = tmp_dir_path / f"{sample}{PSP_EXT}"
        if not tmp_psp_path.exists():
            msg = f"pop_var_caller generate-psps ran, but it did not create the expected psp file: {tmp_psp_path}"
            logging.error(msg)
            raise RuntimeError(msg)

        shutil.move(tmp_psp_path, psp_path)

    return {"psp_path": psp_path, "should_have_run": True}


def generate_psps_for_samples(
    project_dir: Path,
    genome_fasta: Path,
    catalog_path: None | Path = None,
    regions_path: None | Path = None,
    min_copies: str = "",
    min_period: int = 0,
    max_period: int = 0,
    max_str_len: int = 0,
    min_purity: float = 0.0,
    build_index_if_missing: bool = False,
    num_threads: int = 0,
    re_run: bool = False,
    verbose: bool = True,
    num_psps_in_parallel: int = 1,
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

    if regions_path is not None:
        regions_path = Path(regions_path)
        if not regions_path.exists():
            raise FileNotFoundError(
                f"The regions bed file was not found: {regions_path}"
            )
        regions_path = regions_path.resolve()

    psps_dir = get_psps_dir(project_dir)
    psps_dir.mkdir(parents=True, exist_ok=True)

    samples_to_process = get_samples_to_process(project_dir)
    if verbose:
        print(f"Num. samples: {len(samples_to_process)}")

    if re_run:
        samples_todo = samples_to_process
    else:
        samples_todo = [
            sample_info
            for sample_info in samples_to_process
            if not get_psp_path(project_dir, sample_info["sample"]).exists()
        ]
    num_analyses_to_do = len(samples_todo)
    if verbose:
        print(f"Num. psps to generate: {num_analyses_to_do}")

    generate_psp_for_sample = partial(
        _generate_psp_for_sample,
        project_dir=project_dir,
        psps_dir=psps_dir,
        genome_fasta=genome_fasta,
        catalog_path=catalog_path,
        regions_path=regions_path,
        min_copies=min_copies,
        min_period=min_period,
        max_period=max_period,
        max_str_len=max_str_len,
        min_purity=min_purity,
        build_index_if_missing=build_index_if_missing,
        num_threads=num_threads,
        re_run=re_run,
        verbose=verbose,
        num_analyses_to_do=num_analyses_to_do,
    )

    if num_psps_in_parallel > 1:
        with Pool(num_psps_in_parallel) as pool:
            results = pool.map(generate_psp_for_sample, samples_todo)
    else:
        results = map(generate_psp_for_sample, samples_todo)

    psp_paths = []
    num_analyses_done = 0
    for res in results:
        psp_paths.append(res["psp_path"])
        num_analyses_done += int(res["should_have_run"])

    return {"psp_paths": psp_paths, "num_analyses_done": num_analyses_done}


def get_psp_paths(project_dir) -> list[Path]:
    psps_dir = get_psps_dir(project_dir)
    if not psps_dir.exists():
        return []
    return sorted(path for path in psps_dir.iterdir() if path.suffix == PSP_EXT)
