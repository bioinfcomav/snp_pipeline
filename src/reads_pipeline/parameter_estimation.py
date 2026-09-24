"""Estimation of the parameters used as priors by the SNP calling.

pop_var_caller fits the cohort numbers -the per library sequencing error rates
and base quality calibrations, the allele frequency density and the genotype
prior seeded from it, each library contamination and the repeat tract slippage
ladder- from the censuses stored in the psp files, and writes them in a
parameters file that the calling run reads.

The sequencing batches, the sets of experiments that were sequenced together,
are not in the alignments, so they have to be stated.  We assume that the
experiments of one project, the dir in which the reads and the crams are
stored, were sequenced together.
"""

from pathlib import Path
import logging
import shutil
import tempfile
import tomllib

from .paths import (
    get_project_dir,
    get_log_path,
    get_psps_dir,
    get_snp_calling_dir,
    get_parameters_path,
    POP_VAR_CALLER_BIN,
)
from .run_cmd import run_cmd
from .psp import (
    get_cram_paths,
    get_psp_paths,
    get_default_repeat_catalog_path,
    get_rayon_env,
)
from .read_group import get_read_groups_in_cram_header

logger = logging.getLogger(__name__)

SEQUENCING_BATCHES_SECTION_NAME = "sequencing_batches"


def get_batch_of_each_read_group(project_dir) -> dict[str, str]:
    """The sequencing batch of every read group.

    We assume that the experiments of one project, the dir that holds their
    crams, were sequenced together.

    The read groups are taken from the @RG ID tags of the cram headers, which
    is the name that pop_var_caller writes in the parameters file as the
    declared_id of every read group it fits.
    """
    batches = {}
    crams_with_no_read_group = []
    for cram_path in get_cram_paths(project_dir):
        batch = cram_path.parent.name
        read_group_ids = [
            read_group["ID"]
            for read_group in get_read_groups_in_cram_header(cram_path, project_dir)
            if read_group.get("ID", "")
        ]
        if not read_group_ids:
            crams_with_no_read_group.append(cram_path)
            continue
        for read_group_id in read_group_ids:
            if batches.setdefault(read_group_id, batch) != batch:
                msg = (
                    f"The read group {read_group_id} is in the crams of two projects: "
                    f"{batches[read_group_id]} and {batch}, so it is not known which one it was sequenced in"
                )
                logging.error(msg)
                raise RuntimeError(msg)

    if crams_with_no_read_group:
        msg = (
            "The sequencing batch of a read group is the project dir that holds "
            f"its cram, and {len(crams_with_no_read_group)} cram(s) have no @RG header line:\n"
        )
        msg += "".join(f"  {path}\n" for path in crams_with_no_read_group)
        logging.error(msg)
        raise RuntimeError(msg)

    return batches


def _toml_str(string: str) -> str:
    escaped = string.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _get_read_groups_in_parameters_file(parameters: dict, parameters_path: Path):
    try:
        read_groups = parameters["fitted_from"]["read_groups"]
    except KeyError:
        raise RuntimeError(
            f"The parameters file has no fitted_from.read_groups table: {parameters_path}"
        )
    for read_group in read_groups:
        if "read_group" not in read_group or "declared_id" not in read_group:
            raise RuntimeError(
                f"Malformed fitted_from.read_groups row in the parameters file: {read_group}"
            )
    return read_groups


def _create_sequencing_batches_section(
    parameters: dict, batch_of_read_group: dict[str, str], parameters_path: Path
) -> str:
    read_groups = _get_read_groups_in_parameters_file(parameters, parameters_path)

    batch_names = []
    for read_group in read_groups:
        declared_id = read_group["declared_id"]
        if declared_id not in batch_of_read_group:
            raise RuntimeError(
                f"The read group {declared_id}, fitted in the parameters file, has no cram in this project, so its sequencing batch is unknown"
            )
        batch_names.append(batch_of_read_group[declared_id])
    # Only the batches of the read groups that were fitted are numbered, so the
    # batch ids are the dense 0..n index that pop_var_caller sizes its tables for
    batch_ids = {name: idx for idx, name in enumerate(sorted(set(batch_names)))}

    by_read_group_rows = []
    batch_name_of_sample = {}
    for read_group, batch_name in zip(read_groups, batch_names):
        by_read_group_rows.append(
            f"    {{ read_group = {read_group['read_group']}, batch = {batch_ids[batch_name]} }},"
        )

        # A sample sequenced in several batches would be refused by
        # pop_var_caller, because a contaminating read is drawn against one set
        # of neighbours, and a sample has one genotype to draw
        sample = read_group["sample"]
        if batch_name_of_sample.setdefault(sample, batch_name) != batch_name:
            raise RuntimeError(
                f"The sample {sample} has read groups sequenced in different projects: "
                f"{batch_name_of_sample[sample]} and {batch_name}. "
                "A sample can only be in one sequencing batch"
            )

    by_sample_rows = [
        f"    {{ sample = {_toml_str(sample)}, batch = {batch_ids[batch_name]} }},"
        for sample, batch_name in batch_name_of_sample.items()
    ]

    batch_comments = "\n".join(
        f"#   batch {batch_id}: {name}" for name, batch_id in batch_ids.items()
    )
    by_read_group_rows = "\n".join(by_read_group_rows)
    by_sample_rows = "\n".join(by_sample_rows)

    section = f"""[{SEQUENCING_BATCHES_SECTION_NAME}]
# Declared by the reads pipeline: we assume that the experiments of one project
# were sequenced together, so every project is one sequencing batch.
{batch_comments}
batching_was_declared = true
by_read_group = [
{by_read_group_rows}
]
by_sample = [
{by_sample_rows}
]
"""
    return section


def _remove_sequencing_batches_section(text: str) -> str:
    """Removes the sequencing batches section, however it is written.

    The section can be written as one table with inline rows or as several
    [[sequencing_batches.by_read_group]] tables, so every block whose header
    talks about the sequencing batches is removed.
    """
    kept_lines = []
    in_section = False
    for line in text.splitlines(keepends=True):
        if line.startswith("["):
            header = line.strip().strip("[]")
            in_section = header == SEQUENCING_BATCHES_SECTION_NAME or header.startswith(
                SEQUENCING_BATCHES_SECTION_NAME + "."
            )
        if not in_section:
            kept_lines.append(line)
    return "".join(kept_lines)


def declare_sequencing_batches(parameters_path: Path, batch_of_read_group: dict):
    """Writes the sequencing batches in a parameters file.

    The batching is not in the alignments, so pop_var_caller can not infer it,
    it takes it from the parameters file that the calling run reads.
    """
    text = parameters_path.read_text()
    parameters = tomllib.loads(text)

    section = _create_sequencing_batches_section(
        parameters, batch_of_read_group, parameters_path
    )

    text = _remove_sequencing_batches_section(text)
    if not text.endswith("\n"):
        text += "\n"
    text += "\n" + section

    # The written file has to be a parameters file that can still be read
    new_parameters = tomllib.loads(text)
    batches = new_parameters[SEQUENCING_BATCHES_SECTION_NAME]
    if not batches["batching_was_declared"]:
        raise RuntimeError("The sequencing batches were not declared")
    if len(batches["by_read_group"]) != len(parameters["fitted_from"]["read_groups"]):
        raise RuntimeError(
            "A run that names any sequencing batch must name them all, but some read group was left out"
        )

    parameters_path.write_text(text)
    return {
        "num_batches": len({row["batch"] for row in batches["by_read_group"]}),
        "num_read_groups": len(batches["by_read_group"]),
        "num_samples": len(batches["by_sample"]),
    }


def estimate_parameters(
    project_dir: Path,
    genome_fasta: Path,
    catalog_path: None | Path = None,
    ploidy: int = 2,
    inbreeding: None | float = None,
    declare_batches: bool = True,
    num_threads: int = 0,
    re_run: bool = False,
    verbose: bool = True,
):
    project_dir = get_project_dir(project_dir)

    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        force=True,
    )

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

    parameters_path = get_parameters_path(project_dir)
    if parameters_path.exists() and not re_run:
        if verbose:
            print(f"Skipping analysis for existing parameters file: {parameters_path}")
        return {"parameters_path": parameters_path, "should_have_run": False}

    psps_dir = get_psps_dir(project_dir)
    psp_paths = get_psp_paths(project_dir)
    if not psp_paths:
        raise RuntimeError(
            f"There is no psp file to fit the parameters from in {psps_dir}, generate them first"
        )
    if verbose:
        print(f"Num. psps to fit the parameters from: {len(psp_paths)}")

    batch_of_read_group = (
        get_batch_of_each_read_group(project_dir) if declare_batches else None
    )

    snp_calling_dir = get_snp_calling_dir(project_dir)
    snp_calling_dir.mkdir(parents=True, exist_ok=True)

    # The parameters are fitted in a temporary dir inside the snp calling dir
    # and the file is moved, once it is complete, to its final path.  Both paths
    # are in the same file system, so the move is a rename, and an interrupted
    # run leaves either the parameters of the previous run or no parameters at
    # all, but never a partial file.
    with tempfile.TemporaryDirectory(
        prefix="parameters_tmp_dir_", dir=snp_calling_dir
    ) as tmp_dir:
        tmp_parameters_path = Path(tmp_dir) / parameters_path.name

        cmd = [POP_VAR_CALLER_BIN, "estimate-parameters"]
        cmd.extend(["--reference", str(genome_fasta)])
        cmd.extend(["--catalog", str(catalog_path)])
        cmd.extend(["--psp", str(psps_dir)])
        cmd.extend(["--output", str(tmp_parameters_path)])
        cmd.extend(["--ploidy", str(ploidy)])
        if inbreeding is not None:
            cmd.extend(["--inbreeding", str(inbreeding)])

        logging.info("Estimating the SNP calling parameters")
        run_cmd(
            cmd,
            project_dir=project_dir,
            verbose=verbose,
            env=get_rayon_env(num_threads),
        )

        if not tmp_parameters_path.exists():
            msg = f"pop_var_caller estimate-parameters ran, but it did not create the expected parameters file: {tmp_parameters_path}"
            logging.error(msg)
            raise RuntimeError(msg)

        if declare_batches:
            res = declare_sequencing_batches(tmp_parameters_path, batch_of_read_group)
            if verbose:
                print(
                    f"Sequencing batches declared: {res['num_batches']} batches for {res['num_read_groups']} read groups and {res['num_samples']} samples"
                )

        shutil.move(tmp_parameters_path, parameters_path)

    if verbose:
        print(f"Parameters written in: {parameters_path}")
    return {"parameters_path": parameters_path, "should_have_run": True}
