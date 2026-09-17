import logging
from pathlib import Path

import pandas

from reads_pipeline.paths import get_read_group_info_xls, get_log_path, SAMTOOLS_BIN
from reads_pipeline.run_cmd import run_cmd


def get_read_group_info(project_dir) -> dict:
    xls_path = get_read_group_info_xls(project_dir)
    if not xls_path.exists():
        raise RuntimeError(f"Missing excel file with read group info: {xls_path}")

    read_groups_df = pandas.read_excel(xls_path, index_col="id", dtype=str).fillna("")
    if "sample" not in read_groups_df.columns:
        raise RuntimeError("Malformed read group excel file, missing column: id")

    if not read_groups_df.index.is_unique:
        duplicates = read_groups_df.index[read_groups_df.index.duplicated()]
        raise RuntimeError(
            f"Malformed read group excel file, read group ids are not unique: {duplicates}"
        )
    read_groups = read_groups_df.to_dict(orient="index")
    clean_read_groups = {}
    for id_, read_group_info in read_groups.items():
        clean_read_groups[id_] = {
            key: value for key, value in read_group_info.items() if value
        }
    return clean_read_groups


def create_minimap_rg_str(read_id: str, read_group_info: dict, project_dir):
    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        force=True,
    )
    try:
        sample = read_group_info["sample"]
    except KeyError:
        excel_file = get_read_group_info_xls(project_dir)
        msg = f"The read group with id: {read_id} does not have a sample in excel file {excel_file}"
        logging.error(msg)
        raise RuntimeError(msg)

    rg_str = f"@RG\\\\tID:{read_id}\\\\tSM:{sample}"

    if read_group_info.get("platform", ""):
        rg_str += f"\\\\tPL:{read_group_info['platform']}"
    if read_group_info.get("library", ""):
        rg_str += f"\\\\tLB:{read_group_info['library']}"
    return rg_str


def get_read_group_id_from_path(path):
    return path.name.split(".")[0]


def get_read_groups_in_cram_header(cram_path: Path, project_dir) -> list[dict]:
    """The @RG lines of an alignment file, read from its own header.

    Only the header is read, so no reference genome is required and no
    alignment is decoded.
    """
    cmd = [SAMTOOLS_BIN, "view", "-H", str(cram_path)]
    process = run_cmd(cmd, project_dir=project_dir)["process"]

    read_groups = []
    for line in process.stdout.decode().splitlines():
        if not line.startswith("@RG\t"):
            continue
        read_group = {}
        for field in line.split("\t")[1:]:
            tag, sep, value = field.partition(":")
            if sep:
                read_group[tag] = value
        read_groups.append(read_group)
    return read_groups


def get_samples_in_cram(cram_path: Path, project_dir) -> set[str]:
    """The samples of a cram, taken from the SM tags of its @RG header lines.

    The mapping step writes the sample of every read group into the cram
    header, so the read group excel file is not needed to know which sample a
    cram holds.
    """
    read_groups = get_read_groups_in_cram_header(cram_path, project_dir)
    return {
        read_group["SM"] for read_group in read_groups if read_group.get("SM", "")
    }
