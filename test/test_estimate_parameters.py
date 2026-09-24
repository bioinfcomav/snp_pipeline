import os
import shutil
import stat
import tempfile
import tomllib
from pathlib import Path

import pandas
import pytest

from reads_pipeline.paths import get_crams_dir, get_parameters_path, get_psps_dir
from reads_pipeline.parameter_estimation import (
    declare_sequencing_batches,
    estimate_parameters,
    get_batch_of_each_read_group,
)

from .config import TEST_DATA_DIR

PARAMETERS_AS_WRITTEN = TEST_DATA_DIR / "parameters" / "parameters_as_written.toml"
PARAMETERS_AS_SERDE_WRITES_IT = (
    TEST_DATA_DIR / "parameters" / "parameters_as_serde_writes_it.toml"
)
# The read groups and the samples of both parameters file fixtures
GOLDEN_SAMPLE2 = 'Ailsa ‘Craig’ "×2"'


def _create_project(project_dir: Path, read_groups: dict, create_psps=True):
    """A project with the crams, the psps and the read group info in place.

    read_groups maps every read group id to (sample, project dir name)
    """
    project_dir = Path(project_dir)
    reads_dir = project_dir / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    pandas.DataFrame(
        {
            "id": list(read_groups.keys()),
            "sample": [sample for sample, _ in read_groups.values()],
            "library": ["lib1"] * len(read_groups),
        }
    ).to_excel(reads_dir / "reads.xlsx", index=False)

    for read_group_id, (sample, bioproject) in read_groups.items():
        crams_dir = get_crams_dir(project_dir) / bioproject
        crams_dir.mkdir(parents=True, exist_ok=True)
        (crams_dir / f"{read_group_id}.cram").write_text(
            f"@HD\tVN:1.6\n@RG\tID:{read_group_id}\tSM:{sample}\tLB:lib1\n"
        )

    if create_psps:
        psps_dir = get_psps_dir(project_dir)
        psps_dir.mkdir(parents=True, exist_ok=True)
        for sample, _ in read_groups.values():
            (psps_dir / f"{sample}.psp").write_text("psp")

    genome_fasta = project_dir / "genome.fasta"
    genome_fasta.write_text(">chrom1\nACGT\n")
    Path(str(genome_fasta) + ".repeats.parquet").touch()
    return genome_fasta


def test_batch_of_each_read_group():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _create_project(
            project_dir,
            {
                "rg1": ("sample1", "bioproject1"),
                "rg2": ("sample1", "bioproject1"),
                "rg3": ("sample2", "bioproject2"),
            },
        )
        assert get_batch_of_each_read_group(project_dir) == {
            "rg1": "bioproject1",
            "rg2": "bioproject1",
            "rg3": "bioproject2",
        }


@pytest.mark.parametrize(
    "fixture", [PARAMETERS_AS_WRITTEN, PARAMETERS_AS_SERDE_WRITES_IT]
)
def test_declare_sequencing_batches(fixture):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as tmp_dir:
        parameters_path = Path(tmp_dir) / "parameters.toml"
        shutil.copy(fixture, parameters_path)
        before = tomllib.loads(parameters_path.read_text())

        # the two read groups of TS-1 were sequenced in one project and the
        # read group of the other sample in another one
        res = declare_sequencing_batches(
            parameters_path,
            {"HWI.3": "bioproject2", "HWI.4": "bioproject2", "HWI.5": "bioproject1"},
        )
        assert res == {"num_batches": 2, "num_read_groups": 3, "num_samples": 2}

        after = tomllib.loads(parameters_path.read_text())
        batches = after["sequencing_batches"]
        assert batches["batching_was_declared"]
        # the batches are numbered by project name, so bioproject1 is batch 0
        assert batches["by_read_group"] == [
            {"read_group": 0, "batch": 1},
            {"read_group": 1, "batch": 1},
            {"read_group": 2, "batch": 0},
        ]
        assert batches["by_sample"] == [
            {"sample": "TS-1", "batch": 1},
            {"sample": GOLDEN_SAMPLE2, "batch": 0},
        ]

        # nothing else in the file is touched
        for section in before:
            if section == "sequencing_batches":
                continue
            assert after[section] == before[section]


def test_a_sample_sequenced_in_two_projects():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as tmp_dir:
        parameters_path = Path(tmp_dir) / "parameters.toml"
        shutil.copy(PARAMETERS_AS_WRITTEN, parameters_path)
        text_before = parameters_path.read_text()

        # the two read groups of TS-1 are in different projects
        with pytest.raises(RuntimeError, match="only be in one sequencing batch"):
            declare_sequencing_batches(
                parameters_path,
                {
                    "HWI.3": "bioproject1",
                    "HWI.4": "bioproject2",
                    "HWI.5": "bioproject2",
                },
            )
        assert parameters_path.read_text() == text_before


def test_read_group_with_no_cram():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as tmp_dir:
        parameters_path = Path(tmp_dir) / "parameters.toml"
        shutil.copy(PARAMETERS_AS_WRITTEN, parameters_path)
        with pytest.raises(RuntimeError, match="sequencing batch is unknown"):
            declare_sequencing_batches(
                parameters_path, {"HWI.3": "bioproject1", "HWI.4": "bioproject1"}
            )


FAKE_PARAMETERS = """format_version = 1
ploidy = 2

[fitted_from]
reference_digest = "0123456789abcdef0123456789abcdef"
samples = [
    "sample1",
    "sample2",
]
read_groups = [
    { read_group = 0, declared_id = "rg1", library = "lib1", sample = "sample1" },
    { read_group = 1, declared_id = "rg2", library = "lib1", sample = "sample1" },
    { read_group = 2, declared_id = "rg3", library = "lib1", sample = "sample2" },
]

[sequencing_batches]
batching_was_declared = false
by_read_group = [
    { read_group = 0, batch = 0 },
    { read_group = 1, batch = 0 },
    { read_group = 2, batch = 0 },
]
by_sample = [
    { sample = "sample1", batch = 0 },
    { sample = "sample2", batch = 0 },
]

[ordinary_site_prior]
reference_concentration = 1.0
"""

# It writes the parameters file that pop_var_caller estimate-parameters would
# write, and checks the arguments it is given
FAKE_POP_VAR_CALLER = """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
assert args[0] == "estimate-parameters"
assert Path(args[args.index("--reference") + 1]).exists()
assert Path(args[args.index("--catalog") + 1]).exists()
psps_dir = Path(args[args.index("--psp") + 1])
assert list(psps_dir.glob("*.psp"))
assert args[args.index("--ploidy") + 1] == os.environ["EXPECTED_PLOIDY"]
Path(args[args.index("--output") + 1]).write_text(os.environ["FAKE_PARAMETERS"])
Path(os.environ["RAYON_REPORT"]).write_text(
    os.environ.get("RAYON_NUM_THREADS", "unset")
)
"""


def _create_fake_pop_var_caller(
    bin_dir: Path, monkeypatch, content=FAKE_POP_VAR_CALLER
):
    bin_dir.mkdir(parents=True, exist_ok=True)
    fake_bin = bin_dir / "pop_var_caller"
    fake_bin.write_text(content)
    fake_bin.chmod(fake_bin.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_PARAMETERS", FAKE_PARAMETERS)
    monkeypatch.setenv("EXPECTED_PLOIDY", "2")
    monkeypatch.setenv("RAYON_REPORT", str(bin_dir.parent / "rayon_num_threads"))


def test_estimate_parameters(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(
            project_dir,
            {
                "rg1": ("sample1", "bioproject1"),
                "rg2": ("sample1", "bioproject1"),
                "rg3": ("sample2", "bioproject2"),
            },
        )
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)

        res = estimate_parameters(project_dir, genome_fasta=genome_fasta, verbose=False)
        assert res["should_have_run"]
        parameters_path = get_parameters_path(project_dir)
        assert res["parameters_path"] == parameters_path

        batches = tomllib.loads(parameters_path.read_text())["sequencing_batches"]
        assert batches["batching_was_declared"]
        assert batches["by_read_group"] == [
            {"read_group": 0, "batch": 0},
            {"read_group": 1, "batch": 0},
            {"read_group": 2, "batch": 1},
        ]
        assert batches["by_sample"] == [
            {"sample": "sample1", "batch": 0},
            {"sample": "sample2", "batch": 1},
        ]

        # no tmp dir is left behind
        assert [
            path.name for path in parameters_path.parent.iterdir() if path.is_dir()
        ] == ["psps"]

        # a second run does nothing
        res = estimate_parameters(project_dir, genome_fasta=genome_fasta, verbose=False)
        assert not res["should_have_run"]

        res = estimate_parameters(
            project_dir, genome_fasta=genome_fasta, verbose=False, re_run=True
        )
        assert res["should_have_run"]


def test_estimate_parameters_with_no_batches(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, {"rg1": ("sample1", "bioproject1")})
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)

        estimate_parameters(
            project_dir,
            genome_fasta=genome_fasta,
            verbose=False,
            declare_batches=False,
        )
        batches = tomllib.loads(get_parameters_path(project_dir).read_text())[
            "sequencing_batches"
        ]
        assert not batches["batching_was_declared"]


def test_no_psps_to_fit():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(
            project_dir, {"rg1": ("sample1", "bioproject1")}, create_psps=False
        )
        with pytest.raises(RuntimeError, match="no psp file"):
            estimate_parameters(project_dir, genome_fasta=genome_fasta, verbose=False)


def test_a_failed_run_leaves_no_parameters(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, {"rg1": ("sample1", "bioproject1")})
        _create_fake_pop_var_caller(
            project_dir / "bin",
            monkeypatch,
            content="#!/usr/bin/env bash\nexit 1\n",
        )
        with pytest.raises(RuntimeError):
            estimate_parameters(project_dir, genome_fasta=genome_fasta, verbose=False)
        assert not get_parameters_path(project_dir).exists()


def test_num_threads_bounds_the_fit(monkeypatch):
    """pop_var_caller estimate-parameters takes no --threads option, so the only
    way of narrowing the fit is the size of the rayon global thread pool"""
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(
            project_dir,
            {
                "rg1": ("sample1", "bioproject1"),
                "rg2": ("sample1", "bioproject1"),
                "rg3": ("sample2", "bioproject2"),
            },
        )
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)
        rayon_report = project_dir / "rayon_num_threads"

        estimate_parameters(
            project_dir, genome_fasta=genome_fasta, verbose=False, num_threads=4
        )
        assert rayon_report.read_text() == "4"

        # zero threads means every core, so nothing is said about the pool
        monkeypatch.delenv("RAYON_NUM_THREADS", raising=False)
        estimate_parameters(
            project_dir, genome_fasta=genome_fasta, verbose=False, re_run=True
        )
        assert rayon_report.read_text() == "unset"


def test_a_read_group_in_two_projects():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        _create_project(project_dir, {"rg1": ("sample1", "bioproject1")})
        # the same read group, in the crams of another project
        crams_dir = get_crams_dir(project_dir) / "bioproject2"
        crams_dir.mkdir(parents=True, exist_ok=True)
        (crams_dir / "rg1.cram").write_text(
            "@HD\tVN:1.6\n@RG\tID:rg1\tSM:sample1\tLB:lib1\n"
        )
        with pytest.raises(RuntimeError, match="crams of two projects"):
            get_batch_of_each_read_group(project_dir)
