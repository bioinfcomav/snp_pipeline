import os
import stat
import tempfile
from pathlib import Path

import pandas
import pytest

from reads_pipeline.paths import get_crams_dir, get_psp_path, get_psps_dir
from reads_pipeline.psp import (
    generate_psps_for_samples,
    get_cram_paths,
    get_samples_to_process,
    get_psp_paths,
)


def _write_fake_cram(cram_path: Path, read_group_id: str, *samples: str):
    """A fake cram: a sam header, which samtools reads like a cram's one"""
    header = "@HD\tVN:1.6\n"
    for idx, sample in enumerate(samples, start=1):
        rg_id = read_group_id if idx == 1 else f"{read_group_id}_{idx}"
        header += f"@RG\tID:{rg_id}\tSM:{sample}\tLB:lib1\n"
    cram_path.write_text(header)


def _create_project(project_dir: Path, read_groups: dict, create_catalog=True):
    """A project with the crams and the read group info already in place"""
    project_dir = Path(project_dir)
    reads_dir = project_dir / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    read_groups_df = pandas.DataFrame(
        {
            "id": list(read_groups.keys()),
            "sample": list(read_groups.values()),
            "library": ["lib1"] * len(read_groups),
        }
    )
    read_groups_df.to_excel(reads_dir / "reads.xlsx", index=False)

    crams_dir = get_crams_dir(project_dir) / "bioproject1"
    crams_dir.mkdir(parents=True, exist_ok=True)
    for read_group_id, sample in read_groups.items():
        _write_fake_cram(crams_dir / f"{read_group_id}.cram", read_group_id, sample)
    # the stats dir lives besides the bioproject dirs and holds no cram
    (get_crams_dir(project_dir) / "stats").mkdir(exist_ok=True)

    genome_fasta = project_dir / "genome.fasta"
    genome_fasta.write_text(">chrom1\nACGT\n")
    if create_catalog:
        Path(str(genome_fasta) + ".repeats.parquet").touch()
    return genome_fasta


def test_crams_are_grouped_per_sample():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _create_project(
            project_dir,
            {"rg1": "sample1", "rg2": "sample1", "rg3": "sample2"},
        )
        assert len(get_cram_paths(project_dir)) == 3

        samples_to_process = get_samples_to_process(project_dir)
        assert [info["sample"] for info in samples_to_process] == ["sample1", "sample2"]
        assert [info["idx"] for info in samples_to_process] == [1, 2]
        assert [path.name for path in samples_to_process[0]["cram_paths"]] == [
            "rg1.cram",
            "rg2.cram",
        ]
        assert [path.name for path in samples_to_process[1]["cram_paths"]] == [
            "rg3.cram"
        ]


def test_crams_with_no_sample_in_the_header_are_all_listed():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _create_project(project_dir, {"rg1": "sample1"})
        crams_dir = get_crams_dir(project_dir) / "bioproject1"
        # no @RG line at all, and an @RG line with no SM tag
        (crams_dir / "rg2.cram").write_text("@HD\tVN:1.6\n")
        (crams_dir / "rg3.cram").write_text("@HD\tVN:1.6\n@RG\tID:rg3\tLB:lib1\n")
        with pytest.raises(RuntimeError) as excinfo:
            get_samples_to_process(project_dir)
        msg = str(excinfo.value)
        # every offending cram is named, not just the first one
        assert "rg2.cram" in msg
        assert "rg3.cram" in msg
        assert "rg1.cram" not in msg
        assert "SM" in msg


def test_cram_with_several_samples():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _create_project(project_dir, {"rg1": "sample1"})
        crams_dir = get_crams_dir(project_dir) / "bioproject1"
        _write_fake_cram(crams_dir / "rg2.cram", "rg2", "sample2", "sample3")
        with pytest.raises(RuntimeError, match="several samples"):
            get_samples_to_process(project_dir)


def test_psps_do_not_need_the_read_group_excel():
    """The excel lives with the reads, which may not even be mounted"""
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _create_project(project_dir, {"rg1": "sample1", "rg2": "sample2"})
        (Path(project_dir) / "reads" / "reads.xlsx").unlink()

        samples_to_process = get_samples_to_process(project_dir)
        assert [info["sample"] for info in samples_to_process] == ["sample1", "sample2"]


def test_samples_already_done_are_skipped():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        genome_fasta = _create_project(
            project_dir, {"rg1": "sample1", "rg2": "sample2"}
        )
        # psps generated in a previous run
        get_psps_dir(project_dir).mkdir(parents=True, exist_ok=True)
        get_psp_path(project_dir, "sample1").write_text("psp1")
        get_psp_path(project_dir, "sample2").write_text("psp2")

        # pop_var_caller is not run at all, so this works with no binary installed
        res = generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False
        )
        assert res["num_analyses_done"] == 0
        assert not res["psp_paths"]
        assert [path.name for path in get_psp_paths(project_dir)] == [
            "sample1.psp",
            "sample2.psp",
        ]
        # the psps of the previous run are left untouched
        assert get_psp_path(project_dir, "sample1").read_text() == "psp1"


def test_missing_repeat_catalog():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        genome_fasta = _create_project(
            project_dir, {"rg1": "sample1"}, create_catalog=False
        )
        with pytest.raises(FileNotFoundError, match="repeat-catalog"):
            generate_psps_for_samples(
                project_dir, genome_fasta=genome_fasta, verbose=False
            )


# It writes one psp per sample, just like pop_var_caller generate-psps does, but
# it reads the sample from the fake cram instead of walking the alignments
FAKE_POP_VAR_CALLER = """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
assert args[0] == "generate-psps"
alignments = [Path(args[idx + 1]) for idx, arg in enumerate(args) if arg == "--alignment"]
output_dir = Path(args[args.index("--output-dir") + 1])
assert Path(args[args.index("--reference") + 1]).exists()
assert Path(args[args.index("--catalog") + 1]).exists()

samples = set()
for path in alignments:
    for line in path.read_text().splitlines():
        if line.startswith("@RG\t"):
            samples.update(f[3:] for f in line.split("\t") if f.startswith("SM:"))
assert len(samples) == 1, samples
sample = samples.pop()
psp_path = output_dir / (sample + ".psp")
psp_path.write_text(",".join(sorted(path.name for path in alignments)))
(Path(os.environ["RAYON_REPORT_DIR"]) / (sample + ".rayon")).write_text(
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
    monkeypatch.setenv("RAYON_REPORT_DIR", str(bin_dir))


def test_generate_psps(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(
            project_dir,
            {"rg1": "sample1", "rg2": "sample1", "rg3": "sample2"},
        )
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)

        res = generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False
        )
        assert res["num_analyses_done"] == 2
        assert [path.name for path in get_psp_paths(project_dir)] == [
            "sample1.psp",
            "sample2.psp",
        ]
        # the crams of one sample are walked together into one psp
        assert get_psp_path(project_dir, "sample1").read_text() == "rg1.cram,rg2.cram"
        # no tmp dir is left behind
        assert sorted(path.name for path in get_psps_dir(project_dir).iterdir()) == [
            "sample1.psp",
            "sample2.psp",
        ]

        # a second run has nothing to do
        res = generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False
        )
        assert res["num_analyses_done"] == 0

        # a new sample is generated, the ones already done are kept
        _write_fake_cram(
            get_crams_dir(project_dir) / "bioproject1" / "rg4.cram", "rg4", "sample3"
        )
        read_groups_df = pandas.DataFrame(
            {
                "id": ["rg1", "rg2", "rg3", "rg4"],
                "sample": ["sample1", "sample1", "sample2", "sample3"],
                "library": ["lib1"] * 4,
            }
        )
        read_groups_df.to_excel(project_dir / "reads" / "reads.xlsx", index=False)
        res = generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False
        )
        assert res["num_analyses_done"] == 1
        assert [path.name for path in res["psp_paths"]] == ["sample3.psp"]

        res = generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False, re_run=True
        )
        assert res["num_analyses_done"] == 3


def test_generate_psps_in_parallel(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(
            project_dir, {"rg1": "sample1", "rg2": "sample2", "rg3": "sample3"}
        )
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)

        res = generate_psps_for_samples(
            project_dir,
            genome_fasta=genome_fasta,
            verbose=False,
            num_psps_in_parallel=3,
        )
        assert res["num_analyses_done"] == 3
        assert len(get_psp_paths(project_dir)) == 3


def test_a_failed_run_leaves_no_psp(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, {"rg1": "sample1"})
        _create_fake_pop_var_caller(
            project_dir / "bin",
            monkeypatch,
            content="#!/usr/bin/env bash\nexit 1\n",
        )

        with pytest.raises(RuntimeError):
            generate_psps_for_samples(
                project_dir, genome_fasta=genome_fasta, verbose=False
            )
        # neither a psp nor a partial file is left behind
        assert not list(get_psps_dir(project_dir).iterdir())


def test_num_threads_bounds_the_walk(monkeypatch):
    """pop_var_caller generate-psps takes no --threads option, so the only way
    of narrowing every walk is the size of the rayon global thread pool"""
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, {"rg1": "sample1"})
        _create_fake_pop_var_caller(project_dir / "bin", monkeypatch)

        generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False, num_threads=3
        )
        rayon_report = project_dir / "bin" / "sample1.rayon"
        assert rayon_report.read_text() == "3"

        # zero threads means every core, so nothing is said about the pool
        monkeypatch.delenv("RAYON_NUM_THREADS", raising=False)
        generate_psps_for_samples(
            project_dir, genome_fasta=genome_fasta, verbose=False, re_run=True
        )
        assert rayon_report.read_text() == "unset"
