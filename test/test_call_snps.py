import os
import stat
import tempfile
from pathlib import Path

import pytest

from reads_pipeline.paths import (
    get_parameters_path,
    get_psps_dir,
    get_snp_calling_dir,
    get_vcf_path,
)
from reads_pipeline.snp_calling import call_snps, get_vcf_index_path


def _create_project(project_dir: Path, samples=("sample1", "sample2"), parameters=True):
    """A project with the psps, and their parameters, already fitted"""
    project_dir = Path(project_dir)
    psps_dir = get_psps_dir(project_dir)
    psps_dir.mkdir(parents=True, exist_ok=True)
    for sample in samples:
        (psps_dir / f"{sample}.psp").write_text("psp")

    if parameters:
        get_parameters_path(project_dir).write_text("format_version = 1\n")

    genome_fasta = project_dir / "genome.fasta"
    genome_fasta.write_text(">chrom1\nACGT\n")
    Path(str(genome_fasta) + ".repeats.parquet").touch()
    return genome_fasta


# It writes the vcf, the parameters it called with, which pop_var_caller always
# writes beside the vcf, and the temporary files it leaves behind when it stops
FAKE_POP_VAR_CALLER = """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
assert args[0] == "call-from-psps"
assert Path(args[args.index("--reference") + 1]).exists()
assert Path(args[args.index("--catalog") + 1]).exists()
psps_dir = Path(args[args.index("--psp") + 1])
psps = sorted(path.stem for path in psps_dir.glob("*.psp"))
assert psps

output = Path(args[args.index("--output") + 1])
Path(os.environ["CMD_REPORT"]).write_text("\\n".join(args))

output.write_text("##fileformat=VCFv4.3\\n#CHROM\\t" + "\\t".join(psps) + "\\n")
stem = output.name.split(".")[0]
(output.parent / (stem + ".parameters.toml")).write_text("format_version = 1\\n")
(output.parent / (output.name + ".paralog-spill.tmp")).write_text("spill")
"""

FAKE_TABIX = """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
assert args[:2] == ["-p", "vcf"]
vcf_path = Path(args[2])
assert vcf_path.exists()
Path(str(vcf_path) + ".tbi").write_text("index")
"""


def _create_fake_bins(bin_dir: Path, monkeypatch, caller=FAKE_POP_VAR_CALLER):
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name, content in (("pop_var_caller", caller), ("tabix", FAKE_TABIX)):
        fake_bin = bin_dir / name
        fake_bin.write_text(content)
        fake_bin.chmod(fake_bin.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("CMD_REPORT", str(bin_dir / "cmd"))


def _get_cmd(project_dir: Path) -> list[str]:
    return (project_dir / "bin" / "cmd").read_text().splitlines()


def test_call_snps(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir)
        _create_fake_bins(project_dir / "bin", monkeypatch)

        res = call_snps(project_dir, genome_fasta=genome_fasta, verbose=False)
        assert res["should_have_run"]

        vcf_path = get_vcf_path(project_dir)
        assert res["vcf_path"] == vcf_path
        assert vcf_path.read_text().splitlines()[1].endswith("sample1\tsample2")
        # the vcf is indexed and the parameters it was called with are kept
        assert get_vcf_index_path(vcf_path).exists()
        assert (get_snp_calling_dir(project_dir) / "variants.parameters.toml").exists()

        # the numbers come from the fitted parameters, not from the defaults
        cmd = _get_cmd(project_dir)
        assert "--parameters" in cmd
        assert cmd[cmd.index("--parameters") + 1] == str(
            get_parameters_path(project_dir).resolve()
        )
        assert "--defaults" not in cmd

        # neither the temporary dir nor the files that the calling left in it
        # are kept
        assert sorted(
            path.name for path in get_snp_calling_dir(project_dir).iterdir()
        ) == [
            "parameters.toml",
            "psps",
            "variants.parameters.toml",
            "variants.vcf.gz",
            "variants.vcf.gz.tbi",
        ]

        # a second run has nothing to do
        res = call_snps(project_dir, genome_fasta=genome_fasta, verbose=False)
        assert not res["should_have_run"]

        res = call_snps(
            project_dir, genome_fasta=genome_fasta, verbose=False, re_run=True
        )
        assert res["should_have_run"]


def test_the_repeat_criteria_and_the_threads_are_given(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir)
        _create_fake_bins(project_dir / "bin", monkeypatch)

        call_snps(
            project_dir,
            genome_fasta=genome_fasta,
            verbose=False,
            min_copies="8,6,6,6,5,4",
            min_period=2,
            max_period=5,
            max_str_len=90,
            min_purity=0.9,
            num_threads=6,
            paralog_fdr=0,
            ploidy=4,
        )
        cmd = _get_cmd(project_dir)
        for option, value in (
            ("--min-copies", "8,6,6,6,5,4"),
            ("--min-period", "2"),
            ("--max-period", "5"),
            ("--max-str-len", "90"),
            ("--min-purity", "0.9"),
            ("--threads", "6"),
            ("--paralog-fdr", "0"),
            ("--ploidy", "4"),
        ):
            assert cmd[cmd.index(option) + 1] == value


def test_calling_with_the_defaults(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, parameters=False)
        _create_fake_bins(project_dir / "bin", monkeypatch)

        # the defaults are a different claim, so they have to be asked for
        with pytest.raises(FileNotFoundError, match="estimate_parameters"):
            call_snps(project_dir, genome_fasta=genome_fasta, verbose=False)

        call_snps(
            project_dir,
            genome_fasta=genome_fasta,
            verbose=False,
            use_default_parameters=True,
        )
        cmd = _get_cmd(project_dir)
        assert "--defaults" in cmd
        assert "--parameters" not in cmd


def test_no_psps_to_call():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir, samples=())
        with pytest.raises(RuntimeError, match="no psp file"):
            call_snps(project_dir, genome_fasta=genome_fasta, verbose=False)


def test_a_failed_run_leaves_no_vcf(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        genome_fasta = _create_project(project_dir)
        _create_fake_bins(
            project_dir / "bin", monkeypatch, caller="#!/usr/bin/env bash\nexit 1\n"
        )
        with pytest.raises(RuntimeError):
            call_snps(project_dir, genome_fasta=genome_fasta, verbose=False)

        assert not get_vcf_path(project_dir).exists()
        # not even the temporary dir the calling was writing into
        assert sorted(
            path.name for path in get_snp_calling_dir(project_dir).iterdir()
        ) == ["parameters.toml", "psps"]
