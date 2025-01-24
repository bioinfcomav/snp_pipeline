import tempfile
import shutil

from .config import TEST_PROJECT1_DIR
from reads_pipeline import (
    run_fastqc,
    collect_fastqc_stats,
    run_fastp,
    collect_fastp_stats,
)


def test_run_stats():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT1_DIR, project_dir, dirs_exist_ok=True)
        run_fastqc(project_dir, threads=2)

        run_fastqc(project_dir, re_run=False)
        run_fastqc(project_dir, re_run=True, threads=2)

        result = collect_fastqc_stats(project_dir)
        assert "num_seqs" in result["raw"].columns


def test_run_fastp():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT1_DIR, project_dir, dirs_exist_ok=True)
        run_fastp(project_dir)
        stats = collect_fastp_stats(project_dir)
        assert "dir" in stats.columns
        run_fastp(project_dir, re_run=True)
