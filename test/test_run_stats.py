from pathlib import Path
import tempfile
import shutil

from reads_pipeline import run_fastqc

TEST_DATA_DIR = Path(__file__).absolute().parent / "data"
TEST_PROJECT1_DIR = TEST_DATA_DIR / "project1"


def test_run_stats():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT1_DIR, project_dir, dirs_exist_ok=True)
        run_fastqc(project_dir)

        run_fastqc(project_dir, re_run=False)
        run_fastqc(project_dir, re_run=True)
