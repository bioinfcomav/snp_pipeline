from pathlib import Path

from reads_pipeline import run_fastqc

TEST_DATA_DIR = Path(__file__).absolute().parent / "data"
TEST_PROJECT1_DIR = TEST_DATA_DIR / "project1"


def test_run_stats():
    run_fastqc(TEST_PROJECT1_DIR, re_run=True)
