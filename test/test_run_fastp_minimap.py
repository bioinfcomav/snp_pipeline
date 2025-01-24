import shutil
import tempfile

from .config import (
    TEST_PROJECT2_DIR,
    MINIMAP_PROJECT2_TOMATO_INDEX,
    MINIMAP_PROJECT2_TOMATO_FASTA,
)
from reads_pipeline import run_fastp_minimap


def test_run_pipeline():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
        )
