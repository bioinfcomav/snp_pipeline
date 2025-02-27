import tempfile
import shutil
from pathlib import Path

from .config import TEST_PROJECT2_DIR
from reads_pipeline.gatk import create_genome_reference


def test_create_genome_reference():
    with tempfile.TemporaryDirectory(prefix="gatk_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        genome_fasta = TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta"
        snv_calling_dir = project_dir_path / "snv_calling"
        snv_calling_dir.mkdir()

        create_genome_reference(genome_fasta, snv_calling_dir, project_dir_path)
