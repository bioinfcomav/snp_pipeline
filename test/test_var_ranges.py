import tempfile
import shutil
from pathlib import Path
from subprocess import run

from .config import TEST_PROJECT6_DIR


def test_var_ranges():
    with tempfile.TemporaryDirectory(prefix="gatk_db_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT6_DIR, project_dir_path, dirs_exist_ok=True)
        cmd = [
            "uv",
            "run",
            "merge_gvcf_var_ranges",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)
