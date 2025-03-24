import tempfile
import shutil
from pathlib import Path
from subprocess import run

from .config import (
    TEST_PROJECT2_DIR,
)


def test_run_pipeline():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        project_path = Path(project_dir)
        script_path = "run_mapping"
        cmd = ["uv", "run", str(script_path), str(project_path)]
        run(cmd, check=True, cwd=project_dir)
