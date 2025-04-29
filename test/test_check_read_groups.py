import tempfile
from pathlib import Path
import shutil
from subprocess import run

from reads_pipeline.script_check_read_groups_to_map import check_read_groups_to_map

from .config import TEST_PROJECT7_DIR


def _prepare_snv_calling_test_dir(project_dir):
    project_dir_path = Path(project_dir)
    shutil.copytree(TEST_PROJECT7_DIR, project_dir_path, dirs_exist_ok=True)


def test_check_read_groups():
    with tempfile.TemporaryDirectory(prefix="check_read_groups_test") as project_dir:
        _prepare_snv_calling_test_dir(project_dir)
        project_dir_path = Path(project_dir)

        script_path = "check_read_groups_to_map"
        cmd = ["uv", "run", str(script_path), str(project_dir_path)]
        process = run(cmd, check=True, cwd=project_dir, capture_output=True)
        assert b"ERROR" in process.stdout
