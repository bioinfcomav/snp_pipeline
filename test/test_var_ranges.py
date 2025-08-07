import tempfile
import shutil
from pathlib import Path
from subprocess import run

from .config import TEST_PROJECT6_DIR
from reads_pipeline.paths import get_joint_gatk_segments_bed


def test_var_ranges():
    with tempfile.TemporaryDirectory(prefix="gatk_db_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT6_DIR, project_dir_path, dirs_exist_ok=True)
        cmd = [
            "uv",
            "run",
            "prepare_gatk_sample_var_joining_segments",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)
        bed_path = get_joint_gatk_segments_bed(project_dir)
        line = bed_path.open("rt").readline()
        assert line == "SL4.0ch01\t1\t4\n"
