import tempfile
import shutil
from pathlib import Path
from subprocess import run

from .config import TEST_PROJECT6_DIR
from reads_pipeline.paths import (
    get_joint_var_calling_intervals_bed,
)


def test_add_sample_snv_calls_to_db():
    with tempfile.TemporaryDirectory(prefix="gatk_db_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT6_DIR, project_dir_path, dirs_exist_ok=True)

        cmd = [
            "uv",
            "run",
            "prepare_var_joining_segments",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)

        join_segments_bed = get_joint_var_calling_intervals_bed(project_dir)
        text = join_segments_bed.open("rt").read()
        assert text.startswith("SL4")
        assert text.endswith("2001\n")
