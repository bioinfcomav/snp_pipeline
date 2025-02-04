import tempfile
import shutil
from pathlib import Path
from subprocess import run

import reads_pipeline
from .config import (
    TEST_PROJECT2_DIR,
)


def test_run_pipeline():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        config_path = Path(project_dir) / "pipeline.toml"
        script_path = Path(reads_pipeline.__path__[0]) / "genomic_pipeline.py"
        cmd = ["uv", "run", str(script_path), str(config_path)]
        run(cmd, check=True, cwd=project_dir)
