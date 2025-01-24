import logging
from subprocess import run
from pathlib import Path
import tempfile

from .paths import get_log_path

logger = logging.getLogger(__name__)


def run_cmd(cmd, project_dir: Path):
    logging.basicConfig(
        filename=get_log_path(project_dir), filemode="a", level=logging.INFO, force=True
    )
    logging.info("Running cmd: " + " ".join(cmd))

    process = run(cmd, check=False, capture_output=True)
    if process.returncode:
        msg = f"There was a problem running: {cmd[0]}\n"
        msg += "cmd: " + " ".join(cmd) + "\n"
        msg += "stderr:\n" + process.stderr.decode()
        msg += "stdout:\n" + process.stdout.decode()
        logging.error(msg)
        raise RuntimeError("There was an error running the command: {msg}")

    return {"process": process}


def run_bash_script(script_content: str, project_dir: Path):
    logging.basicConfig(filename=get_log_path(project_dir), level=logging.INFO)

    with tempfile.NamedTemporaryFile(suffix=".sh", mode="wt") as shell_fhand:
        shell_fhand.write(script_content)
        shell_fhand.flush()
        cmd = ["bash", shell_fhand.name]
        process = run(cmd, check=False, capture_output=True)
        if process.returncode:
            msg = "There was a problem running a bash script\n"
            msg += "script\n{script_content}\n"
            msg += "stderr:\n" + process.stderr.decode()
            msg += "stdout:\n" + process.stdout.decode()
            logging.error(msg)
            raise RuntimeError("There was an error running a bash script")
