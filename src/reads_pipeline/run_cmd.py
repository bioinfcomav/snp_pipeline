import logging
from subprocess import run

from .paths import get_log_path

logger = logging.getLogger(__name__)


def run_cmd(cmd, project_dir):
    logging.basicConfig(filename=get_log_path(project_dir), level=logging.INFO)
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
