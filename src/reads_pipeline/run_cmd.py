import logging
import os
from subprocess import run, Popen, PIPE, STDOUT, CompletedProcess
from pathlib import Path
import tempfile

from .paths import get_log_path

logger = logging.getLogger(__name__)

# Every line of the log says when it was written, because working out what a
# step that ran for hours was doing, and when it stopped, is done with this file
LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(project_dir: Path):
    logging.basicConfig(
        filename=get_log_path(project_dir),
        filemode="a",
        level=logging.INFO,
        format=LOG_FORMAT,
        datefmt=LOG_DATE_FORMAT,
        force=True,
    )


def _run_streaming_output(cmd, env, verbose) -> CompletedProcess:
    """Runs a command writing everything it says, as it says it.

    A step that runs for hours says nothing at all if its output is kept until
    it ends, and it says nothing ever if it is killed, so the lines are written
    to the log, and to the terminal, as they arrive.  Both streams are merged
    into one, which is the order in which they happened.
    """
    lines = []
    with Popen(
        cmd, stdout=PIPE, stderr=STDOUT, env=env, bufsize=1, text=True
    ) as process:
        for line in process.stdout:
            line = line.rstrip("\n")
            lines.append(line)
            logging.info(f"{cmd[0]}: {line}")
            if verbose:
                print(line, flush=True)

    output = "\n".join(lines)
    return CompletedProcess(cmd, process.returncode, stdout=output.encode(), stderr=b"")


def run_cmd(
    cmd,
    project_dir: Path,
    verbose=False,
    env: None | dict = None,
    stream_output=False,
):
    setup_logging(project_dir)
    msg = "Running cmd: " + " ".join(cmd)
    if env:
        msg += " with env: " + " ".join(f"{key}={value}" for key, value in env.items())
    logging.info(msg)
    if verbose:
        print(f"Running: {cmd}")

    # The environment given is added to the one of this process, the command
    # still needs the PATH and everything else to run
    env = os.environ | env if env else None
    if stream_output:
        process = _run_streaming_output(cmd, env=env, verbose=verbose)
    else:
        process = run(cmd, check=False, capture_output=True, env=env)

    if process.returncode:
        msg = f"There was a problem running: {cmd[0]}\n"
        msg += "cmd: " + " ".join(cmd) + "\n"
        if not stream_output:
            # The output of a streamed command is in the log already
            msg += "stderr:\n" + process.stderr.decode()
            msg += "stdout:\n" + process.stdout.decode()
        if verbose:
            print(msg)
        logging.error(msg)
        raise RuntimeError(f"There was an error running the command: {msg}")

    return {"process": process}


def run_bash_script(script_content: str, project_dir: Path, verbose=False):
    setup_logging(project_dir)

    with tempfile.NamedTemporaryFile(suffix=".sh", mode="wt") as shell_fhand:
        shell_fhand.write(script_content)
        shell_fhand.flush()
        cmd = ["bash", shell_fhand.name]
        if verbose:
            print(f"Running: {cmd}")
        # input(f"{shell_fhand.name}")
        process = run(cmd, check=False, capture_output=True)
        if process.returncode:
            msg = "There was a problem running a bash script\n"
            msg += "script\n{script_content}\n"
            msg += "stderr:\n" + process.stderr.decode()
            msg += "stdout:\n" + process.stdout.decode()
            if verbose:
                print(msg)
            logging.error(msg)
            raise RuntimeError("There was an error running a bash script")
    return {"process": process}
