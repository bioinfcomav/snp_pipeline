import re
import stat
import tempfile
from pathlib import Path

import pytest

from reads_pipeline.paths import get_log_path
from reads_pipeline.run_cmd import run_cmd

# Every line of the log starts with the moment it was written
LOG_LINE = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} (INFO|ERROR|WARNING) ")

TALKATIVE_CMD = """#!/usr/bin/env python3
import sys
print("walking sample1")
print("something went to stderr", file=sys.stderr)
print("walking sample2")
"""


def _create_cmd(dir_path: Path, content=TALKATIVE_CMD) -> Path:
    path = dir_path / "talkative"
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_the_log_lines_say_when_they_were_written():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        cmd_path = _create_cmd(project_dir)

        run_cmd([str(cmd_path)], project_dir=project_dir)

        log_lines = get_log_path(project_dir).read_text().splitlines()
        assert log_lines
        assert all(LOG_LINE.match(line) for line in log_lines)
        assert "Running cmd:" in log_lines[0]


def test_a_streamed_command_writes_what_it_says_as_it_says_it():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        cmd_path = _create_cmd(project_dir)

        res = run_cmd([str(cmd_path)], project_dir=project_dir, stream_output=True)

        log = get_log_path(project_dir).read_text()
        # both streams are in the log, in the order in which they happened
        for said in ("walking sample1", "something went to stderr", "walking sample2"):
            assert said in log
        assert all(LOG_LINE.match(line) for line in log.splitlines())

        # and the caller still gets everything the command said
        stdout = res["process"].stdout.decode()
        assert "walking sample1" in stdout
        assert "walking sample2" in stdout


def test_a_command_that_is_not_streamed_is_not_in_the_log():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        cmd_path = _create_cmd(project_dir)

        res = run_cmd([str(cmd_path)], project_dir=project_dir)

        assert "walking sample1" not in get_log_path(project_dir).read_text()
        # the two streams are kept apart, which is what the commands that are
        # parsed, like samtools view -H, need
        assert res["process"].stdout.decode().startswith("walking sample1")
        assert res["process"].stderr.decode().startswith("something went to stderr")


@pytest.mark.parametrize("stream_output", [True, False])
def test_a_failed_command_is_logged(stream_output):
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        project_dir = Path(project_dir)
        cmd_path = _create_cmd(
            project_dir,
            content="#!/usr/bin/env python3\nprint('almost there')\nraise SystemExit(3)\n",
        )

        with pytest.raises(
            RuntimeError, match="There was an error running the command"
        ):
            run_cmd(
                [str(cmd_path)], project_dir=project_dir, stream_output=stream_output
            )

        log = get_log_path(project_dir).read_text()
        assert "There was a problem running" in log
        # whatever it said before it failed is in the log either way
        assert "almost there" in log
