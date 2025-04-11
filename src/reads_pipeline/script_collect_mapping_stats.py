import sys

from reads_pipeline.script_check_mapping import get_args
from reads_pipeline.paths import get_project_dir


def main():
    args = get_args()
    project_dir = get_project_dir(args.project_dir)

    if not project_dir.exists():
        msg = f"The project directorory {project_dir} does not exist"
        print(msg)
        sys.exit(2)

    print("hola")
