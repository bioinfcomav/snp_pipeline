import os
import shutil


def move_files_and_dirs(src_dir: str, dest_dir: str):
    if not os.path.isdir(src_dir):
        raise ValueError(f"Source directory does not exist: {src_dir}")
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    for filename in os.listdir(src_dir):
        src_file = os.path.join(src_dir, filename)
        dest_file = os.path.join(dest_dir, filename)

        shutil.move(src_file, dest_file)
