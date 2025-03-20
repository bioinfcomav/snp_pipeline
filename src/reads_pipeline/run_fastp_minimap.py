import os
import argparse

from reads_pipeline import run_fastp_minimap


def get_args():
    parser = argparse.ArgumentParser(
        description="Cleans and aligns reads using fastp and minimap2."
    )

    # project_dir: Defaults to the current working directory
    parser.add_argument(
        "project_dir",
        nargs="?",  # Makes it optional
        default=os.getcwd(),
        help="Project directory (default: current working directory).",
    )

    # minimap_index: Required path argument
    parser.add_argument(
        "--minimap_index",
        type=str,
        required=True,
        help="Path to the minimap2 index file.",
    )

    # genome_fasta: Required path argument
    parser.add_argument(
        "--genome_fasta", type=str, required=True, help="Path to the genome FASTA file."
    )

    # deduplicate: Required boolean argument
    parser.add_argument(
        "--deduplicate",
        type=lambda x: (str(x).lower() in ["true", "1", "yes"]),
        required=True,
        help="Boolean flag to enable or disable deduplication (true/false).",
    )

    return parser.parse_args()


def main():
    args = get_args()
    run_fastp_minimap(
        project_dir=args.project_dir,
        genome_fasta=args.genome_fasta,
        minimap_index=args.minimap_index,
        deduplicate=args.deduplicate,
    )


if __name__ == "__main__":
    main()
