import shutil
import tempfile

from .config import (
    TEST_PROJECT2_DIR,
    MINIMAP_PROJECT2_TOMATO_INDEX,
    MINIMAP_PROJECT2_TOMATO_FASTA,
)
from reads_pipeline import (
    run_fastp_minimap,
    collect_cram_stats,
    plot_mapq_distributions,
)


def test_run_pipeline():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
            deduplicate=True,
        )
        collect_cram_stats(project_dir)
        plot_mapq_distributions(project_dir)

    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
            deduplicate=False,
        )
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
            deduplicate=False,
            re_run=False,
        )
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
            deduplicate=False,
            re_run=True,
        )
        collect_cram_stats(project_dir)
        plot_mapq_distributions(project_dir)
