import shutil
import tempfile
from subprocess import run

from .config import (
    TEST_PROJECT2_DIR,
    MINIMAP_PROJECT2_TOMATO_INDEX,
    MINIMAP_PROJECT2_TOMATO_FASTA,
    TEST_PROJECT3_DIR,
    MINIMAP_PROJECT3_GENOME_INDEX,
    MINIMAP_PROJECT3_GENOME_FASTA,
    TEST_PROJECT4_DIR,
    MINIMAP_PROJECT4_GENOME_INDEX,
    MINIMAP_PROJECT4_GENOME_FASTA,
)
from reads_pipeline import (
    run_fastp_minimap,
    collect_cram_stats,
    plot_mapq_distributions,
    plot_coverage_distributions,
)


def test_run_pipeline():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT2_TOMATO_INDEX,
            genome_fasta=MINIMAP_PROJECT2_TOMATO_FASTA,
            deduplicate=True,
            fastp_trim_front1=1,
            fastp_trim_front2=1,
            fastp_trim_tail1=1,
            fastp_trim_tail2=2,
        )
        plot_coverage_distributions(project_dir)
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


def test_minimap_ref_path():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT3_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT3_GENOME_INDEX,
            genome_fasta=MINIMAP_PROJECT3_GENOME_FASTA,
            deduplicate=False,
            fastp_trim_front1=10,
            fastp_trim_front2=10,
            fastp_trim_tail1=10,
            fastp_trim_tail2=10,
        )

        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT3_GENOME_INDEX,
            genome_fasta=MINIMAP_PROJECT3_GENOME_FASTA,
            deduplicate=False,
            fastp_trim_front1=10,
            fastp_trim_front2=10,
            fastp_trim_tail1=10,
            fastp_trim_tail2=10,
            trim_quals_num_bases=0,
            re_run=True,
        )


def test_trim_quals_order():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        shutil.copytree(TEST_PROJECT4_DIR, project_dir, dirs_exist_ok=True)
        run_fastp_minimap(
            project_dir,
            minimap_index=MINIMAP_PROJECT4_GENOME_INDEX,
            genome_fasta=MINIMAP_PROJECT4_GENOME_FASTA,
            deduplicate=False,
            fastp_trim_front1=10,
            fastp_trim_front2=10,
            fastp_trim_tail1=10,
            fastp_trim_tail2=10,
        )


def test_run_mapping():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        print("testing")
        shutil.copytree(TEST_PROJECT4_DIR, project_dir, dirs_exist_ok=True)
        cmd = [
            "uv",
            "run",
            "run_mapping",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)
