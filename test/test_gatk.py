import tempfile
import shutil
from pathlib import Path
from subprocess import run

from .config import TEST_PROJECT2_DIR, TEST_PROJECT5_DIR
from reads_pipeline.gatk import (
    create_genome_reference,
    do_sample_snv_calling_basic_germline,
    create_db_with_independent_sample_snv_calls,
    get_samples_in_gatk_db,
    GATKDBFileMode,
    do_svn_joint_genotyping_for_all_samples_together,
    do_snv_calling_per_sample,
)
from reads_pipeline.fastp_minimap import run_fastp_minimap_for_fastqs

from reads_pipeline.paths import MINIMAP2_BIN, get_vcfs_per_sample_dir


def test_snv_calling_per_sample():
    with tempfile.TemporaryDirectory(prefix="gatk_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT5_DIR, project_dir, dirs_exist_ok=True)
        genome_fasta = TEST_PROJECT5_DIR / "SL4.0ch01_13486265-13488265.fasta"
        snv_calling_dir = project_dir_path / "snv_calling"
        snv_calling_dir.mkdir()

        res = create_genome_reference(genome_fasta, snv_calling_dir, project_dir_path)
        genome_reference_path = res["genome_path"]

        minimap_index = genome_reference_path.with_suffix(".mmi")
        cmd = [MINIMAP2_BIN, "-d", str(minimap_index), str(genome_reference_path)]
        run(cmd, capture_output=True, check=True)

        res = run_fastp_minimap_for_fastqs(
            project_dir=project_dir_path,
            minimap_index=minimap_index,
            genome_fasta=genome_reference_path,
            deduplicate=False,
        )

        res = do_snv_calling_per_sample(
            project_dir=project_dir_path, genome_fasta=genome_fasta
        )
        assert res["out_vcf_paths_done"]["sample1"].exists()
        assert len(res["out_vcf_paths_done"]) == 1

        res = do_snv_calling_per_sample(
            project_dir=project_dir_path, genome_fasta=genome_fasta
        )
        assert len(res["out_vcf_paths_done"]) == 0

        res = do_snv_calling_per_sample(
            project_dir=project_dir_path,
            genome_fasta=genome_fasta,
            re_run=True,
        )
        assert len(res["out_vcf_paths_done"]) == 1


def test_create_genome_reference():
    with tempfile.TemporaryDirectory(prefix="gatk_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT2_DIR, project_dir, dirs_exist_ok=True)
        genome_fasta = TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta"
        snv_calling_dir = project_dir_path / "snv_calling"
        snv_calling_dir.mkdir()

        res = create_genome_reference(genome_fasta, snv_calling_dir, project_dir_path)
        genome_reference_path = res["genome_path"]
        fai_path = res["fai_path"]

        minimap_index = genome_reference_path.with_suffix(".mmi")
        cmd = [MINIMAP2_BIN, "-d", str(minimap_index), str(genome_reference_path)]
        run(cmd, capture_output=True, check=True)

        res = run_fastp_minimap(
            project_dir=project_dir_path,
            minimap_index=minimap_index,
            genome_fasta=genome_reference_path,
            deduplicate=False,
        )
        cram_paths = res["cram_paths"]

        vcf_bams_dir = get_vcfs_for_bams_dir(project_dir)
        vcf_bams_dir.mkdir(parents=True)
        vcf_paths = []
        for cram_path in cram_paths:
            vcf_path = vcf_bams_dir / cram_path.with_suffix(".vcf.gz").name
            do_sample_snv_calling_basic_germline(
                genome_reference_path,
                [cram_path],
                out_vcf=vcf_path,
                project_dir=project_dir,
            )
            vcf_paths.append(vcf_path)

        create_db_with_independent_sample_snv_calls(
            vcf_paths,
            project_dir=project_dir,
            genome_fai_path=fai_path,
            mode=GATKDBFileMode.CREATE,
        )
        samples = get_samples_in_gatk_db(project_dir=project_dir)
        assert samples == ["sample1"]

        out_vcf = snv_calling_dir / "joint_genotyping.vcf.gz"
        do_svn_joint_genotyping_for_all_samples_together(
            project_dir=project_dir, genome_fasta=genome_reference_path, out_vcf=out_vcf
        )
