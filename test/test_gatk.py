import tempfile
import shutil
from pathlib import Path
from subprocess import run


from .config import (
    TEST_PROJECT2_DIR,
    TEST_PROJECT5_DIR,
    TEST_PROJECT6_DIR,
    PROJECT6_GENOME_FAI,
    PROJECT6_GENOME_FASTA,
)
from reads_pipeline.gatk import (
    create_genome_reference,
    do_sample_snv_calling_basic_germline,
    create_db_with_independent_sample_snv_calls,
    get_samples_in_gatk_db,
    GATKDBFileMode,
    do_svn_joint_genotyping_for_all_samples_together,
    do_snv_calling_per_sample,
    filter_vcf_with_gatk,
)
from reads_pipeline.fastp_minimap import run_fastp_minimap_for_fastqs
from reads_pipeline.paths import (
    MINIMAP2_BIN,
    get_vcfs_per_sample_dir,
    get_joint_gatk_segments_bed,
)


def _prepare_snv_calling_test_dir(project_dir):
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
    return {
        "genome_fasta": genome_fasta,
    }


def test_snv_calling_per_sample():
    with tempfile.TemporaryDirectory(prefix="gatk_test") as project_dir:
        res = _prepare_snv_calling_test_dir(project_dir)
        genome_fasta = res["genome_fasta"]
        project_dir_path = Path(project_dir)
        res = do_snv_calling_per_sample(
            project_dir=project_dir_path,
            genome_fasta=genome_fasta,
            num_snvs_in_parallel=2,
        )
        assert res["samples_done"]["sample1"]["out_vcf"].exists()
        assert len(res["samples_done"]) == 1

        res = do_snv_calling_per_sample(
            project_dir=project_dir_path, genome_fasta=genome_fasta
        )
        assert len(res["samples_done"]) == 0

        res = do_snv_calling_per_sample(
            project_dir=project_dir_path,
            genome_fasta=genome_fasta,
            re_run=True,
        )
        assert len(res["samples_done"]) == 1


def test_per_sample_snv_calling_script():
    with tempfile.TemporaryDirectory(prefix="snp_pipeline_test") as project_dir:
        _prepare_snv_calling_test_dir(project_dir)
        cmd = [
            "uv",
            "run",
            "do_per_sample_snv_calling",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)


def test_add_sample_snv_calls_to_db():
    with tempfile.TemporaryDirectory(prefix="gatk_db_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT6_DIR, project_dir_path, dirs_exist_ok=True)
        vcfs = [
            path
            for path in get_vcfs_per_sample_dir(project_dir_path).iterdir()
            if str(path).endswith(".vcf.gz")
        ]
        create_db_with_independent_sample_snv_calls(
            vcfs,
            project_dir=project_dir_path,
            genome_fai_path=PROJECT6_GENOME_FAI,
            mode=GATKDBFileMode.CREATE,
        )
        samples_in_db = get_samples_in_gatk_db(project_dir)
        assert samples_in_db == ["sample1"]

        joint_vcf = project_dir_path / "joint.vcf.gz"
        do_svn_joint_genotyping_for_all_samples_together(
            project_dir_path, genome_fasta=PROJECT6_GENOME_FASTA, out_vcf=joint_vcf
        )
        filtered_vcf = project_dir_path / "joint.filtered.vcf.gz"
        filters = {
            "QD2": "QD < 2.0",  # Quality by Depth
            "FS60": "FS > 60.0",  # Fisher Strand Bias
            "MQ40": "MQ < 40.0",  # Mapping Quality
            "MQRankSum": "MQRankSum < -12.5",  # Mapping Quality Rank Sum
            "ReadPosRankSum": "ReadPosRankSum < -20.0",  # Mapping Quality Rank Sum
        }

        filter_vcf_with_gatk(
            joint_vcf,
            filtered_vcf,
            genome_fasta=PROJECT6_GENOME_FASTA,
            filters=filters,
            project_dir=project_dir,
        )

        bed = get_joint_gatk_segments_bed(project_dir)
        with bed.open("wt") as fhand:
            fhand.write("SL4.0ch01\t1\t1000\n")
            fhand.write("SL4.0ch01\t1001\t2000\n")
            fhand.flush()

        cmd = [
            "uv",
            "run",
            "create_gatk_joint_vcfs_per_segment",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)


def test_create_gatk_db_script():
    with tempfile.TemporaryDirectory(prefix="gatk_db_test") as project_dir:
        project_dir_path = Path(project_dir)
        shutil.copytree(TEST_PROJECT6_DIR, project_dir_path, dirs_exist_ok=True)
        cmd = [
            "uv",
            "run",
            "create_gatk_db",
            project_dir,
        ]
        run(cmd, cwd=project_dir, check=True)
