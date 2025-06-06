from pathlib import Path

TEST_DATA_DIR = Path(__file__).absolute().parent / "data"
TEST_PROJECT1_DIR = TEST_DATA_DIR / "project1"
TEST_PROJECT2_DIR = TEST_DATA_DIR / "project2"
MINIMAP_PROJECT2_TOMATO_INDEX = (
    TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta.mmi"
)
MINIMAP_PROJECT2_TOMATO_FASTA = TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta"
TEST_PROJECT3_DIR = TEST_DATA_DIR / "project3"
MINIMAP_PROJECT3_GENOME_FASTA = TEST_PROJECT3_DIR / "genomexyztest.fasta"
MINIMAP_PROJECT3_GENOME_INDEX = TEST_PROJECT3_DIR / "genomexyztest.fasta.mmi"
TEST_PROJECT4_DIR = TEST_DATA_DIR / "project4"
MINIMAP_PROJECT4_GENOME_FASTA = TEST_PROJECT3_DIR / "genomexyztest.fasta"
MINIMAP_PROJECT4_GENOME_INDEX = TEST_PROJECT3_DIR / "genomexyztest.fasta.mmi"
TEST_PROJECT5_DIR = TEST_DATA_DIR / "project5"
TEST_PROJECT6_DIR = TEST_DATA_DIR / "project6"
PROJECT6_GENOME_FASTA = TEST_PROJECT6_DIR / "SL4.0ch01_13486265-13488265.fasta"
PROJECT6_GENOME_FAI = TEST_PROJECT6_DIR / "SL4.0ch01_13486265-13488265.fasta.fai"
TEST_PROJECT7_DIR = TEST_DATA_DIR / "project7"
