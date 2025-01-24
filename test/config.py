from pathlib import Path

TEST_DATA_DIR = Path(__file__).absolute().parent / "data"
TEST_PROJECT1_DIR = TEST_DATA_DIR / "project1"
TEST_PROJECT2_DIR = TEST_DATA_DIR / "project2"
MINIMAP_PROJECT2_TOMATO_INDEX = (
    TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta.mmi"
)
MINIMAP_PROJECT2_TOMATO_FASTA = TEST_PROJECT2_DIR / "tomato_SL12_400000_700000.fasta"
