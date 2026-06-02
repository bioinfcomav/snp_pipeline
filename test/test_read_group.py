import tempfile

from reads_pipeline.read_group import create_minimap_rg_str


def _tags(rg_str):
    return rg_str.replace("\\\\t", "\t").split("\t")


def test_rg_platform_and_library_use_distinct_tags():
    project_dir = tempfile.mkdtemp()
    rg_str = create_minimap_rg_str(
        "read1",
        {"sample": "sample1", "platform": "ILLUMINA", "library": "lib1"},
        project_dir,
    )
    tags = _tags(rg_str)
    # the library must be tagged LB, not a second PL (regression)
    assert "PL:ILLUMINA" in tags
    assert "LB:lib1" in tags
    assert sum(tag.startswith("PL:") for tag in tags) == 1


def test_rg_only_library():
    project_dir = tempfile.mkdtemp()
    rg_str = create_minimap_rg_str(
        "read1", {"sample": "sample1", "library": "lib1"}, project_dir
    )
    tags = _tags(rg_str)
    assert "LB:lib1" in tags
    assert not any(tag.startswith("PL:") for tag in tags)
