import sys
from pathlib import Path
import tomllib
import copy

import reads_pipeline

DEFAULTS = {
    "general": {"re_run": False},
    "fastp": {
        "min_read_len": 30,
        "num_threads": 3,
        "trim_front1": 0,
        "trim_front2": 0,
        "trim_tail1": 0,
        "trim_tail2": 0,
    },
    "minimap": {"num_threads": 3},
    "samtools": {
        "sort_num_threads": 8,
        "duplicates_num_threads": 4,
        "calmd_num_threads": 2,
    },
    "fastqc": {"num_threads": 6},
}


class PipelineConfig:
    def __init__(self, config: dict, default_project_dir: Path):
        self._config = copy.deepcopy(config)
        self._default_project_dir = default_project_dir

    def __getitem__(self, key):
        value = None

        if key == "project_dir":
            value = Path(self._config.get(key, self._default_project_dir))

        minimap_config = self._config["minimap"]
        if "index_path" not in minimap_config:
            raise RuntimeError(
                "index_path is a required argument in the minimap section of the config file"
            )
        else:
            minimap_config["index_path"] = Path(minimap_config["index_path"])

        samtools_config = self._config["samtools"]
        if "deduplicate" not in samtools_config:
            raise RuntimeError(
                "deduplicate is a required bool argument in the samtools section of the config file"
            )

        for tool in DEFAULTS.keys():
            tool_config = self._config[tool]
            for one_key in DEFAULTS[tool].keys():
                tool_config[one_key] = tool_config.get(one_key, DEFAULTS[tool][one_key])

        if "general" not in self._config:
            raise RuntimeError["A general section is required in the config file"]
        if "genome_path" not in self._config["general"]:
            raise RuntimeError(
                "genome_path is a required argument in the general section of the config file"
            )
        else:
            self._config["general"]["genome_path"] = Path(
                self._config["general"]["genome_path"]
            )

        if value is None:
            value = self._config[key]

        return value


def run_pipeline(config):
    reads_pipeline.run_fastqc(
        project_dir=config["project_dir"], re_run=config["general"]["re_run"], threads=1
    )
    reads_pipeline.collect_fastqc_stats(project_dir=config["project_dir"])
    reads_pipeline.run_fastp_minimap(
        project_dir=config["project_dir"],
        minimap_index=config["minimap"]["index_path"],
        genome_fasta=config["general"]["genome_path"],
        deduplicate=config["samtools"]["deduplicate"],
        min_read_len=config["fastp"]["min_read_len"],
        fastp_num_threads=config["fastp"]["num_threads"],
        fastp_trim_front1=config["fastp"]["trim_front1"],
        fastp_trim_tail1=config["fastp"]["trim_tail1"],
        fastp_trim_front2=config["fastp"]["trim_front2"],
        fastp_trim_tail2=config["fastp"]["trim_tail2"],
        minimap_num_threads=config["minimap"]["num_threads"],
        sort_num_threads=config["samtools"]["sort_num_threads"],
        duplicates_num_threads=config["samtools"]["duplicates_num_threads"],
        calmd_num_threads=config["samtools"]["calmd_num_threads"],
        re_run=config["general"]["re_run"],
    )
    reads_pipeline.collect_fastp_stats(project_dir=config["project_dir"])
    reads_pipeline.collect_cram_stats(project_dir=config["project_dir"])


if __name__ == "__main__":
    args = sys.argv
    if len(args) == 1:
        config_path = "pipeline.toml"
    elif len(args) == 2:
        config_path = args[1]
    else:
        msg = f"Usage: {args[0]} [pipeline.toml]"
        print(msg)
        sys.exit(1)
    config_path = Path(config_path)

    if not config_path.exists():
        msg = f"The config file {config_path} does not exist"
        print(msg)
        sys.exit(2)

    default_project_dir = config_path.parent

    with config_path.open("rb") as fhand:
        config = PipelineConfig(
            tomllib.load(fhand), default_project_dir=default_project_dir
        )

    run_pipeline(config)
