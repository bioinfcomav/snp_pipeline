from pathlib import Path
import tomllib

from reads_pipeline.paths import get_config_path

DEFAULTS = {
    "general": {
        "re_run": False,
        "verbose": True,
        "num_mappings_in_parallel": 1,
        "num_snvs_in_parallel": 1,
    },
    "fastp": {
        "min_read_len": 30,
        "num_threads": 3,
        "trim_front1": 0,
        "trim_front2": 0,
        "trim_tail1": 0,
        "trim_tail2": 0,
    },
    "trim_quals": {"num_bases": 3, "qual_reduction": 20},
    "minimap": {"num_threads": 3},
    "samtools": {
        "sort_num_threads": 8,
        "duplicates_num_threads": 4,
        "calmd_num_threads": 2,
        "samtools_stats_num_threads": 4,
    },
    "fastqc": {"num_threads": 6},
    "gatk": {"per_sample_calling_min_mapq": 10},
}


class PipelineConfig:
    def __init__(self, project_dir: Path):
        config_path = get_config_path(project_dir)
        if not config_path.exists():
            raise ValueError(f"The config file {config_path} does not exist")

        with config_path.open("rb") as fhand:
            config = tomllib.load(fhand)

        if "gatk" not in config:
            config["gatk"] = {}

        self._config = config

    def __getitem__(self, key):
        value = None

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
