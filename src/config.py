from dataclasses import dataclass
from pathlib import Path

@dataclass
class AreaConfig:
    boolean_file: Path
    rank_file: Path
    common_column: str
    outdir: Path
    processes: int = 4
    permutations: int = 1000
    seed: int = 42
    gpu: bool = False
    verbose: bool = False
    include_rank_columns_file: Path | None = None
    include_boolean_columns_file: Path | None = None
    include_sample_file: Path | None = None

