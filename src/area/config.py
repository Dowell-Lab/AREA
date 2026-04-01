"""
Configuration file loader for AREA.

Reads a YAML config file and returns a flat dictionary whose keys match
the argparse dest names (underscored, e.g. ``boolean_file``).
"""

import yaml


# Maps YAML key → argparse dest name.  Keys that already match are not
# listed here (e.g. "threads", "gpu", "verbose").
_KEY_TO_DEST = {
    "boolean_file":      "boolean_file",
    "rank_file":         "rank_file",
    "join_column":       "join_column",
    "out_dir":           "out_dir",
    "threads":           "threads",
    "gpu":               "gpu",
    "keep_rank_columns": "keep_rank_columns",
    "keep_bool_columns": "keep_bool_columns",
    "keep_samples":      "keep_samples",
    "verbose":           "verbose",
}


def load_config(path: str) -> dict:
    """Read a YAML config file and return a validated dict.

    Parameters
    ----------
    path : str
        Path to the YAML file.

    Returns
    -------
    dict
        Keys are argparse dest names, values are the parsed settings.
        ``None``/``null`` values are preserved so they can be
        distinguished from "not set".

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file contains an unrecognised key.
    """
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}

    config = {}
    for key, value in raw.items():
        if key not in _KEY_TO_DEST:
            raise ValueError(
                f"Unrecognised config key '{key}'. "
                f"Valid keys: {sorted(_KEY_TO_DEST.keys())}"
            )
        config[_KEY_TO_DEST[key]] = value

    return config
