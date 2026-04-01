from __future__ import annotations
import importlib
import numpy as np

#this script just handles the use_gpu flag. If its true we have to use cupy instead of numpy for a bunch of stuff. 
#GPUs are only useful when you have lots of arrtibute vs bool comparisons 


def get_array_backend(use_gpu: bool = False):
    if not use_gpu:
        return np

    cupy_spec = importlib.util.find_spec("cupy")
    if cupy_spec is None:
        raise ImportError(
            "GPU mode requested, but CuPy is not installed. "
            "Install cupy or run with gpu=False."
        )

    import cupy as cp
    return cp


def is_gpu_backend(xp) -> bool:
    return xp.__name__ == "cupy"


def to_numpy(x, xp=None):
    if x is None:
        return None
    if xp is not None and is_gpu_backend(xp):
        return x.get() if hasattr(x, "get") else x
    return np.asarray(x)


def to_scalar(x, xp=None):
    if x is None:
        return None
    if xp is not None and is_gpu_backend(xp):
        x = x.get() if hasattr(x, "get") else x
    if hasattr(x, "item"):
        return x.item()
    return x


def as_xparray(x, xp):
    return xp.array(x)
