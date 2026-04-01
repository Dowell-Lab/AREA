"""
Array backend selection (CPU via NumPy, GPU via CuPy).

Every module that needs array operations imports these helpers rather
than choosing numpy/cupy on its own.
"""

import numpy as np


def load_array_backend(use_gpu: bool):
    """Return the appropriate array library.

    Parameters
    ----------
    use_gpu : bool
        If True, import and return ``cupy``.  Otherwise return ``numpy``.

    Returns
    -------
    module
        ``numpy`` or ``cupy``.
    """
    if use_gpu:
        import cupy
        return cupy
    return np


def to_numpy(value, use_gpu: bool):
    """Move a scalar from GPU memory to a plain Python/NumPy value.

    Parameters
    ----------
    value
        A numeric scalar (possibly a CuPy 0-d array).
    use_gpu : bool
        If True, call ``value.get()`` to transfer to CPU.

    Returns
    -------
    float or int
    """
    return value.get() if use_gpu else value


def trapz(y, xp):
    """Version-safe trapezoidal integration.

    NumPy >= 2.0 renamed ``trapz`` to ``trapezoid``; CuPy may still
    use the old name.
    """
    if hasattr(xp, "trapezoid"):
        return xp.trapezoid(y)
    return xp.trapz(y)
