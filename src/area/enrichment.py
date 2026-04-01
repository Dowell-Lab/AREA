"""
Core AREA math: enrichment score, permutation null distribution, and
Normalized Enrichment Score (NES) with p-value.
"""

import numpy as np
from scipy import stats as scipy_stats

from .backend import load_array_backend, to_numpy, trapz


# ---------------------------------------------------------------------------
# Enrichment score
# ---------------------------------------------------------------------------

def compute_enrichment_score(binary_vector, xp=np, verbose=False):
    """Compute the GSEA-like area enrichment score.

    Parameters
    ----------
    binary_vector : array-like
        1-D sequence of 0/1 (or truthy/falsy) values indicating
        comorbidity presence along a rank ordering.
    xp : module
        Array backend (``numpy`` or ``cupy``).
    verbose : bool
        Print intermediate diagnostics.

    Returns
    -------
    enrichment_score : float
    normalized_scores : array
    trend : array
    cumulative_scores : array
    """
    binary = xp.array([1 if v > 0 else 0 for v in binary_vector])
    n = len(binary)
    total_hits = float(xp.sum(binary))
    bin_width = 1.0 / n

    normalized_scores = xp.multiply(xp.divide(binary, total_hits), bin_width)
    cumulative_scores = xp.cumsum(normalized_scores)

    trend = xp.append(xp.arange(0, 1, 1.0 / (n - 1)), 1.0)
    trend = xp.multiply(trend, bin_width)

    enrichment_score = (trapz(cumulative_scores, xp) - trapz(trend, xp)) * 2

    if verbose:
        print(f"  bin_width               : {bin_width}")
        print(f"  sum(normalized_scores)  : {xp.sum(normalized_scores)}")
        print(f"  sum(cumulative_scores)  : {xp.sum(cumulative_scores)}")
        print(f"  sum(trend)              : {xp.sum(trend)}")
        print(f"  len(trend)              : {len(trend)}")
        print(f"  enrichment_score        : {enrichment_score}")

    return enrichment_score, normalized_scores, trend, cumulative_scores


# ---------------------------------------------------------------------------
# Permutation null distribution
# ---------------------------------------------------------------------------

def permute_enrichment_scores(binary_vector, n_permutations=1000, seed=42,
                              xp=np, verbose=False):
    """Generate a null distribution of enrichment scores by permutation.

    Parameters
    ----------
    binary_vector : array-like
        Original binary hit vector.
    n_permutations : int
        Number of random permutations.
    seed : int
        Random seed for reproducibility.
    xp : module
        Array backend.
    verbose : bool
        Passed through to ``compute_enrichment_score``.

    Returns
    -------
    null_scores : list[float]
    """
    xp.random.seed(seed=seed)
    binary_array = xp.array(binary_vector)

    null_scores = []
    for _ in range(n_permutations):
        shuffled = xp.random.permutation(binary_array)
        es, *_ = compute_enrichment_score(shuffled, xp=xp, verbose=verbose)
        null_scores.append(es)

    return null_scores


# ---------------------------------------------------------------------------
# Normalized Enrichment Score + p-value
# ---------------------------------------------------------------------------

def compute_nes_pvalue(observed_es, null_scores, use_gpu=False, xp=np):
    """Compute the Normalized Enrichment Score and one-sided p-value.

    Parameters
    ----------
    observed_es : float
        Enrichment score from the real data.
    null_scores : list[float]
        Null distribution from permutations.
    use_gpu : bool
        Whether values need GPU → CPU transfer before the scipy call.
    xp : module
        Array backend.

    Returns
    -------
    nes : float
    pvalue : float
    """
    if observed_es > 0:
        subset = xp.array([x for x in null_scores if x > 0])
        mu = xp.mean(subset)
        sigma = xp.std(subset)
        nes = -(observed_es / mu)
        pvalue = 1 - scipy_stats.norm.cdf(
            to_numpy(observed_es, use_gpu),
            to_numpy(mu, use_gpu),
            to_numpy(sigma, use_gpu),
        )
    else:
        subset = xp.array([x for x in null_scores if x < 0])
        mu = xp.mean(subset)
        sigma = xp.std(subset)
        nes = observed_es / mu
        pvalue = scipy_stats.norm.cdf(
            to_numpy(observed_es, use_gpu),
            to_numpy(mu, use_gpu),
            to_numpy(sigma, use_gpu),
        )

    return nes, pvalue
