from typing import Tuple

import numpy as np
from scipy.special import betaln

MAX_BARRIER = 5000
MAX_CONVERSIONS = 800000


def search_for_barrier(
    z_low: int,
    z_high: int,
    alpha: float,
    power_level: float,
    null_p: float,
    alt_p: float,
) -> int:
    """
    Uses a binary search to choose a value of z (aka d) that satisfies significance
    and power constraints.

    Parameters
    ----------
    z_low: int
        Lower bound for the search
    z_high: int
        Upper bound for the search
    alpha: float
        Significance constraint
    power_level: float
        Power constraint
    null_p: float
        Probability of positive conversion under H_0
    alt_p : float
        Probability of positive conversion under H_1

    Returns
    -------
    int
        cutoff value

    """
    log_null_p = np.log(null_p)
    log_null_1_p = np.log(1 - null_p)

    log_alt_p = np.log(alt_p)
    log_alt_1_p = np.log(1 - alt_p)

    z = z_low + 2 * np.floor((z_high - z_low) / 4)

    while z_low < z_high:
        null_cdf = 0
        alt_cdf = 0
        old_low = z_low
        old_high = z_high

        for n in range(int(z), MAX_CONVERSIONS + 1, 1000):

            k = 0.5 * (n + z)
            prefix = z / n / k
            lbeta_k = betaln(k, n + 1 - k)

            null_cdf += prefix * np.exp(
                -lbeta_k + (k - z) * log_null_1_p + k * log_null_p
            )
            alt_cdf += prefix * np.exp(-lbeta_k + (k - z) * log_alt_1_p + k * log_alt_p)

            if np.isnan(null_cdf) | np.isnan(alt_cdf):
                break

            if alt_cdf > power_level:
                if null_cdf < alpha:
                    z_high = z
                else:
                    z_low = z + 2

                break
            elif null_cdf > alpha:
                z_low = z + 2

        if (np.isnan(null_cdf)) | (np.isnan(alt_cdf)) | (n >= MAX_CONVERSIONS):
            print("NaN...")
            break
        print(f"High: {z_high}, Low: {z_low}, Z: {z}, null: {null_cdf}, alt: {alt_cdf}")
        # import ipdb; ipdb.set_trace()
        z = z_low + 2 * np.floor((z_high - z_low) / 4)

    return z


def get_conversions_for_specified_barrier(
    z: int, alpha: float, power_level: float, null_p: float, alt_p: float
) -> int:
    """
    Iterates through conversion bounds ("N" from the notes) to find minimum
    conversion bound for a specified barrier (z) that meets given power and significance
    constraints

    Parameters
    ----------
    z: int
        Cutoff barrier for test
    alpha: float
        Significance constraint
    power_level: float
        Power constraint
    null_p: float
        Positive conversion rate under H_0
    alt_p: float
        Positive conversion rate under H_1

    Returns
    -------
    int
        The max number of conversions for the test satisfying power/significance/barrier constraints

    """

    null_cdf = 0
    alt_cdf = 0

    log_null_p = np.log(null_p)
    log_null_1_p = np.log(1 - null_p)

    log_alt_p = np.log(alt_p)
    log_alt_1_p = np.log(1 - alt_p)

    for n in range(int(z), MAX_CONVERSIONS, 2):
        k = 0.5 * (n + z)

        prefix = z / n / k
        lbeta_k = betaln(k, n + 1 - k)

        null_cdf += prefix * np.exp(-lbeta_k + (k - z) * log_null_1_p + k * log_null_p)
        alt_cdf += prefix * np.exp(-lbeta_k + (k - z) * log_alt_1_p + k * log_alt_p)

        if np.isnan(null_cdf) | np.isnan(alt_cdf):
            return np.nan
        if alt_cdf > power_level:
            if null_cdf < alpha:
                return n

        return np.nan


def get_test_constraints(
    alpha: float, power_level: float, null_p: float, alternative_p: float
) -> Tuple[int, int]:
    """
    Computes test boundaries for given constraints

    Parameters
    ----------
    alpha: float
        Significance constraint
    power_level: float
        Power constraint
    null_p: float
        Probability of conversion under H_0
    alternative_p: float
        Probability of conversion under H_1

    Returns
    -------
    tuple
        tuple of (conversion constraint, barrier constraint)

    """

    best_odd_z = search_for_barrier(
        1, MAX_BARRIER, alpha, power_level, null_p, alternative_p
    )
    best_even_z = search_for_barrier(
        2, MAX_BARRIER + 1, alpha, power_level, null_p, alternative_p
    )

    odd_n = get_conversions_for_specified_barrier(
        best_odd_z, alpha, power_level, null_p, alternative_p
    )
    even_n = get_conversions_for_specified_barrier(
        best_even_z, alpha, power_level, null_p, alternative_p
    )

    if np.isnan(odd_n) | (even_n < odd_n):
        return even_n, best_even_z
    else:
        return odd_n, best_odd_z
