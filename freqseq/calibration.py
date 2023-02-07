from typing import Tuple

import numpy as np

from freqseq.hypothesis import (
    compute_transformed_walk_parameters,
    get_barrier_crossing_rate,
    get_p_sucess,
)
from freqseq.search import get_test_constraints


def build_test(
    p: float,
    delta: float,
    alpha: float,
    beta: float,
    J: int = 5000,
) -> Tuple[float, int]:
    """

    Builds the random walk statistical test and verifies its false positive
    and true positive rates using simulation

    Parameters
    ----------
    p: float
        Treatment probability
    delta: float
        Effect size
    alpha: float
        Desired false positive rate
    beta: float
        Desired true positive rate
    J: int
        Number of simulation iterations

    Returns
    -------
    tuple
        Tuple containing N, d, the variance conversion factor, the empirical false positive rate,
        and the empirical true positive rate

    """

    p_success = get_p_sucess(p, delta)

    p_star, sigma = compute_transformed_walk_parameters(p, delta)

    N, d = get_test_constraints(alpha, beta, 0.5, p_star)

    if np.isnan(N) or np.isnan(d):
        return np.nan, np.nan, np.nan, 5, 5


    if delta < 0:
        p = 1 - p

    X = np.arange(N)
    expectation = X * (2 * p - 1)

    fpr = get_barrier_crossing_rate(p, N, d, J, mu=expectation, sigma=sigma)
    tpr = get_barrier_crossing_rate(p_success, N, d, J, mu=expectation, sigma=sigma)

    print(fpr, tpr)


    return N, d, sigma, fpr, tpr


def objective_function(
    p: float,
    delta: float,
    alpha: float,
    beta: float,
    calibrated_alpha: float,
    calibrated_beta: float,
    J: int = 5000,
) -> float:
    """
    For given treatment probability, minimum detect size, desired true positive rate,
    and desired false positive rate, calculates the error between the actual
    false positive and true positive error rates and empirical error rates, from simulation.

    Parameters
    ----------
    p: float
        Treatment probability
    delta: float
        Effect size
    alpha: float
        Desired false positive rate
    beta: float
        Desired true positive rate
    calibrated_alpha: float
        The calibrated false positive constraint
    calibrated_beta: float
        The calibrated true positive constraint
    J: int
        Number of simulation iterations

    Returns
    -------
    float
        The error between alpha and beta and the empirical values used with the test

    """

    N, d, sigma, fpr, tpr = build_test(p, delta, calibrated_alpha, calibrated_beta, J)
    return np.abs(tpr - beta) + np.abs(fpr - alpha)
