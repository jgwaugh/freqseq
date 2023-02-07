import numpy as np

from freqseq.hypothesis import (
get_p_sucess,
get_barrier_crossing_rate,
compute_transformed_walk_parameters
)

from freqseq.search import get_test_constraints


def objective_function(
        p: float,
        delta: float,
        alpha: float,
        beta: float,
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
    J: int
        Number of simulation iterations

    Returns
    -------
    float
        The error between alpha and beta and the empirical values used with the test

    """

    p_success = get_p_sucess(p, delta)

    p_star, sigma = compute_transformed_walk_parameters(p, delta)


    N, d = get_test_constraints(alpha, beta, 0.5, p_star)

    if np.isnan(N) or np.isnan(d):
        return 5


    X = np.arange(N)
    expectation = X * (2 * p - 1)


    fpr = get_barrier_crossing_rate(p, N, d, J, mu=expectation, sigma=sigma)
    tpr = get_barrier_crossing_rate(p_success, N, d, J, mu=expectation, sigma=sigma)

    return (tpr - beta)**2 + (fpr - alpha)**2
