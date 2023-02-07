from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

RowVector = NDArray


def simulate_walk(p_success: float, N: int, J: int) -> NDArray:
    """
    Simulates a random walk with probability of stepping up `p_success`

    Parameters
    ----------
    p_success: float
        The probability of the walk stepping up
    N: int
        The length of the walks to simulate
    J: int
        The number of walks to simulate
    Returns
    -------
    NDArray
        Matrix of position for each simulation iteration

    """
    steps = 2 * (np.random.uniform(0, 1, (J, N)) <= p_success).astype(int) - 1
    walk = np.cumsum(steps, axis=1)

    return walk


def get_barrier_crossing_rate(
    p_success: float,
    N: int,
    d: Union[RowVector, float],
    J: int,
    mu: Optional[RowVector] = None,
    sigma: Optional[float] = None,
    crosses_upper: bool = True,
) -> float:
    """
    Simulates a random walk and returns the number of times it crosses the barrier
    which can be used to compute the false positive and true positive rates.

    mu and sigma refer to mean and variance parameters, which can be used
    to reshape the walk.


    Parameters
    ----------
    p_success: float
        The probability of the walk stepping up
    N: int
        The conversion test boundary, aka the max number of conversions
    d: NDArray
        Row vector containing the boundary at a given number of conversions
    J: int
        The number of walks to simulate
    mu : NDArray
        Mean of the random walk as a function of conversions
    sigma: float
        Variance of the random walk as a function of conversions
    crosses_upper: bool
        Boolean indicating if the test ends when the walk goes above
        or below the boundary


    Returns
    -------
    float
        The rate at which the walk crosses boundary, over all simulations

    """

    walk = simulate_walk(p_success, N, J)

    if not isinstance(mu, type(None)):
        mu = np.zeros(N)
    if not isinstance(sigma, type(None)):
        sigma = 1

    walk = (walk - mu) / sigma

    if crosses_upper:
        crosses = (walk >= d).astype(int)
    else:
        crosses = (walk <= d).astype(int)
    return np.mean(crosses.argmax(axis=1) != 0)


def get_p_sucess(p: float, delta: float) -> float:
    """
    Gets the probaiblity of stepping up given a relative
    difference between two rates

    Parameters
    ----------
    p: float
        Probability of being assigned to the treatment group
    delta: float
        The relative difference between treatment and control
        p_t = (1 + delta)p_c

    Returns
    -------
    float
        p(treatment = 1 | treatment + control = 1)

    """

    return p * (1 + delta) / (1 + p * delta)
