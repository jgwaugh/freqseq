import pytest

from freqseq import get_test_constraints

"""
All tests are based off of Evan's online calculator, 
found here: https://www.evanmiller.org/ab-testing/sequential.html

Note that Evan's online app is only for positive lift with unbiased
treatment probability - need to use simulation on top of this to 
handle the other cases. 
"""


def get_p_sucess(delta: float) -> float:
    """
    Gets the probaiblity of stepping up given a relative
    difference between two rates

    Parameters
    ----------
    delta: float
        The relative difference between treatment and control
        p_t = (1 - delta)p_c

    Returns
    -------
    float
        p(treatment = 1 | treatment + control = 1)

    """

    return 0.5 * (1 - delta) / (1 - 0.5 * delta)


@pytest.mark.parametrize(
    "delta, alpha, power, expected_N, expected_d",
    [
        (-0.4, 0.05, 0.8, 243, 31),
        (-0.5, 0.05, 0.8, 170, 26),
        (-0.15, 0.06, 0.85, 1491, 73),
        (-0.9, 0.03, 0.6, 52, 16),
    ],
)
def test_build_test(delta, alpha, power, expected_N, expected_d):
    p_success = get_p_sucess(delta)
    N, d = get_test_constraints(alpha, power, 0.5, p_success)

    assert (N == expected_N) & (d == expected_d)
