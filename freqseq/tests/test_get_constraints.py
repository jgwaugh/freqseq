import pytest

from freqseq import get_test_constraints
from freqseq.hypothesis import get_p_sucess

"""
All tests are based off of Evan's online calculator, 
found here: https://www.evanmiller.org/ab-testing/sequential.html

Note that Evan's online app is only for positive lift with unbiased
treatment probability - need to use simulation on top of this to 
handle the other cases. 
"""


@pytest.mark.parametrize(
    "delta, alpha, power, expected_N, expected_d",
    [
        (0.4, 0.05, 0.8, 243, 31),
        (0.5, 0.05, 0.8, 170, 26),
        (0.15, 0.06, 0.85, 1491, 73),
        (0.9, 0.03, 0.6, 52, 16),
    ],
)
def test_build_test(delta, alpha, power, expected_N, expected_d):
    p_success = get_p_sucess(0.5, delta)
    N, d = get_test_constraints(alpha, power, 0.5, p_success)

    assert (N == expected_N) & (d == expected_d)
