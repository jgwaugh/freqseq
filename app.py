import numpy as np
import streamlit as st
from skopt import gp_minimize

from freqseq.calibration import build_test, objective_function


@st.cache
def get_test_results(p, delta, alpha, beta):
    """Caches test constrains for use in app"""
    N, d, sigma, fpr, tpr = build_test(p, delta, alpha, beta, 5000)
    initial_error = np.abs(tpr - beta) + np.abs(fpr - alpha)

    if initial_error <= 0.04:
        return N, d, sigma, fpr, tpr

    def error_function(x):
        return objective_function(p, delta, alpha, beta, x[0], x[1])

    alpha_min = alpha / 20
    alpha_max = min(1, 6 * alpha)
    beta_min = beta - beta / 8
    beta_max = min(1, beta + beta / 8)

    res = gp_minimize(
        error_function,
        [(alpha_min, alpha_max), (beta_min, beta_max)],
        x0=[alpha, beta],
        n_calls=20,
        random_state=777,
        verbose=True,
    )

    best_alpha = res.x[0]
    best_beta = res.x[1]

    return build_test(p, delta, best_alpha, best_beta)


st.write(
    """
    # Frequentist Sequential Testing
    
    This app solves for $N$, $d$, and $\sigma$ required to perform
    a frequentist sequential test. 
    
    To start, specify the parameters below. 
    
    Note that the app assumes
    
    $$ p_t = (1 + \delta)p_c $$
    
    Where $p_t$ and $p_c$ are probability of treatment and control events
    occuring, respectively. Call $T$ and $C$ the total number of conversions in 
    the treatment and control groups at any given time. 
    
    Reject $H_0 : p_c = p_t$ for $\delta > 0$ when $S > d$
    
    $$ S = (T - C - (T + C)(2p -1))\sigma^{-1} $$
    
    For $\delta < 0$, reject when $S' > d$
    
    $$ S' = (C - T - (T + C)(2(1 - p) -1)\sigma^{-1} $$
    
    If $T + C > N$ and $d$ has yet to be crossed, the test fails. 
    """
)

p = st.number_input(
    "Insert the probability of treatment, $p$",
    min_value=0.01,
    max_value=0.99,
    value=0.5,
)

delta = st.slider(
    "Insert the minimum detectable effect",
    min_value=-2.0,
    max_value=2.0,
    step=0.1,
    value=0.3,
)
if delta == 0:
    raise ValueError(f"Delta cannot be equal to zero!")

alpha = st.number_input(
    "Insert the maximum false positive rate", min_value=0.0, max_value=1.0, value=0.05
)

beta = st.number_input(
    "Insert the minimum true positive rate", min_value=0.0, max_value=1.0, value=0.8
)


N, d, sigma, fpr, tpr = get_test_results(p, delta, alpha, beta)


st.write(
    """
# Test Specifications
"""
)
col1, col2, col3 = st.columns(3)
col1.metric("N", int(N))
col2.metric("d", int(d))
col3.metric("sigma", np.round(sigma, 3))
col1.metric("Empirical false positive rate", np.round(fpr, 2))
col2.metric("Empirical True Positive Rate", np.round(tpr, 2))
