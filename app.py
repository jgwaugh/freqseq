from skopt import gp_minimize

from freqseq.calibration import build_test, objective_function

p = 0.13
delta = -0.1
alpha = 0.05
beta = 0.8


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


best_results = build_test(p, delta, best_alpha, best_beta)
