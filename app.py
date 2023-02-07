import matplotlib.pyplot as plt
import numpy as np

from freqseq import get_test_constraints
from freqseq.hypothesis import (
    compute_transformed_walk_parameters,
    get_barrier_crossing_rate,
    get_p_sucess,
    simulate_walk,
)

# parameters
p = 0.1
delta = 0.8
# assert delta < 1/p - 2
J = 5000
p_success = get_p_sucess(p, delta)

p_star, sigma = compute_transformed_walk_parameters(p, delta)

print(f"modified probabiliy value is {p_star}")

N, d = get_test_constraints(0.15, 0.88, 0.5, p_star)
print(f"the original (N, d) is {(N, d)}")
print(f"modified probabiliy value is {p_star}")


print(f"the new N is {N}")
print(f"the new d is {d}")

d_0 = [d] * N


X = np.array(list(range(N)))
expectation = X * (2 * p - 1)


fpr = get_barrier_crossing_rate(p, N, d, J, mu=expectation, sigma=sigma)
tpr = get_barrier_crossing_rate(p_success, N, d, J, mu=expectation, sigma=sigma)


print(f"True Positives: {tpr}")
print(f"False Positives: {fpr}")


walks = simulate_walk(p_success, N, 5)


f = plt.figure(dpi=100)
plt.title("Rejection Region With Bias")
plt.plot(X, d_0, label="d")
plt.plot(X, X * (2 * p_star - 1), label="E[X]")

for j in range(len(walks)):
    plt.plot(X, (walks[j, :] - expectation) / sigma)

plt.legend()
plt.show()
# f.savefig("images/walk_with_bias.png")
