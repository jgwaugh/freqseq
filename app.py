import numpy as np

from freqseq import get_test_constraints
from freqseq.hypothesis import get_barrier_crossing_rate, get_p_sucess, simulate_walk

# parameters
p = 0.1
delta = 0.8
# assert delta < 1/p - 2
J = 5000
p_success = get_p_sucess(p, delta)

w = 2 * p - 1
v = 2 * p_success - 1

sigma = 1 - v**2 + (v - w) ** 2
sigma = sigma ** (1 / 2)
u = (v - w) / sigma
p_star = (u + 1) / 2


print(f"modified probabiliy value is {p_star}")

N, d = get_test_constraints(0.15, 0.88, 0.5, p_star)
print(f"the original (N, d) is {(N, d)}")
print(f"modified probabiliy value is {p_star}")
# d -= 10
# N += 100

# import ipdb; ipdb.set_trace()
# d = d * DELTA/ cos_theta
# N = int(N * s)
# d = 38
# d -= 10


print(f"the new N is {N}")
print(f"the new d is {d}")

d_0 = [d] * N
# d = np.arange(0, N) * (2 * p - 1) + d


X = np.array(list(range(N)))
expectation = X * w


fpr1 = get_barrier_crossing_rate(p, N, d, J, mu=expectation, sigma=sigma)
tpr1 = get_barrier_crossing_rate(p_success, N, d, J, mu=expectation, sigma=sigma)


print(f"True Positives: {tpr}")
print(f"False Positives: {fpr}")
# print(set(d.tolist()[0]))
# print(N)

walks = simulate_walk(p_success, N, 5)
X = np.arange(0, N)

import matplotlib.pyplot as plt

f = plt.figure(dpi=100)
plt.title("Rejection Region With Bias")
plt.plot(X, d_0, label="d")
plt.plot(X, X * (2 * p_star - 1), label="E[X]")

for j in range(len(walks)):
    plt.plot(X, (walks[j, :] - expectation) / sigma)

plt.legend()
plt.show()
# f.savefig("images/walk_with_bias.png")
