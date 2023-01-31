import matplotlib.pyplot as plt

import numpy as np


np.random.seed(4583)

# parameters
p = 0.4
delta = -0.1
N = 500
p_success = p * (1- delta) / (1 - p * delta)


# walk
steps = 2 * (np.random.uniform(0, 1, N) <= p_success).astype(int) - 1
walk_path = np.cumsum(steps)
X = np.array(list(range(N)))
expectation = X * (2*p - 1)

# plot
f = plt.figure(dpi=100)
plt.plot(X, walk_path, label='Walk')
plt.plot(X, expectation, label="Expectation: y = x(2p-1)")
plt.title(f"Random Walk with Treatment Probability = {p} \n and delta = {delta}")
plt.legend()
plt.xlabel('Total Conversions')
plt.ylabel('Position')
f.savefig(f'images/p_{p}_delta_{delta}_random_walk.png')

