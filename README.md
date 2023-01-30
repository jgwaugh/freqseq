# Frequentist Sequential Testing

## Intro

These notes cover frequentist sequential hypothesis testing with bias treatment allocation. Credit for the initial derivation goes to [Evan Miller](https://www.evanmiller.org/sequential-ab-testing.html#notes). I've simply gone over his derivation in more granularity, added calculations dealing with treatment bias, and written a `python` implementation. 


## Problem

Suppose we are running a random experiment to determine the efficacy of some intervention. We would like to end this experiment early if results look promising without "peaking" (peaking incurs bias by effectively testing multiple hypotheses). These notes give a derivation of a testing approach. 


### Random Walk

We construct a sequential test using a random walk. More specifically, suppose we have to groups, $C$ and $T$, or control and treatment. We assign individuals to $C$ with probability $1- p$ and individuals to $T$ with probability $p$. 

Individuals in each group "convert" (could mean they die, or they buy your product - whatever the context) with rates $p_c$ and $p_t$. 

We are interested in 

$$ H_0 : p_c = p_t $$

and 

$$ H_1: p_c < p_t $$


Each individual is then a bernoulli random variable, $X_i$, with $P(X_i = 1) = pp_t$ and $P(X_i = -1) = (1-p)p_c$. 



$$S_k = \sum_{i=1}^k{X_i}$$

 is then a biased random walk. Under $H_0$, $p_c = p_t$, so we can effectively estimate $S_k$ with bias $s$.

Note that the following is an unbiased random walk:

$$ \tilde{S_k} = \sum_{i=1}^k{X_i} - k(2p -1) $$

As $X_i$ follows a  $Bernoulli(p)$ distribution. 


#### Test Statistic

The test is defined by choosing a bound, $d$, and a number of conversions, $N$, such that the probability of the walk escaping the region under $H_0$ is less thatn $\alpha$ for some predefined false positive rate. 

More specifically, define $r_{n, k}$ as

$$ r_{n, k} = \frac{k}{n} {n \choose \frac{n + k}{2}} p ^ {\frac{n + k}{2}} (1 - p)^{\frac{n - k}{2}} $$

$r_{n, k}$ is the probability of reaching $k$ for the very first time after $n$ iterations of the random walk. The basic idea is that this requires $k$ treatment conversions and then a balance of $\frac{n - k}{2}$ treatment converisons and $\frac{n - k}{2}$ control conversions (so a total of $\frac{n + k}{2}$ treatment conversions). The combinatorial handles the order and the term $\frac{k}{n}$ controls for the fact that only $\frac{k}{n}$ of the ${n \choose \frac{n + k}{2}}$ paths arrive at $k$ conversions at exactly time $n$. For more information, see Chapter 3 of [this book](https://bitcoinwords.github.io/assets/papers/an-introduction-to-probability-theory-and-its-applications.pdf).





