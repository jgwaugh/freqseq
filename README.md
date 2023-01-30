# Frequentist Sequential Testing

## Intro

These notes cover frequentist sequential hypothesis testing with bias treatment allocation. Credit for the initial derivation goes to [Evan Miller](https://www.evanmiller.org/sequential-ab-testing.html#notes). I've simply gone over his derivation in more granularity, added calculations dealing with treatment bias, and written a `python` implementation. 


## Problem

Suppose we are running a random experiment to determine the efficacy of some intervention. We would like to end this experiment early if results look promising without "peaking" (peaking incurs bias by effectively testing multiple hypotheses). These notes give a derivation of a testing approach. 


### Random Walk

We construct a sequential test using a random walk. More specifically, suppose we have two groups, $C$ and $T$, or control and treatment. We assign individuals to $C$ with probability $1- p$ and individuals to $T$ with probability $p$. 

Individuals in each group "convert" with rates $p_c$ and $p_t$. 

We are interested in 

$$ H_0 : p_c = p_t $$

and 

$$ H_1: p_c > p_t $$


Each individual is then a bernoulli random variable, $X_i$, with $P(X_i = 1) = pp_t$ and $P(X_i = -1) = (1-p)p_c$. 



$$S_k = \sum_{i=1}^k{X_i}$$

 is then a biased random walk. Under $H_0$, $p_c = p_t$, so we can effectively estimate $S_k$ with bias $s$.

Note that the following is an unbiased random walk:

$$ \tilde{S_k} = \sum_{i=1}^k{X_i} - k(2p -1) $$

as $X_i$ follows a  $Bernoulli(p)$ distribution. 


### Test Statistic

The test is defined by choosing a bound, $d$, and a number of conversions, $N$, such that the probability of the walk escaping the region under $H_0$ is less thatn $\alpha$ for some predefined false positive rate. 

More specifically, define $r_{n, d}$ as

$$ r_{n, d} = \frac{d}{n} {n \choose \frac{n + d}{2}} p ^ {\frac{n + d}{2}} (1 - p)^{\frac{n - d}{2}} $$

$r_{n, d}$ is the probability of reaching $k$ for the very first time after $n$ iterations of the random walk. The basic idea is that this requires $d$ treatment conversions and then a balance of $\frac{n - d}{2}$ treatment converisons and $\frac{n - d}{2}$ control conversions (so a total of $\frac{n + d}{2}$ treatment conversions). The combinatorial handles the order and the term $\frac{d}{n}$ controls for the fact that only $\frac{d}{n}$ of the ${n \choose \frac{n + d}{2}}$ paths arrive at $k$ conversions at exactly time $n$. For more information, see Chapter 3 of [this book](https://bitcoinwords.github.io/assets/papers/an-introduction-to-probability-theory-and-its-applications.pdf).

Next, define $R_{N, d}$ as 

$$ R_{N, d} = \sum_{n = 1} ^Nr_{n, d} $$

This is the probabilty of escaping the boundary $d$ in less than $N$ iterations. 

We can then choose $N$ and $d$ such that for some $\alpha$, 

$$ R_{N, d} < \alpha $$

Then if $S_k$ crosses $d$ for any $k \leq N$, we reject $H_0$. 

### Power

There are an infinite number of pairs $(N, d)$ that satisfy the significance equation. 

We can the pair to use by adding the following constraint:

$$ P(S_k > d, k \leq N | H_1)  > \beta$$

Where $\beta$ is the probability of rejecting the null under the alternative hypothesis. 


Under $H_1$, we need to solve for $p_c$ and $p_t$. For example, we want customers who take our drug to die at a 10% lower rate than customers in the control group. First, specify some minimum amount to detect, $\delta$. We can then write 

$$ p_t = (1 - \delta)p_c$$

$S_k$ steps up when a conversion takes place in the treatment group. $S_k$ only steps in *either* direction when a conversion takes place. If at time $t$, a customer is assigned to the control group and does **not** convert, the walk does not move. Thus, group assignement **and** conversion rate dictate how the walk moves. 


Define $Y \sim Bernoulli(p)$ indicating assignment to treatment or control. Define $Z_t \sim Bernoulli(p * p_t)$  and $Z_c \sim Bernoulli((1-p)p_c)$. $Z_t$ and $Z_c$ correspond to conversion in the treatment and control groups, respectively. 


To be more precise, the walk steps up under the event $Z_t  = 1 | Z_t + Z_c = 1$. The probability of this event occuring, $p^*$, is defined as 

$$ p^* = P(Z_t = 1 | Z_t + Z_c = 1) = 
\frac{P(Z_t + Z_c = 1 | Z_t = 1)P(Z_t = 1)}
{P(Z_t + Z_c = 1)} = 
\frac{p * p_t}{ p * p_t + (1 - p) * p_c}



