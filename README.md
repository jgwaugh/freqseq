# Frequentist Sequential Testing

## Intro

These notes cover frequentist sequential hypothesis testing with bias treatment allocation. Credit for the initial derivation goes to [Evan Miller](https://www.evanmiller.org/sequential-ab-testing.html#notes). I've simply gone over his derivation in more granularity, added calculations dealing with treatment bias, and written a `python` implementation. 


## Problem

Suppose we are running a random experiment to determine the efficacy of some intervention. We would like to end this experiment early if results look promising without "peaking" (peaking incurs bias by effectively testing multiple hypotheses). These notes give a derivation of a testing approach. 


### Random Walk

We construct a sequential test using a random walk. More specifically, suppose we have to groups, $C$ and $T$, or control and treatment. We assign individuals to $C$ with probability $1- s$ and individuals to $T$ with probability $s$. 

Individuals in each group "convert" (could mean they die, or they buy your product - whatever the context) with rates $p_c$ and $p_t$. 

We are interested in 

$$ H_0 : p_c = p_t $$

and 

$$ H_1: p_c < p_t $$


Each individual is then a bernoulli random variable, $X_i$, with $P(X_i = 1) = sp_t$ and $P(X_i = -1) = (1-s)p_c$. 



$$S_k = \sum_{i=1}^k{X_i}$$

 is then a biased random walk. 

Under $H_0$, $p_c = p_t$, so we can effectively estimate $S_k$ with bias $s$.



