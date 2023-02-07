# Notes

The code for this is based on [Evan's online calculator](
    https://www.evanmiller.org/ab-testing/sequential.html
). 

To evaluate the power and significance equations from [before](../notes.md), we need to compute the following term on our machines 
without hitting overflow errors as $n$ gets large:

$$ r_{n, d} = \frac{d}{n} {n \choose \frac{n + d}{2}} p ^ {\frac{n + d}{2}} q^{\frac{n - d}{2}} $$

where $p + q = 1$. 

To do this, define $k  = \frac{n + d}{2}$. 

Then, we have 

$$ \frac{d}{n} \frac{n!}{(n - k!)k!} \space p^{k} q^{n-k}  = 
\frac{d}{nk} \frac{n!}{(n - k)!(k - 1)!} \space p^{k} q^{n-k}
$$

Observe that the beta function, $B(m, n)$ has the following property when $m$ or $n$ is a rational number:

$$ B(m, n) = \frac{(m - 1)!(n - 1)!}{(m + n - 1)!} $$

Thus, from above we have 

$$ \frac{d}{n k} \frac{n!}{(n - k!)(k - 1)!} \space p ^ {k} q^{n-k}
= \frac{d}{n k} B(k, n + 1 - k)^{-1}\space p ^ {k} q^{n-k}
$$

Since `scipy` has [log beta](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaln.html) functions built in, we can decompose the right hand side with a logorithm. 

## Testing

[Evan's web app](https://www.evanmiller.org/ab-testing/sequential.html) only 
tests positive effects in the treatment group with equal assignment. 

This package uses unit tests seeded with outputs from Evan's app.

To test, simply run `pytest` in the command line.
