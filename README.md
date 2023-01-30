# freqseq

## Intro

`freqseq` stands for "Frequentist Sequential" and is an implementation of one sided frequentist sequential hypothesis testing in  `python`. Credit for the initial derivation goes to [Evan Miller](https://www.evanmiller.org/sequential-ab-testing.html#notes). I've simply gone over his derivation in more granularity, added calculations dealing with treatment assignment bias, and written a `python` implementation. 


## Problem

Suppose we are running a random experiment to determine the efficacy of some intervention. We would like to end this experiment early if results look promising without "peaking" (peaking incurs bias by effectively testing multiple hypotheses). 

What we want here is a **sequential test**, a test that allows for early stopping if results look promising without incurring a heightened false positive rate. 

## Repo Contents

1. `notes.md` - a derivation of the testing appraoch
2. `freqseq` - a `python` package containing implementation of the testing strategy (to be added)


**Note: this is a working branch and development is under progress, hence some missing code**


