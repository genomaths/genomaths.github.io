# Information Hellinger and J-divergence of four Arabidopsis Mutants

## The J-divergence 


Given a the methylation levels from two individuals at a given
cytosine site, this function computes the J information divergence
(_JD_) between methylation levels. The motivation to introduce
_JD_ in Methyl-IT is founded on:

1. It is a symmetrised form of [Kullback–Leibler divergence:](https://is.gd/oS8Bwv) 
<img src="https://render.githubusercontent.com/render/math?math=D_{KL}(P || Q)">. 
Kullback and Leibler themselves actually defined the divergence as: 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=JD=D_{KL}(P || Q) %2B D_{KL}(Q || P)"></div>

which is symmetric and nonnegative, where the probability distributions
_P_ and _Q_ are defined on the same probability space
(see reference (1) and [Wikipedia](https://is.gd/oS8Bwv)).

2. In general, _JD_ is highly correlated with Hellinger divergence, 
which is the main divergence currently used in Methyl-IT (see examples for function 
[estimateDivergence](https://genomaths.github.io/methylit/reference/estimateDivergence.html).

3. By construction, the unit of measurement of _JD_ is given in units of bit of information, 
which set the basis for further information-thermodynamics analyses.



