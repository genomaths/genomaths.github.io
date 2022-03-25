# Information Hellinger and J-divergence of Methylation levels in four Arabidopsis Mutant Samples
The data set includes the Hellinger and J-divergence of Methylation levels in four Arabidopsis Mutant Sample analyzed
in (<a href="#1">1</a>): 
1. MSH1
2. MET1
3. DDM1

The methylation levels of mutant MET1 and DDM1 were reported in (<a href="#2">2</a>)


## Hellinger divergence

The methylation level <img src="https://render.githubusercontent.com/render/math?math=p_{ij}"> for an individual _i_ at
cytosine site _j_ corresponds to a probability vector 
<img src="https://render.githubusercontent.com/render/math?math=p^{ij}=(p_{ij},1 - p_{ij})">. 
Then, the information divergence between methylation levels
The methylation level <img src="https://render.githubusercontent.com/render/math?math=p^1j"> and
<img src="https://render.githubusercontent.com/render/math?math=p^2j"> from individuals 1 and 2 at site _j_
is the divergence between the vectors 
<img src="https://render.githubusercontent.com/render/math?math=p^1j = (p_{1j}, 1 - p_{1j})"> and 
<img src="https://render.githubusercontent.com/render/math?math=p^1j = (p_{2j}, 1 - p_{2j})">.
If the vector of coverage is supplied, then the information divergence is estimated according to the formula:

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=HD = \dfrac{2 (n_1 %2B 1) (n_2  %2B 1)}{n_1  %2B n_2 %2B 2} ((\sqrt(p_{1j}) - \sqrt(p_{2j}))^2  %2B (\sqrt(1 - p_{1j}) - \sqrt(1 - p_{2j}))^2)"></div>

This formula corresponds to Hellinger divergence as given in the first
formula from Theorem 1 from reference 1. Otherwise:

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=HD = (\sqrt(p_{1j}) - \sqrt(p_{2j}))^2  %2B (\sqrt(1 - p_{1j}) - \sqrt(1 - p_{2j}))^2"></div>

which is formula applied to the current set of samples.

## The J-divergence 

Given a the methylation levels from two individuals at a given
cytosine site, this function computes the J information divergence
(_JD_) between methylation levels. The motivation to introduce
_JD_ in Methyl-IT is founded on:

1. It is a symmetrised form of [Kullback–Leibler divergence:](https://is.gd/oS8Bwv) 
<img src="https://render.githubusercontent.com/render/math?math=D_{KL}(P || Q)">. 
Kullback and Leibler themselves actually defined the divergence as (<a href="#4">4</a>): 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=JD=D_{KL}(P || Q) %2B D_{KL}(Q || P)"></div>

which is symmetric and nonnegative, where the probability distributions
_P_ and _Q_ are defined on the same probability space
(see reference (1) and [Wikipedia](https://is.gd/oS8Bwv)).

2. In general, _JD_ is highly correlated with Hellinger divergence, 
which is the main divergence currently used in Methyl-IT (see examples for function 
[estimateDivergence](https://genomaths.github.io/methylit/reference/estimateDivergence.html).

3. By construction, the unit of measurement of _JD_ is given in units of bit of information, 
which set the basis for further information-thermodynamics analyses.


# References
1. <a name="1"></a>Sanchez R, Mackenzie SA. Genome-Wide Discriminatory Information Patterns of Cytosine DNA Methylation. International Journal of Molecular Sciences. 2016; 17(6):938. <a href="https://doi.org/10.3390/ijms17060938">DOI: doi.org/10.3390/ijms17060938</a>. 
2. <a name="2"></a>Stroud, H., et all (2012). DNA methyltransferases are required to induce heterochromatic re-replication in Arabidopsis. 
<a href="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002808">PLoS genetics</a>, 8(7), e1002808.
4. Basu  A., Mandal  A., Pardo L (2010) Hypothesis testing for two discrete populations based on the Hellinger distance. Stat Probab Lett 80: 206-214. 
<a href="https://doi.org/10.1016/j.spl.2009.10.008">DOI: doi.org/10.1016/j.spl.2009.10.008</a>.
6. <a name="2"></a>Kullback, S.; Leibler, R.A. (1951). "On information and sufficiency". Annals of Mathematical Statistics. 22 (1): 79–86. <a href="https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-22/issue-1/On-Information-and-Sufficiency/10.1214/aoms/1177729694.full#:~:text=10.1214/aoms/1177729694">doi:10.1214/aoms/1177729694</a>.


