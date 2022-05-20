#  Hellinger and J Information Divergences of Methylation levels in four Arabidopsis Mutant Samples
The data set includes the Hellinger and J-divergence of Methylation levels in four Arabidopsis Mutant Sample analyzed
in (<a href="#1">1</a>): 
1.	Memory 1rt generation
2.	Non-memory 1rt generation
3.	Memory 3rd generation
4.	Col-0 control 3rd generation
5.	Graph control Col-0/Col-0
6.	MSH1 mutant
7.	Quadruple mutant dcl-2-3-4-msh1 (we do not have dcl-2-3-4 from the lab)

The Arabidopsis thaliana methylome datasets (with results reported in Table 1) from bisulfite (methylome) sequencing of msh1 memory line and non-memory (normal looking) sibling plants with isogenic Col-0 wild-type control in Arabidopsis was downloaded from the Gene Expression Omnibus (GEO) Series GSE129303a and GSE118874. Mutant and Grafted raw data is available with serie accession number GSE152570.

The raw data for mutant MET1 was taken from GEO serie GSE122394:
  1. Mutants: GSM3465886, GSM3465887, and GSM3465888
  2. Wildtype: GSM3465878, GSM3465879, GSM3465880, and GSM3465881.

Information divergences RData files are in the folder named 'idiv', which were estimated with MethylIT function [estimateDivergence](https://genomaths.github.io/methylit/reference/estimateDivergence.html). 
The best fitted models are in the folder 'gofs', which were estimated with MethylIT function [gofReport](https://genomaths.github.io/methylit/reference/gofReport.html). 

**Since GitHub does not permit files with more 25MB, these datasets are available at [PSU GitLab](https://git.psu.edu/genomath/datasets/-/tree/main/at_mutants). Hence,
this repository at GItHub is only informative**. You can see and download the data at: <https://git.psu.edu/genomath/datasets/-/tree/main/at_mutants>. 

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
Kullback and Leibler themselves actually defined the divergence as (<a href="#2">2</a>): 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=JD=D_{KL}(P || Q) %2B D_{KL}(Q || P)"></div>

which is symmetric and nonnegative, where the probability distributions
_P_ and _Q_ are defined on the same probability space
(see reference (1) and [Wikipedia](https://is.gd/oS8Bwv)).

2. In general, _JD_ is highly correlated with Hellinger divergence, 
which is the main divergence currently used in Methyl-IT (see examples for function 
[estimateDivergence](https://genomaths.github.io/methylit/reference/estimateDivergence.html).

3. By construction, the unit of measurement of _JD_ is given in units of bit of information, 
which set the basis for further information-thermodynamics analyses.     

Methylation levels <img src="https://render.githubusercontent.com/render/math?math=p_{ij}"> 
at a given cytosine site _i_ from an individual _j_, lead to the probability vectors
<img src="https://render.githubusercontent.com/render/math?math=p^{ij}=(p_{ij}, 1 - p_{ij})">. Then
, the J-information divergence between the methylation levels <img src="https://render.githubusercontent.com/render/math?math=p_{ij}">
and <img src="https://render.githubusercontent.com/render/math?math=q^i = (q_i, 1 - q_i)"> (used as reference individual),
is given by the expression : 

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=JD(p^{ij}, q^i) = \frac{1}{2}((p_{ij} - q_i) log(\frac{p_{ij}}{q_i}) %2B (q_i - p_{ij}) log(\frac{(1 - p_{ij})}{(1 - q_i)})"></div>   

The statistic with asymptotic Chi-squared distribution is based on the statistic suggested by Kupperman (1957) (<a href="#3">3</a>) for 
<img src="https://render.githubusercontent.com/render/math?math=D_{KL}"> and commented in references (<a href="#4">4-5</a>) . That is:

<div align="center"><img src="https://render.githubusercontent.com/render/math?math=2\frac{(n_{1} %2B 1) (n_{2} %2B 1)}{(n_1 %2B n_2 %2B 2)} JD(p^{ij}, q^i)"></div>

Where <img src="https://render.githubusercontent.com/render/math?math=n_1"> and <img src="https://render.githubusercontent.com/render/math?math=n_2"> 
are the total counts (coverage in the case of methylation) used to compute the probabilities <img src="https://render.githubusercontent.com/render/math?math=p_{ij}"> and
<img src="https://render.githubusercontent.com/render/math?math=q_i">. A basic Bayesian correction is added to prevent zero counts.

# References
1. <a name="1"></a>Kundariya H, Yang X, Morton K, Sanchez R, Axtell MJ, Hutton SF, et al. MSH1-induced heritable enhanced growth vigor through grafting is associated with the RdDM pathway in plants. Nat Commun. 2020;11: 5343. <a href="doi:10.1038/s41467-020-19140-x">doi:10.1038/s41467-020-19140-x</a>.
2. <a name="2"></a>Kullback, S.; Leibler, R.A. (1951). "On information and sufficiency". Annals of Mathematical Statistics. 22 (1): 79–86. <a href="https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-22/issue-1/On-Information-and-Sufficiency/10.1214/aoms/1177729694.full#:~:text=10.1214/aoms/1177729694">doi:10.1214/aoms/1177729694</a>.
3. <a name="3"></a>Kupperman, M., 1957. Further application to information theory to multivariate analysis and statistical inference. Ph.D. Dissertation, George Washington University</a>.
4. <a name="4-5"></a>Salicrú M, Morales D, Menéndez ML, Pardo L. On the applications of divergence type measures in testing statistical hypotheses. Journal of Multivariate Analysis. 1994. pp. 372–391. doi:10.1006/jmva.1994.1068</a>.
5. <a name="4-5"></a>Basu  A., Mandal  A., Pardo L (2010) Hypothesis testing for two discrete populations based on the Hellinger distance. Stat Probab Lett 80: 206-214. 
<a href="https://doi.org/10.1016/j.spl.2009.10.008">DOI: doi.org/10.1016/j.spl.2009.10.008</a>.
