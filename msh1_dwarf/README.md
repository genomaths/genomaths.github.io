# Differentially Methylated Position on _msh1_ TDNA Mutant with respect to its Wildtype Control

This dataset carries the DMPs on _msh1_ TDNA Mutant with respect to its Wildtype Control on WIG files format.
The msh1 T-DNA mutant was obtained from Arabidopsis Biological Resource Center (SAIL_877_F01, stock number CS877617).

These _mash1_ TDNA Mutant samples have _*dwarf*_ phenotype as described in reference ([1](#1)). Dwarf samples labeled 
dmps_dw1, dmps_dw2, and dmps_dw3 are second-generation siblings from a single first-generation msh1 (see Fig. 1 from
[1](#1)); while dwarf samples labeled as dmps_dw_1 and dmps_dw_2 derive from seeds from the second generation (
which is the dwarf material used in reference [2](#2)).


# Methylation analysis
Raw sequencing reads were quality-controlled with FastQC (version 0.11.5), trimmed with TrimGalore! (version 0.4.1) 
and Cutadapt (version 1.15), then aligned to the TAIR10 reference genome using Bismark (version 0.19.0) with bowtie2
(version 2.3.3.1). The deduplicate_bismark function in Bismark with default parameters was used to remove duplicated
reads and reads with coverage >500 were also removed to control PCR bias. Whole-genome bisulfite conversion rate was
computed based on chloroplast genome read counts for every sample, with conversion rate >99% for all samples. DMPs 
were identified using Methyl-IT (version 0.3.2; https://github.com/genomaths/MethylIT) R package as described previously
with some parameters modified for this study. Briefly, cytosine with minimum read coverage of 4 and minimum methylated 
reads of 3 were used. Hellinger Divergence (HD) was calculated with a pool of control (wild type) samples as reference. 
Cytosines with methylation level difference >20% in the treatment vs. reference comparison were selected and further 
filtered by estimating the optimal cutoff for HD based on Youden index to obtain DMPs.

The WIG files only reports DMPs with methylation levels greater than 0.5, which are Bayesian corrected methylation 
levels as described in MethylIT pipeline (see https://genomaths.github.io/methylit/). 

# References
1.<a name="1"></a> Virdi, K., Laurie, J., Xu, YZ. et al. Arabidopsis MSH1 mutation alters the epigenome and produces 
                   heritable changes in plant growth. Nat Commun 6, 6386 (2015). https://doi.org/10.1038/ncomms7386.    
2.<a name="2"></a> Kenchanmane Raju SK, Shao MR, Wamboldt Y, Mackenzie S. Epigenomic plasticity of Arabidopsis msh1 
                   mutants under prolonged cold stress. Plant direct. 2018 Aug;2(8):e00079.
