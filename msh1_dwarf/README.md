# Differentially Methylated Positions on _msh1_ TDNA Mutant with respect to its Wildtype Control

This dataset carries the DMPs on _msh1_ TDNA Mutant with respect to its Wildtype Control on WIG files format.
The msh1 T-DNA mutant was obtained from Arabidopsis Biological Resource Center (SAIL_877_F01, stock number CS877617).

These _mash1_ TDNA Mutant samples have _*dwarf*_ phenotype as described in reference ([1](#1)). Dwarf samples labeled 
dmps_dw1, dmps_dw2, and dmps_dw3 are second-generation siblings from a single first-generation msh1 (see Fig. 1 from
([1](#1))); while dwarf samples labeled as dmps_dw_1 and dmps_dw_2 derive from seeds from the second generation (
which is the dwarf material used in reference ([2](#2))).


# Methylation analysis
Raw sequencing reads were quality-controlled with FastQC (version 0.11.5), trimmed with TrimGalore! (version 0.4.1) 
and Cutadapt (version 1.15), then aligned to the TAIR10 reference genome using Bismark (version 0.19.0) with bowtie2
(version 2.3.3.1). The deduplicate_bismark function in Bismark with default parameters was used to remove duplicated
reads and reads with coverage >500 were also removed to control PCR bias. Whole-genome bisulfite conversion rate was
computed based on chloroplast genome read counts for every sample, with conversion rate >99% for all samples. DMPs 
were identified using Methyl-IT (version 0.3.2; https://github.com/genomaths/MethylIT) R package as described previously
with some parameters modified for this study. Briefly, cytosine with minimum read coverage of 4 and minimum methylated 
reads of 3 were used. Hellinger Divergence (HD) was calculated with a pool of control (wild type) samples as reference. 

Only cytosines with methylation levels greater than 0.3 are reported in the WIG files, which are Bayesian corrected methylation 
levels as described in MethylIT pipeline (see https://genomaths.github.io/methylit/). 

Next, depending on the methylation context, cytosines site with the following methylation level difference were considere in 
further downstream analysis:
  1. CG: $\ge 0.34$
  2. CHG: $\ge 0.33$
  3. CHH: $\ge 0.28$
   
The optimal cutoff for HD to obtain DMPs was applied according to the machine-learning approach suggested in 
MethylIT pipeline (also see https://genomaths.github.io/methylit/ and reference ([3](#3)) using MethylIT's
[_estimateCutPoint_](https://genomaths.github.io/methylit/reference/estimateCutPoint.html) function:

```{r cuts}
## Cutpoint estimation for CG methylation context
cut_cg = estimateCutPoint(LR = ps_cg, simple = FALSE,
                          control.names = c( "col1", "col2", "col3",
                                             "col_1", "col_2" ),
                          treatment.names = c( "dw1", "dw2", "dw3",
                                               "dw_1", "dw_2"),
                          column = c(hdiv = TRUE, bay.TV = TRUE,
                                     wprob = TRUE, pos = TRUE),
                          div.col = 9,
                          classifier1 = "pca.logistic",
                          classifier2 = "pca.qda", n.pc = 4,
                          center = TRUE, scale = TRUE,
                          verbose = FALSE)
```
DMGs were estimated with *countTest2* function:

```{r dmgs}
dmgs <- countTest2(ds, num.cores = 4L, minCountPerIndv = 7,
                   maxGrpCV = c(1, 1), Minlog2FC = 1, test = "LRT",
                   CountPerBp = 0.001, verbose = TRUE)
```


# References
1.<a name="1"></a> Virdi, K., Laurie, J., Xu, YZ. et al. Arabidopsis MSH1 mutation alters the epigenome and produces 
                   heritable changes in plant growth. Nat Commun 6, 6386 (2015). https://doi.org/10.1038/ncomms7386.    
2.<a name="2"></a> Kenchanmane Raju SK, Shao MR, Wamboldt Y, Mackenzie S. Epigenomic plasticity of Arabidopsis msh1 
                   mutants under prolonged cold stress. Plant direct. 2018 Aug;2(8):e00079.             
3.<a name="3"></a> Kundariya, H., Yang, X., Morton, K. et al. MSH1-induced heritable enhanced growth vigor through 
                   grafting is associated with the RdDM pathway in plants. Nat Commun 11, 5343 (2020). 
                   https://doi.org/10.1038/s41467-020-19140-x
               
