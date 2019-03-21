---
title: "Methyl-IT R Package"
subtitle: <h1> Manual for Methyl-IT version available at https://github.com/genomaths/MethylIT </h1>
author: |
 | Robersy Sanchez
 | rus547@psu.edu
 | ORCID: orcid.org/0000-0002-5246-1453
 | Mackenzie's lab    
 |
 |
 | Department of Biology and Plant Science. 
 | Pennsylvania State University, University Park, PA 16802
 | Description: Methylation Analysis based on Information Theory.
date: "21 March 2019"
fontsize: 11pt
fontfamily: "serif"
output:
  rmarkdown::html_document: 
    toc: true
    toc_depth: 3
    toc_float: 
      collapsed: false
      smooth_scroll: true
    number_sections: false
    highlight: tango
    theme: united
    geometry: margin=1in
    keep_md: yes

---

<style type="text/css">

body{ /* Normal  */
      font-size: 18px;
      font-family: "Times New Roman", Times, serif;
  }
td {  /* Table  */
  font-size: 8px;
}

h1.title {
  font-size: 38px;
  font-family: "Times New Roman", Times, serif;
  color: DarkRed;
  .subTitle {
  font-size: 24px;
  font-family: "Times New Roman", Times, serif;
  color: DarkRed;
  }
}

h1 { /* Header 1 */
  font-size: 28px;
  font-family: "Times New Roman", Times, serif;
  color: DarkBlue;
}
h2 { /* Header 2 */
    font-size: 22px;
    color: DarkBlue;
    font-family: "Times New Roman", Times, serif;
}
h3 { /* Header 3 */
   font-size: 18px;
   color: DarkBlue;
   font-family: "Times New Roman", Times, serif;
}
code.r{ /* Code block */
    font-size: 12px;
}
pre { /* Code block - determines code spacing between lines */
    font-size: 14px;
}
</style>


# Methyl-IT documented functions:

## [AICmodel](https://genomaths.github.io/MethylIT_HTML_Manual/AICmodel.html)  
## [BICmodel](https://genomaths.github.io/MethylIT_HTML_Manual/BICmodel.html)                                    
## [countTest](https://genomaths.github.io/MethylIT_HTML_Manual/countTest.html)                                  
## [countTest2](https://genomaths.github.io/MethylIT_HTML_Manual/countTest2.html)                                
## [estimateCutPoint](https://genomaths.github.io/MethylIT_HTML_Manual/estimateCutPoint.html)                    
## [estimateDivergence](https://genomaths.github.io/MethylIT_HTML_Manual/estimateDivergence.html)                
## [estimateECDF](https://genomaths.github.io/MethylIT_HTML_Manual/estimateECDF.html)                                                    
## [estimateHellingerDiv](https://genomaths.github.io/MethylIT_HTML_Manual/estimateHellingerDiv.html)            
## [evaluateDIMPclass](https://genomaths.github.io/MethylIT_HTML_Manual/evaluateDIMPclass.html)                  
## [filterByCoverage](https://genomaths.github.io/MethylIT_HTML_Manual/filterByCoverage.html)                    
## [filterGRange](https://genomaths.github.io/MethylIT_HTML_Manual/filterGRange.html)                            
## [FisherTest](https://genomaths.github.io/MethylIT_HTML_Manual/FisherTest.html)                                
## [fitGammaDist](https://genomaths.github.io/MethylIT_HTML_Manual/fitGammaDist.html)                            
## [fitGGammaDist](https://genomaths.github.io/MethylIT_HTML_Manual/fitGGammaDist.html)                          
## [getDIMPatGenes](https://genomaths.github.io/MethylIT_HTML_Manual/getDIMPatGenes.html)                        
## [getGEOSuppFiles](https://genomaths.github.io/MethylIT_HTML_Manual/getGEOSuppFiles.html)                      
## [getPotentialDIMP](https://genomaths.github.io/MethylIT_HTML_Manual/getPotentialDIMP.html)                    
## [ggamma](https://genomaths.github.io/MethylIT_HTML_Manual/ggamma.html)                                        
## [lapply](https://genomaths.github.io/MethylIT_HTML_Manual/lapply.html)                                        
## [MethylIT](https://genomaths.github.io/MethylIT_HTML_Manual/MethylIT.html)                                    
## [nonlinearFitDist](https://genomaths.github.io/MethylIT_HTML_Manual/nonlinearFitDist.html)                    
## [pcaLDA](https://genomaths.github.io/MethylIT_HTML_Manual/pcaLDA.html)                                        
## [pcaLogisticR](https://genomaths.github.io/MethylIT_HTML_Manual/pcaLogisticR.html)                            
## [pcaQDA](https://genomaths.github.io/MethylIT_HTML_Manual/pcaQDA.html)                                        
## [poolFromGRlist](https://genomaths.github.io/MethylIT_HTML_Manual/poolFromGRlist.html)                        
## [predict.LogisticR](https://genomaths.github.io/MethylIT_HTML_Manual/predict.LogisticR.html)                  
## [print.CutPoint](https://genomaths.github.io/MethylIT_HTML_Manual/print.CutPoint.html)                        
## [pweibull3P](https://genomaths.github.io/MethylIT_HTML_Manual/pweibull3P.html)                                
## [readCounts2GRangesList](https://genomaths.github.io/MethylIT_HTML_Manual/readCounts2GRangesList.html)        
## [selectDIMP](https://genomaths.github.io/MethylIT_HTML_Manual/selectDIMP.html)                                
## [sortBySeqnameAndStart](https://genomaths.github.io/MethylIT_HTML_Manual/sortBySeqnameAndStart.html)          
## [uniqueGRanges](https://genomaths.github.io/MethylIT_HTML_Manual/uniqueGRanges.html)                          
## [uniqueGRfilterByCov](https://genomaths.github.io/MethylIT_HTML_Manual/uniqueGRfilterByCov.html)              
## [unlist](https://genomaths.github.io/MethylIT_HTML_Manual/unlist.html)                                        
## [weibull3P](https://genomaths.github.io/MethylIT_HTML_Manual/weibull3P.html) 


