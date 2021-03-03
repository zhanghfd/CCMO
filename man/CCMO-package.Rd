\name{CCMO-package}
\alias{CCMO-package}
\docType{package}
\title{
Powerful statistical methods for disease-SNP association analysis using case-control mother-offspring pair data.
}
\description{
Two statistical methods (one without parent-of-origin effect and the other with parent-of-origin effect) are implemented in this package, which are designed for detecting SNPs associated with offspring disease using the genotype data from mother-offspring pairs by incorporting maternal covariates. The first method allows for main SNP effects of mother and offspring, the interaction between maternal and offspring genotypes, and the interactions between the SNP genotypes and maternal covariates. The disease prevalence formula should be specified by users. In this method, Hardy-Weiberg equilibrium can be incorprated and mode of inheritance (recessive, additive, and dominant) and be specified. The maternal genotype can be specified to be either dependent (through the so called daLOG model) or independent of the maternal covariates. The second method further incorporates parent-of-origin effects (POEs), which utilizes multi-locus genotypes to boost POE analysis effeciency.
}
\references{
Zhang H, Mukherjee B, Arthur V, Hu G, Hochner H, Chen J (2020). An Efficient and Computationally Robust Statistical Method for Analyzing Case-Control Mother-Offspring Pair Genetic Association Studies. Annals of Applied Statistics 14: 560–-584.

Zhang K, Zhang H, Hochner H, Chen J (2021). Covariate Adjusted Inference of Parent-of-Origin Effects Using Case-Control Mother-Child Paired Multi-Locus Genotype Data. Manuscript.
}

\author{
Hong Zhang, Jinbo Chen, Kai Zhang

Maintainer: Hong Zhang <zhangh@ustc.edu.cn>
}
\keyword{ package }
