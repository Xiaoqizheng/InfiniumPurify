\name{getPurity}
\alias{getPurity}
\title{
Estimate the tumor purity for 450K methylation data 
}
\description{
Estimate the percentage of tumor cells in cancer samples which are mixtures of cancer and normal cells.added a sentence
}
\usage{
getPurity(tumor.data,normal.data = NULL,tumor.type = NULL)
}
\arguments{
  \item{tumor.data}{numeric vector/matrix of beta values for tumor samlpes. The names/rownames of tumor.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples.}
  \item{normal.data}{numeric matrix of beta values for normal samlpes. The rownames of normal.data
  should be probe names of Infinium 450k array, and colnames should be names of normal samples.}
  \item{tumor.type}{cancer type (in abbreviation) of tumor and normal samlpes. Options are "LUAD", "BRCA" and so
  on. See \code{CancerTypeAbbr} for detail.}
}
\details{
Arguments normal.data and tumor.type could be null. If either the number of tumor samples or number of normal smaples is less than 20, the tumor.type argument should be specified according to \code{CancerTypeAbbr}. If the numbers of tumor and normal samples are both more than 20, tumor.type could be null. In such case, \link{getPurity} first identify 1000 iDMCs by Wilcox rank-sum test, then tumor purity for each sample is estimated as the density mode of adjusted methylation levels of iDMCs. 
}
\value{
A vector of tumor purities for each tumor sample. 
}
\references{
N. Zhang, H.J. Wu, W. Zhang, J. Wang, H. Wu and X. Zheng (2015) Predicting tumor purity from methylation
microarray data. \emph{Bioinformatics} \bold{31(21)}, 3401-3405. 

X. Zheng, N. Zhang, H.J. Wu and H. Wu, Estimating and accounting for tumor purity in the analysis of DNA methylation data from cancer studies. \emph{Genome biology} \bold{18}, 17. 
}
\author{
Xiaoqi Zheng \email{xqzheng@shnu.edu.cn}.
}
\examples{

## load example data
data(beta.emp)

normal.data <- beta.emp[,1:21]
tumor.data <- beta.emp[,22:61]

## call purity for single tumor sample 
purity <- getPurity(tumor.data = tumor.data[,1],normal.data = NULL,tumor.type= "LUAD")

## call purity for less than 20 tumor samples
purity <- getPurity(tumor.data = tumor.data[,1:10],normal.data = NULL,tumor.type= "LUAD")

## call purity for more than 20 tumor samples with matched normal samples
purity <- getPurity(tumor.data = tumor.data[,1:40],normal.data = normal.data)
}