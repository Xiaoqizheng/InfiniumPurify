\name{InfiniumDM}
\alias{InfiniumDM}
\title{
Differentially Methylation Calling accounting for tumor purity.
}
\description{
Infer differentially methylated CpG sites with the consideration of tumor purities.
}
\usage{
InfiniumDM(mode,normal.data,tumor.data,purity,threshold,
           tumor1.data,tumor2.data,
           tumor1.purity,tumor2.purity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mode}{Should be either 'tumorVSnormal' or 'tumorVStumor', indicating DM calling between tumor and normal samples or two groups of tumor samples. Default is 'tumorVStumor'}
  \item{normal.data}{Numeric matrix of beta values for normal samples. The rownames of normal.data
  should be probe names of Infinium 450k array, and colnames should be names of normal samples.}
  \item{tumor.data}{Numeric matrix of beta values for normal samples. The rownames of tumor.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples.}
  \item{purity}{Purities for tumor samples. Could be estimated by \link{getPurity}, or user specified 
  purities from other tools.}
  \item{threshold}{Probability threshold in control-free DM calling. Default is 0.1.}
  \item{tumor1.data}{Numeric matrix of beta values for tumor samples in subtype 1. The rownames of tumor1.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples in subtype 1.}
  \item{tumor2.data}{Numeric matrix of beta values for tumor samples in subtype 2. The rownames of tumor2.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples in subtype 2.}
  \item{tumor1.purity}{Purities for tumor samples in subtype 1.}
  \item{tumor2.purity}{Purities for tumor samples in subtype 2.}
}
\details{
When 'mode' is set as tumorVSnormal, the function performs differential methylation analysis between tumor and normal samples,
otherwise, the function performs differential methylation analysis between two groups of tumor samples. For tumor versus normal comparison, if normal.data is provided, the function tests each CpG site for differential methylation between tumor and normal samples with the consideration of tumor purities by a generalized linear regression, if normal.data is not provided, the function computes posterior probability to rank CpG sites. For tumor and tumor comparison, if normal.data is provided, the function takes beta values of normal and two groups of tumor samples as input data, otherwise, only beta values of two groups of tumor samples are used as input data. The function includes tumor purities as an experimental design factor in the linear regression system, and computes the Wald test statistics to determine the DM status of all CpG sites.
}
\value{
A data frame of statistics, p-values and q-values for all CpG sites.
}
\references{
X. Zheng, N. Zhang, H.J. Wu and H. Wu, Estimating and accounting for tumor purity in the analysis of DNA methylation data from cancer studies. \emph{Genome biology} \bold{18}, 17. 

W. Zhang, Z. Li, N. Wei, H. Wu and X. Zheng, Detection of differentially methylated CpG sites between tumor samples with uneven tumor purities. in submission.
}

\author{
Xiaoqi Zheng \email{xqzheng@shnu.edu.cn}.
}


\examples{

## load example data
data(beta.emp)
normal.data <- beta.emp[,1:21]
tumor.data <- beta.emp[,22:61]

## estimate tumor purity
purity <- getPurity(tumor.data = tumor.data,normal.data = normal.data)

## DM calling between tumor and normal samples with normal controls
DMC_ctl = InfiniumDM(mode = "tumorVSnormal",normal.data=normal.data,tumor.data=tumor.data,
purity=purity,threshold=0.1)

## DM calling between tumor and normal samples without normal controls
DMC = InfiniumDM(mode = "tumorVSnormal",normal.data=NULL,tumor.data=tumor.data,
purity=purity,threshold=0.1)

## DM calling between two groups of tumor samples
tumor1.data = tumor.data[,1:20]
tumor2.data = tumor.data[,21:40]
tumor1.purity = purity[1:20]
tumor2.purity = purity[21:40]

## DM calling between two groups of tumor samples without normal controls
DMC = InfiniumDM(mode = "tumorVStumor",normal.data=NULL,
                 tumor1.data=tumor1.data,tumor2.data=tumor2.data,
                 tumor1.purity=tumor1.purity,tumor2.purity=tumor2.purity)

## DM calling between two groups of tumor samples with normal controls
DMC_ctl = InfiniumDM(mode = "tumorVStumor",normal.data=normal.data,
                     tumor1.data=tumor1.data,tumor2.data=tumor2.data,
                     tumor1.purity=tumor1.purity,tumor2.purity=tumor2.purity)

}