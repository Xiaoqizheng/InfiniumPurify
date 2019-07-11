\name{InfiniumPurify}
\alias{InfiniumPurify}
\title{
Purify tumor methylomes caused by normal cell contamination.
}
\description{
Deconvolute purified tumor methylomes accounting for tumor purity.
}
\usage{
InfiniumPurify(tumor.data,normal.data,purity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumor.data}{numeric matrix of beta values for tumor samlpes. The rownames of tumor.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples.}
  \item{normal.data}{numeric matrix of beta values for normal samlpes. The rownames of normal.data
  should be probe names of Infinium 450k array, and colnames should be names of normal samples.}
  \item{purity}{purities for tumor samples. Could be estimated by \link{getPurity}, or user specified 
  purities from other tools.}
}
\details{
The function deconvolutes purified tumor methylomes by a linear regression model. 
}
\value{
A matrix of purified beta values for all CpG sites (row) and tumor samples (column).
}
\references{
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

## estimate tumor purity
purity <- getPurity(tumor.data = tumor.data,normal.data = NULL,tumor.type= "LUAD")

## correct tumor methylome by tumor purity
tumor.purified = InfiniumPurify(tumor.data = tumor.data[1:100,],
                                normal.data = normal.data[1:100,],
                                purity = purity)

}