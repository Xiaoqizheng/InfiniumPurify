\name{InfiniumClust}
\alias{InfiniumClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tumor sample clustering from Infinium 450k array data  
}
\description{
Clustering of tumor samples into subtypes accounting for tumor purity. 
}
\usage{
InfiniumClust(tumor.data, purity, K, maxiter = 100, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{tumor.data}{numeric matrix of beta values for tumor samlpes. The rownames of tumor.data
  should be probe names of Infinium 450k array, and colnames should be names of tumor samples.}
  \item{purity}{purities for tumor samples. Could be estimated by \code{getPurity}, or user specified 
  purities from other tools.}
  \item{K}{the number of clusters.}
  \item{maxiter}{the maximum number of iterations allowed. Default is 100.}
  \item{tol}{tolerance for convergence of EM iterations. Default is 0.001.}
}
\details{
An EM based statistical method for subtype classification based on DNA methylation data, while adjusting for tumor purity. 
}
\value{
InfiniumClust returns a list consisting oflikelihood \code{tol.ll} and membership matrix \code{Z}.
\item{tol.ll}{total log-likelihood of converged EM algorithm.}
\item{Z}{the membership matrix, where row corresponds to tumor samples and column corresponds to K clusters.}

}
\references{
W. Zhang, H. Feng, H. Wu and X. Zheng. Tumor purity improves cancer subtype classification from DNA methylation data.
\emph{Bioinformatics} \bold{33(17)}, 2651-2657. 
}
\author{
Xiaoqi Zheng \email{xqzheng@shnu.edu.cn} and Hao Wu \email{hao.wu@emory.edu}
}

\examples{
## load example data
data(beta.emp)
normal.data <- beta.emp[,1:21]
tumor.data <- beta.emp[,22:31]

## estimate tumor purity
purity <- getPurity(tumor.data = tumor.data,tumor.type= "LUAD")

## cluster tumor samples accounting for tumor purity
out <- InfiniumClust(tumor.data,purity,K=3, maxiter=5, tol=0.001)
}