# InfiniumPurify
R package for DNA methylation analysis
The proposition of cancer cells in a tumor sample, named as tumor purity, is an intrinsic factor of tumor samples and has potentially great influence in variety of analyses including differential methylation, subclonal deconvolution and subtype clustering. InfiniumPurify is an integrated R pa ckage for est imatin g and accoun ting for tum or puri ty based on DNA methylation Infinium 450 k array data. InfiniumPurify has three main functions getPurity, InfiniumDMC and InfiniumClust, which could infer tumor purity, differential methylation analysis and tumor sample c luster accounting for estimated or user-provided tumor purities, respectively. The InfiniumPurify package provides a comprehensive analysis of tumor purity in cancer methylation research.

How to install?
1. Install the devtools package if needed.
install.packages("devtools")
2. Load the devtools package.
library(devtools)
3. Install InfiniumPurify from GitHub. 
install_github("Xiaoqizheng/InfiniumPurify")

