###################################
## functions for DM calling 
##################################
.myasin <- function(x) asin(2*x-1)

### my R implementation of InfiniumPurify
InfiniumDM = function(mode = "tumorVStumor",
                      normal.data,
                      tumor.data,
                      purity,
                      threshold = 0.1,
                      tumor1.data = NULL,
                      tumor2.data = NULL,
                      tumor1.purity = NULL,
                      tumor2.purity = NULL) {
  if (is.null(mode)) {
    stop("Please specify 'mode', should be either 'tumorVSnormal' or 'tumorVStumor'.")
  }
  if (mode == "tumorVSnormal") {
    cat("DM calling between tumor and normal samples ...\n")
    res = .InfiniumDMC(normal.data, tumor.data, purity, threshold)
  } else if (mode == "tumorVStumor") {
    cat("DM calling between two groups of tumor samples ...\n")
    res = .InfiniumDM1(normal.data,
                       tumor1.data,
                       tumor2.data,
                       tumor1.purity,
                       tumor2.purity)
  } else{
    stop("'mode' should be either 'tumorVSnormal' or 'tumorVStumor'!")
  }
  return(res)
}

## function for DM calling between tumor samples
.InfiniumDM1 = function(normal.data = NULL, tumor1.data, tumor2.data, tumor1.purity, tumor2.purity){
  
  if(is.null(normal.data)){
    # DM calling without normal control samples
    res = .CtrlFreeDM(tumor1.data,tumor2.data,tumor1.purity,tumor2.purity)
  }
  else{
    # DM calling with normal control samples
    res = .WithCtrlDM(normal.data,tumor1.data,tumor2.data,tumor1.purity,tumor2.purity)
  }
  return(res)
}

## function for DM calling between tumor and normal samples
.InfiniumDMC = function(normal.data = NULL,tumor.data,purity,threshold=0.1) {
  
  if (is.null(normal.data)){
    # control-free DM calling
    res = .CtrlFreeDMC(tumor.data,purity,threshold=threshold)
  }
  else{
    # with control DM calling
    res = .WithCtrlDMC(X = normal.data,Y = tumor.data,purity)
  }
  return(res)
}


.WithCtrlDMC = function(X, Y, purity) {
  rawdata = na.omit(cbind(Y,X[rownames(Y),]))
  dat.transform = .myasin(rawdata)
  mat = cbind(1, c(purity[colnames(Y)], rep(0, ncol(X))))
  mat.trans = t(mat)
  invXX = solve(mat.trans%*%mat)
  H = invXX %*% mat.trans
  coefs = H %*% t(dat.transform)
  Ypred = t(mat %*% coefs)
  resi = dat.transform - Ypred
  ## estimate  variances in case/cntl seperately
  ncase = ncol(Y)
  ncntl = ncol(X)
  s2.case = rowSums(resi[,1:ncase]^2) / (ncase - ncol(mat))# + const
  s2.cntl = rowSums(resi[,ncase+1:ncntl]^2) / (ncntl - ncol(mat))# + const
  ## shrink the variances a bit. Use an ad hoc shrinkage procedure
  shrinker = function(vv) {
    tmp = log(vv)
    exp((tmp+mean(tmp))/2)
  }
  s2.case = shrinker(s2.case)
  s2.cntl = shrinker(s2.cntl)
  
  ## restrict that s2.case>s2.cntl
  ix = s2.case<s2.cntl
  s2.case[ix] = s2.cntl[ix] + 0.001
  
  
  H1 = H[,1:ncase]; H2 = H[,ncase+1:ncntl]
  ## compute standard error for coefficient.
  se.new = sqrt((H1%*%t(H1))[2,2]*s2.case + (H2%*%t(H2))[2,2]*s2.cntl)
  ## compute test statistics and p-values
  stats = coefs[2,]/se.new
  df = ncol(X) + ncol(Y) - ncol(mat)
  pval = 2*pt(-abs(stats), df=df)
  out = as.data.frame(stats)
  out$pval = pval
  out$qval = p.adjust(out$pval,method = "BH")
  out = out[order(out$qval),]
  out
}


.CtrlFreeDMC <- function(tumor.data,purity,threshold=0.1){
  ## control free DM calling
  sample.comm = intersect(colnames(tumor.data),names(purity))
  tumor.data = na.omit(tumor.data[,sample.comm])
  design.matrix = cbind(1,purity[sample.comm]) # design matrix
  prob = .get_ctlFree_probability(design.matrix,tumor.data,threshold = threshold)
  prob = sort(prob,decreasing=TRUE)
  out = as.data.frame(prob)
  out
}


.get_ctlFree_probability <- function(X,Y,threshold){
  ###### start estimation procedure
  ## estimate coefficient.
  invXX = solve(t(X) %*% X)
  XY = t(X) %*% t(Y)
  
  beta.hat = invXX %*% XY
  
  est.slope = beta.hat[2,]
  
  ## estimate residual
  nsample = dim(X)[1]
  res = Y - t(X%*% beta.hat)
  sd.CG = sqrt(rowSums(res^2) / (nsample-2))  ## this is estimated residual variance
  
  se.beta = as.numeric(invXX[2,2] * sd.CG) ## this is SE of estimated slope
  
  
  ## test for beta=0. this will reject even when effect size is small
  pval0 = 2*(1-pnorm(abs(est.slope)/se.beta)) ## two-sided p-value for testing the null
  
  ## compute P(|beta|>threshold)
  
  postprob = pnorm(est.slope-threshold, sd=se.beta, lower.tail=TRUE) +
    pnorm(est.slope+threshold, sd=se.beta, lower.tail=FALSE)
  
  return(postprob) ## this P(|beta|>threshold)
}

.WithCtrlDM = function(X, Y1, Y2, purity1, purity2) {
  rawdata = na.omit(cbind(X,Y1,Y2))
  dat.transform = .myasin(rawdata)
  mat = cbind(1,c(rep(0,dim(X)[2]),purity1,purity2),c(rep(0,dim(X)[2]+length(purity1)),purity2))
  mat.trans = t(mat)
  invXX = solve(mat.trans%*%mat)
  H = invXX %*% mat.trans
  coefs = H %*% t(dat.transform)
  Ypred = t(mat %*% coefs)
  resi = dat.transform - Ypred
  ## estimate  variances in case/cntl seperately
  ncntl = ncol(X)
  ncase1 = ncol(Y1)
  ncase2 = ncol(Y2)
  s2.cntl = rowSums(resi[,1:ncntl]^2) / (ncntl - ncol(mat))
  s2.case1 = rowSums(resi[,ncntl+1:ncase1]^2) / (ncase1 - ncol(mat))# + const
  s2.case2 = rowSums(resi[,ncntl+ncase1+1:ncase2]^2) / (ncase2 - ncol(mat))# + const
  ## shrink the variances a bit. Use an ad hoc shrinkage procedure
  shrinker = function(vv) {
    tmp = log(vv)
    exp((tmp+mean(tmp))/2)
  }
  s2.cntl = shrinker(s2.cntl)
  s2.case1 = shrinker(s2.case1)
  s2.case2 = shrinker(s2.case2)
  
  ## restrict that s2.case1>s2.cntl
  ix = s2.case1<s2.cntl
  s2.case1[ix] = s2.cntl[ix] + 0.001
  
  ## restrict that s2.case2>s2.cntl
  ix = s2.case2<s2.cntl
  s2.case2[ix] = s2.cntl[ix] + 0.001  
  
  H0 = H[,1:ncntl]; H1 = H[,ncntl+1:ncase1]; H2 = H[,ncntl+ncase1+1:ncase2]
  
  ## compute standard error for coefficient
  se.new = sqrt((H0%*%t(H0))[3,3]*s2.cntl + (H1%*%t(H1))[3,3]*s2.case1 + (H2%*%t(H2))[3,3]*s2.case2) 
  
  ## compute test statistics and p-values
  stats = coefs[3,]/se.new
  df = ncol(X) + ncol(Y1) + ncol(Y2) - ncol(mat)
  pval = 2*pt(-abs(stats), df=df)
  out = as.data.frame(stats)
  out$pval = pval
  out$qval = p.adjust(out$pval,method = "BH")
  out = out[order(out$qval),]
  out
}

.CtrlFreeDM <- function(Y1, Y2, purity1, purity2){
  rawdata = na.omit(cbind(Y1,Y2))
  dat.transform = .myasin(rawdata)
  mat = cbind(1,c(purity1,purity2),c(rep(0,length(purity1)),purity2))
  mat.trans = t(mat)
  invXX = solve(mat.trans%*%mat)
  H = invXX %*% mat.trans
  coefs = H %*% t(dat.transform)
  Ypred = t(mat %*% coefs)
  resi = dat.transform - Ypred
  ## estimate  variances in case/cntl seperately
  ncase1 = ncol(Y1)
  ncase2 = ncol(Y2)
  s2.case1 = rowSums(resi[,1:ncase1]^2) / (ncase1 - ncol(mat))# + const
  s2.case2 = rowSums(resi[,ncase1+1:ncase2]^2) / (ncase2 - ncol(mat))# + const
  ## shrink the variances a bit. Use an ad hoc shrinkage procedure
  shrinker = function(vv) {
    tmp = log(vv)
    exp((tmp+mean(tmp))/2)
  }
  s2.case1 = shrinker(s2.case1)
  s2.case2 = shrinker(s2.case2)
  
  H1 = H[,1:ncase1]; H2 = H[,ncase1+1:ncase2]
  
  ## compute standard error for coefficient
  se.new = sqrt((H1%*%t(H1))[3,3]*s2.case1 + (H2%*%t(H2))[3,3]*s2.case2) 
  ## compute test statistics and p-values
  stats = coefs[3,]/se.new
  df = ncol(Y1) + ncol(Y2) - ncol(mat)
  pval = 2*pt(-abs(stats), df=df)
  out = as.data.frame(stats)
  out$pval = pval
  out$qval = p.adjust(out$pval,method = "BH")
  out = out[order(out$qval),]
  out
}


