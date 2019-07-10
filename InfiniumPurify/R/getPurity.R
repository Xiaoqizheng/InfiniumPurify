
getPurity <- function(tumor.data,normal.data = NULL,tumor.type = NULL)
{
  purity <- c()
  data(iDMC, envir=environment())
  
  if(missing(tumor.data)){
    stop("'tumor.data' is required.")
  }
  
  # how many tumor and normal samples do we have?
  if (is.vector(tumor.data)){
    n.tumor <- 1
  } else {
    n.tumor <- ncol(tumor.data)
  }
  
  if (is.null(normal.data)){
    n.normal <- 0
  } else {
    n.normal <- ncol(normal.data)
  }
  
  if (is.vector(tumor.data)){
    # one tumor sample, load our iDMC data
    if (is.null(tumor.type)){
      stop("should specify 'tumor.type' if having only one tumor sample.")
    }
    tumor.type <- toupper(tumor.type)
    #load("/home/users/zhengxq/dna_methylation/Rpackage_infiniumPurify/data/iDMC.RData")
    if ( ! any(names(iDMC) == tumor.type)) {
      ## tumor.type is not in iDMC
      stop(paste0(tumor.type," is not found! See CancerTypeAbbr() for detail."))
    }
    
    tumor.CpGs <- names(tumor.data)
    comm.CpGs <- intersect(tumor.CpGs,names(iDMC[[tumor.type]]))
    
    idmc.dat <- data.frame(tumor.data[comm.CpGs])
    idmc.dat$hyper <- iDMC[[tumor.type]][comm.CpGs]
    
    cat("Calculating tumor purity ...\n")
    beta.adj <- c(idmc.dat[idmc.dat$hyper == TRUE,1],1-idmc.dat[idmc.dat$hyper == FALSE,1])
    purity <- .get_peak(beta.adj)
  }
  else if (n.tumor < 20 | n.normal < 20) {
    # load our iDMC data
    cat(paste0("The number of tumor samples is: ",ncol(tumor.data),"\n"))
    tumor.sample <- colnames(tumor.data)
    tumor.CpGs <- rownames(tumor.data)
    
    if (is.null(tumor.type)){
      stop("should specify 'tumor.type' if less than 20 tumor samples")
    }
    #load("/home/users/zhengxq/dna_methylation/Rpackage_infiniumPurify/data/iDMC.RData")
    if ( ! any(names(iDMC) == toupper(tumor.type))) {
      ## tumor.type is not in iDMC
      stop(paste0(tumor.type," is not found! See abbr.txt for detail."))
    }
    
    comm.CpGs <- intersect(tumor.CpGs,names(iDMC[[tumor.type]]))

    idmc.dat <- data.frame(tumor.data[comm.CpGs,])
    colnames(idmc.dat) <- colnames(tumor.data)
    idmc.dat$hyper <- iDMC[[tumor.type]][comm.CpGs]
    
    cat("Calculating tumor purity ...\n")
    for(t in tumor.sample){
      beta.adj <- c(idmc.dat[idmc.dat$hyper == TRUE,t],1-idmc.dat[idmc.dat$hyper == FALSE,t])
      pu <- .get_peak(beta.adj)
      purity[t] <- pu
    }
  }
  
  else{ 
    # infer iDMC data by rank-sum test 
    cat(paste0("The number of tumor samples is: ",ncol(tumor.data),"\n"))

    normal.sample <- colnames(normal.data)
    normal.CpGs <- rownames(normal.data)
    
    tumor.sample <- colnames(tumor.data)
    tumor.CpGs <- rownames(tumor.data)
    
    comm.CpGs <- intersect(normal.CpGs,tumor.CpGs)
    
    cat("Getting iDMCs ...\n")
    iDMC <- get_iDMC(tumor.data[comm.CpGs,],normal.data[comm.CpGs,])
    
    idmc.dat <- data.frame(cbind(tumor.data[iDMC,tumor.sample],normal.data[iDMC,normal.sample]))
    colnames(idmc.dat) <- c(colnames(tumor.data),colnames(normal.data))
    idmc.dat$hyper <- rowMeans(idmc.dat[,tumor.sample],na.rm = TRUE) > rowMeans(idmc.dat[,normal.sample],na.rm = TRUE)
    
    cat("Calculating tumor purity ...\n")

    for(t in tumor.sample){
      beta.adj <- c(idmc.dat[idmc.dat$hyper == TRUE,t],1-idmc.dat[idmc.dat$hyper == FALSE,t])
      pu <- .get_peak(beta.adj)
      purity[t] <- pu
    }
  }

  return(purity)
}


.get_peak <- function(dat){
  d <- density(dat,na.rm = TRUE,kernel = "gaussian")
  d$x[which.max(d$y)]
}


.myVar <- function(x){
  var(x,na.rm = TRUE)
}

get_iDMC <- function(tumor.data,normal.data){
  
  all.dat <- cbind(tumor.data,normal.data)
  tumor.sample <- colnames(tumor.data)
  normal.sample <- colnames(normal.data)
  
  .myRanksum <- function(x){
    ## only works within function get_iDMC 
    if (length(na.omit(x[tumor.sample])) == 0 | length(na.omit(x[normal.sample])) == 0){
      return(NA)
    }
    else{
      pval <- suppressWarnings(wilcox.test(as.numeric(x[tumor.sample]),as.numeric(x[normal.sample]))$p.value)
      return(pval)
    }
  }

  ranksum.pval <- apply(all.dat,1,.myRanksum)
  tumor.var <- apply(all.dat[,tumor.sample],1,.myVar)
  
  out <- data.frame(ranksum.pval)
  out <- cbind(out,tumor.var)
  
  cDMC <- out[out$tumor.var >= 0.005,]
  idx <- order(cDMC$ranksum.pval,decreasing=FALSE)[1:1000]
  iDMC <- rownames(cDMC)[idx]
  
  return(iDMC)
}

CancerTypeAbbr <- function(){
  abbr <- NULL
  data(abbr, envir=environment())
  print(abbr)
}
