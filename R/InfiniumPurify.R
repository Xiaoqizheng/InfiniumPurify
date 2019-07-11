
InfiniumPurify <- function(tumor.data,normal.data,purity)
{
  ### get purified tumor methylome from tumor.data, normal.data and purity. 
  
  if(missing(tumor.data) | missing(normal.data) | missing(purity)){
    stop("'tumor.data', 'normal.data' and 'purity' are required.")
  }
  
  probes = intersect(rownames(tumor.data),rownames(normal.data))
  tumor.sample = intersect(colnames(tumor.data),names(purity))
  normal.sample = colnames(normal.data)
  purity = purity[tumor.sample]
  
  if(length(normal.sample) < 20 | length(tumor.sample) < 20){
    stop("tumor and normal samples should be more than 20!")
  }
  
  .get_corrBeta <- function(input){
    x=as.numeric(input[tumor.sample])
    y=as.numeric(input[normal.sample])
    type = c(rep("Tumor",length(x)), rep( "Normal", length(y)))
    data = data.frame(beta=c(x,y), type=type)

    ## fit a regression on arcsin transformed values
    Y = .myasin(data$beta)
    X = c(1-purity, rep(0,length(y)))
    
    fit=lm(Y~X)
    tmp = resid(fit)+coef(fit)[1]
    ## get predicted beta values
    beta.pred = (sin(tmp)+1)/2
    beta.pred
  }
  
  ## correct tumor methylomes for each probe
  all.data = cbind(tumor.data[probes,tumor.sample],normal.data[probes,])
  probes.rmna = probes[rowSums(is.na(all.data)) == 0]

  all.data.corr = t(apply(all.data[probes.rmna,],1,.get_corrBeta))
  tumor.data.corr = all.data.corr[,1:length(tumor.sample)]
  colnames(tumor.data.corr) = tumor.sample
  
  tumor.data.corr
}

.myasin <- function(x) asin(2*x-1)
