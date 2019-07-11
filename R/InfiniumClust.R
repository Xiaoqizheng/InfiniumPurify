
#######################################################
## run EM algorithm, now works for K clusters (user-defined).
#######################################################

InfiniumClust = function(tumor.data, purity, K=3, maxiter=100, tol=1e-3)
{
  # select top 100 cpg site by Variance
  tumor.data = na.omit(as.matrix(tumor.data))
  vars = rowVars(tumor.data)
  names(vars) = rownames(tumor.data)
  dmp = names(sort(vars,decreasing=TRUE)[1:100])
  sample = intersect(colnames(tumor.data),names(purity))
  tumor.data = tumor.data[dmp,sample]
  purity = as.vector(purity[sample])
  
  # initial parameters
  par.intial = .initializeParameter(tumor.data,K)
  p = par.intial$p 
  mu = par.intial$mu
  mu.n = par.intial$mu.n
  sd = par.intial$sd
  sd.n = par.intial$sd.n 
    
  diff=1
  iter=0
  N_cpg = nrow(tumor.data)
  N_sample = ncol(tumor.data)
  
  
  ## Z is the class membership, a matrix of N_sample*K  
  Z = matrix(0, nrow = N_sample, ncol = K)

  while (diff>tol & iter<maxiter) {

    tol.ll2 = 0
    
    for(j in 1:nrow(Z)){
      tmp1 <- rep(0,K)
      for(s in 1:K){
        ## calculate each component of denominator in log scale
        tmp1[s] <- log(p[s])+sum(dnorm(tumor.data[,j], mean = purity[j]*mu[,s]+(1-purity[j])*mu.n,
                                       sd = sqrt(purity[j]^2*sd[,s]^2+(1-purity[j])^2*sd.n^2), log = T))
      }
      ## sum up in log scale
      tmp1.sum.log <- .Rsumlog(tmp1)
      tol.ll2 <- tol.ll2 + tmp1.sum.log
      Z[j,] = exp(tmp1 - tmp1.sum.log)
    }
    
    
    ## M-step:
    ## 1, update p (including p1, p2 and p3)
    p.new = colSums(Z)/sum(Z)
    
    ## 2, update mu.n and mu
    ## start updating mu.n
    mu.n.new <- unlist(lapply(1:N_cpg,function(i){
      tmp1 = 0
      tmp2 = 0
      for(j in 1:length(purity)){
        for(k in 1:K){
          tmp1 = tmp1 + Z[j,k]*(1-purity[j])*(tumor.data[i,j]-purity[j]*mu[i,k])/(purity[j]^2*sd[i,k]^2+(1-purity[j])^2*sd.n[i]^2)
          tmp2 = tmp2 + Z[j,k]*(1-purity[j])^2/(purity[j]^2*sd[i,k]^2+(1-purity[j])^2*sd.n[i]^2)
        }
      }
      return(tmp1/tmp2)     
    }))
    ## end of updating mu.n
    
    ## start updating mu
    mu.new <- t(sapply(1:N_cpg,function(i){
      mu.new.row = c()
      for(k in 1:ncol(Z)){
        tmp1 <- sum(Z[,k]*purity*tumor.data[i,]/(purity^2*sd[i,k]^2+(1-purity)^2*sd.n[i]^2))-
          sum(Z[,k]*purity*(1-purity)/(purity^2*sd[i,k]^2+(1-purity)^2*sd.n[i]^2))*mu.n.new[i]
        tmp2 <- sum(Z[,k]*purity^2/(purity^2*sd[i,k]^2+(1-purity)^2*sd.n[i]^2))
        mu.new.row = append(mu.new.row,tmp1/tmp2)
      }
      return(mu.new.row)
    }))
    ## end of updating mu
    
    
    ## 3, update sd
    sd.new <- t(sapply(1:N_cpg,function(i){
      sd.new.row = c()
      for(k in 1:K){
        func_var <- function(x){
          sum(Z[,k]*(log(p[k])+dnorm(tumor.data[i,], mean = purity*mu[i,k]+(1-purity)*mu.n[i],
                                     sd = sqrt(purity^2*x+(1-purity)^2*sd.n[i]^2), log = T)))
        }
        temp = optimize(func_var, interval = c(0.001,100), maximum=T, tol=tol)
        sd.new.row = append(sd.new.row,sqrt(temp$maximum))
      }
      return(sd.new.row)
    }))
    
    
    ## 4, update sd.n
    sd.n.new <- unlist(lapply(1:N_cpg,function(i){
      func_sd.n <- function(x){
        obj <- 0
        for(j in 1:N_sample){
          for(k in 1:ncol(Z)){
            obj <- obj + Z[j,k]*(log(p[k])+dnorm(tumor.data[i,j], mean = purity[j]*mu[i,k]+(1-purity[j])*mu.n[i],
                                                 sd = sqrt(purity[j]^2*sd[i,k]^2+(1-purity[j])^2*x), log = T))
          }
        }
        return(obj)
      }
      temp = optimize(func_sd.n, interval = c(0.001,100), maximum=T, tol=tol)
      return(sqrt(temp$maximum))
    }))
    
    ## calculate diff to check convergence - I modified a bit
    diff <- sqrt(mean((mu.new-mu)^2+(sd.new-sd)^2+(mu.n.new-mu.n)^2+(sd.n.new-sd.n)^2))
    if(is.nan(diff))
      browser()
    p=p.new;
    mu=mu.new;
    mu.n=mu.n.new;
    sd = sd.new;
    sd.n = sd.n.new;
    
    
    iter=iter+1;
    
    cat("Iter", iter,  ", diff=", diff, "\n")
  }
  
  
  cat("# of clusters chosen:", K, "\n")
  cat("Probability of each cluster = ", p, "\n")
  cat("Total log-likelihood = ", tol.ll2, "\n")
  
  ## should return 1: total log-likelihood; 2: membership matrix of N_sample*K
  res = list(tol.ll=tol.ll2, Z=Z)
  return(res)
}

.Raddlog <- function(a, b) {
  result <- rep(0, length(a))
  idx1 <- a>b+200
  result[idx1] <- a[idx1]
  
  idx2 <- b>a+200
  result[idx2] <- b[idx2]
  
  idx0 <- !(idx1|idx2)
  result[idx0] <- a[idx0] + log1p(exp(b[idx0]-a[idx0]))
  result
}

.Rsumlog <- function(a) {
  s <- a[1]
  for(i in 2:length(a))
    s <- .Raddlog(s, a[i])
  s
}

.initializeParameter <- function(tumor.data,K){
  
  ## obtain initial values for EM
  ## This could be important. I'll obtain that from k-means

  result.kmeans = kmeans(t(tumor.data), K)
  N_cpg = nrow(tumor.data)
  mu = sd = matrix(0,nrow = N_cpg, ncol = K)
  for(i in 1:max(result.kmeans$cluster)) {
    ii = result.kmeans$cluster == i
    
    if(sum(ii) > 1){
      mu[,i] = rowMeans(tumor.data[,ii])
      sd[,i] = sqrt(rowVars(tumor.data[,ii]))
    }
    else{
      mu[,i] = tumor.data[,ii]
      sd[,i] = 0
    }

  }
  p = as.numeric(table(result.kmeans$cluster) / length(result.kmeans$cluster))
  mu.n = rowMeans(tumor.data)
  sd.n = sqrt(rowVars(tumor.data))
  
  return(list(p=p,mu=mu,mu.n=mu.n,sd=sd,sd.n=sd.n))
}



