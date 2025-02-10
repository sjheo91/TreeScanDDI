
# Tree Structure ----------------------------------------------------------

tree_split <- function(tree){
  tree_df <- data.frame(do.call(rbind,strsplit(sort(tree),'_')))
  names(tree_df) <- paste0('lv_',1:ncol(tree_df))
  tree_ind <- list()
  k <- 1
  for(i in 1:ncol(tree_df)){
    u_node <- unique(tree_df[,ncol(tree_df)-i+1])
    for(j in 1:length(u_node)){
      tree_ind[[k]] <- which(tree_df[,ncol(tree_df)-i+1]%in%u_node[j])
      k <- k + 1
    }
  }
  return(list(tree=c(sort(tree),unique(tree_df[,1])), tree_df=tree_df, tree_ind=tree_ind))
}

tree_signal <- function(signal, H0, H1, tree_str, margin_AE, type='additive'){
  
  signal_mat <- matrix(H0, nrow=length(tree_str$tree), ncol=4, byrow=T)
  signal_mat[signal,] <- matrix(H1, nrow=length(signal), ncol=4, byrow=T)
  
  for(i in (nrow(tree_str$tree_df)+1):length(tree_str$tree_ind)){
    if(length(tree_str$tree_ind[[i]])>1){
      signal_mat[i,] <- apply(signal_mat[tree_str$tree_ind[[i]],], 2, 
                              weighted.mean, margin_AE[tree_str$tree_ind[[i]]])
    }else{
      signal_mat[i,] <- signal_mat[tree_str$tree_ind[[i]],]
    }
  }
  
  if(type=='additive'){
    signal_vec <- signal_mat[,1]+signal_mat[,4]-signal_mat[,2]-signal_mat[,3]
  }else{
    signal_vec <- (signal_mat[,1]*signal_mat[,4])/(signal_mat[,2]*signal_mat[,3]) 
  }

  return(signal_vec)
}

# Omega -------------------------------------------------------------------

calc_omega <- function(data, tree_str, tree, alpha=0.5){
  
  tree_ind <- which(tree_str$tree==tree)
  nodes <- tree_str$tree_ind[[tree_ind]]
  margin <- as.vector(apply(data[,c('n111','n101','n011','n001')],2,sum))
  
  f11 <- sum(data$n111[nodes])/margin[1]
  f10 <- sum(data$n101[nodes])/margin[2]
  f01 <- sum(data$n011[nodes])/margin[3]
  f00 <- sum(data$n001[nodes])/margin[4]
  g11 <- 1-1/(max(f00/(1-f00),f10/(1-f10))+max(f00/(1-f00),f01/(1-f01))-f00/(1-f00)+1)
  
  n111 <- sum(data$n111[nodes])
  E111 <- g11*margin[1]
  
  omega <- log2((n111+alpha)/(E111+alpha))
  #omega025_f <- omega-qnorm(0.975)/(log(2)*sqrt(n111))
  omega025 <- log2(qgamma(0.025, shape=n111+alpha, scale=1/(E111+alpha)))
  
  return(data.frame(NodeID=tree_str$tree[tree_ind], n111, E111, RR=n111/E111, omega, omega025))
}

omega_shrinkage <- function(data, tree_str, alpha=0.5){
  do.call(rbind,lapply(tree_str$tree, FUN=calc_omega, data=data, tree_str=tree_str, alpha=alpha))
}

# Chisq -------------------------------------------------------------------

calc_chisq <- function(data, tree_str, tree){
  
  tree_ind <- which(tree_str$tree==tree)
  nodes <- tree_str$tree_ind[[tree_ind]]
  margin <- as.vector(apply(data[,c('n111','n101','n011','n001')],2,sum))
  
  f11 <- sum(data$n111[nodes])/margin[1]
  f10 <- sum(data$n101[nodes])/margin[2]
  f01 <- sum(data$n011[nodes])/margin[3]
  f00 <- sum(data$n001[nodes])/margin[4]
  g11 <- 1-1/(max(f00/(1-f00),f10/(1-f10))+max(f00/(1-f00),f01/(1-f01))-f00/(1-f00)+1)
  
  n111 <- sum(data$n111[nodes])
  E111 <- g11*margin[1]
  E111 <- ifelse(E111==0,1e-6,E111)
  
  chisq <- (n111-E111-0.5)/sqrt(E111)
  chisq <- ifelse(is.infinite(chisq),NA,chisq)
  return(data.frame(NodeID=tree_str$tree[tree_ind], n111, E111, RR=n111/E111, chisq))
}

chisq_method <- function(data, tree_str){
  do.call(rbind,lapply(tree_str$tree, calc_chisq, data=data, tree_str=tree_str))
}

# CRR & CSS ---------------------------------------------------------------
# 72  1500_0694  156 3.7619445 3.093453e+02 1.8250998      0
calc_PRR <- function(x, nodes, drugs){
  
  row_g <- rep(0, nrow(x))
  row_g[nodes] <- 1
  row_g <- factor(row_g, levels=c(1,0))
  
  col_g <- rep(0, ncol(x))
  col_g[drugs] <- 1
  col_g <- factor(col_g, levels=c(1,0))
  
  n <- tapply(x, list(row_g[row(x)], col_g[col(x)]), sum)
  
  N11 <- n[1,1]
  N10 <- n[2,1]
  N01 <- n[1,2]
  N00 <- n[2,2]
  
  PRR <- (N11/(N11+N10))/(N01/(N01+N00))
  PRR025 <- exp(log(PRR)+qnorm(0.025)*sqrt(1/N11-1/(N11+N10)+1/N01-1/(N01+N00)))
  PRR975 <- exp(log(PRR)+qnorm(0.975)*sqrt(1/N11-1/(N11+N10)+1/N01-1/(N01+N00)))
  chisq <- (sum(n)*(abs(N11*N00-N10*N01)-sum(n)/2)^2)/((N11+N10)*(N11+N01)*(N01+N00)*(N10+N00))
  
  return(data.frame(PRR, PRR025, PRR975, chisq))
}

calc_CSS <- function(data, tree_str, tree){
  
  tree_ind <- which(tree_str$tree==tree)
  nodes <- tree_str$tree_ind[[tree_ind]]
  
  x <- as.matrix(data[,c('n111','n101','n011','n001')])
  
  PRR11 <- calc_PRR(x, nodes, drugs=1)
  PRR10 <- calc_PRR(x, nodes, drugs=c(1,2))
  PRR01 <- calc_PRR(x, nodes, drugs=c(1,3))
  
  CRR <- PRR11$PRR/max(PRR10$PRR,PRR01$PRR)
  CRR <- ifelse(is.na(CRR),0,CRR)
  CSS <- PRR11$PRR025/max(PRR10$PRR975,PRR01$PRR975)
  CSS <- ifelse(is.na(CSS),0,CSS)
  
  chisq <- PRR11$chisq
  chisq <- ifelse(is.na(CSS)&is.infinite(chisq),0,chisq)
  
  return(data.frame(NodeID=tree_str$tree[tree_ind], n111=sum(data$n111[nodes]), 
                    PRR=PRR11$PRR, chisq, CRR, PRR025=PRR11$PRR025, CSS))
}

CSS <- function(data, tree_str){
  do.call(rbind,lapply(tree_str$tree, calc_CSS, data=data, tree_str=tree_str))
}

# Proposed Method ---------------------------------------------------------

param_multinom_to_norm <- function(n){
  n <- ifelse(n==0, n+0.5, n)
  mean <- log(n[1])-log(n[2])-log(n[3])+log(n[4])
  sd <- sqrt(1/n[1]+1/n[2]+1/n[3]+1/n[4])
  return(c(mean, sd))
}

param_extract <- function(n_mat, margin_drug, node){
  param_G <- param_multinom_to_norm(n=as.vector(n_mat[node,]))
  param_R <- param_multinom_to_norm(n=as.vector(margin_drug-n_mat[node,]))
  
  mean <- param_G[1]-param_R[1]
  sd <- sqrt(param_G[2]^2+param_R[2]^2)
  return(c(mean, sd))
}

calc_LLR <- function(n_mat, params, node){
  
  param <- params[node,]
  
  if(param[1]<=0|n_mat[node,1]==0){
    LLR <- 0
  }else{
    LLR <- -dnorm(0,param[1],param[2],log=T)
  }
  
  return(LLR)
}

LLR_method <- function(data, tree_str){
  
  x <- data[,c('n111','n101','n011','n001')]
  
  margin_drug <- as.numeric(apply(x,2,sum))
  margin_AE <- as.numeric(apply(x,1,sum))
  n_mat <- as.matrix(x)
  
  for(i in (nrow(tree_str$tree_df)+1):length(tree_str$tree_ind)){
    if(length(tree_str$tree_ind[[i]])>1){
      n_mat <- rbind(n_mat, apply(n_mat[tree_str$tree_ind[[i]],],2,sum))
    }else{
      n_mat <- rbind(n_mat, n_mat[tree_str$tree_ind[[i]],])
    }
  }
    
  params <- t(sapply(1:nrow(n_mat), param_extract, n_mat=n_mat, margin_drug=margin_drug))
  
  LLR <- sapply(1:nrow(n_mat), calc_LLR, n_mat=n_mat, params=params)
  LLR_s <- LLR_MCMC(n_mat, margin_AE, margin_drug, iter=999)
  
  p_value <- sapply(LLR, calc_p, LLR_s)
  
  return(data.frame(NodeID=tree_str$tree, n111=n_mat[,1], LRMIS=params[,1], LLR, p_value))
}

MCMC_filter <- function(n_mat, margin_drug){
  p_mat <- sweep(n_mat,2,margin_drug,'/')
  index <- which((p_mat[,1]*p_mat[,4])/(p_mat[,2]*p_mat[,3])>1&n_mat[,1]>0)
  return(index)
}

LLR_MCMC <- function(n_mat, margin_AE, margin_drug, w, iter=999){
  
  LLR <- rep(NA,iter)
  n_mat_s <- n_mat
  
  for(i in 1:iter){
    
    n_mat_s<- t(sapply(margin_AE, rmultinom, n=1, prob=margin_drug))
    ind <- MCMC_filter(n_mat_s, margin_drug)
    
    params_s <- t(sapply(1:nrow(n_mat_s[ind,]), param_extract, n_mat=n_mat_s[ind,], margin_drug=margin_drug))
    LLR[i] <- max(sapply(1:nrow(n_mat_s[ind,]), calc_LLR, n_mat=n_mat_s[ind,], params=params_s), na.rm=T)
  }
  return(LLR)
}

calc_p <- function(x, y){
  rank(-c(x,y))[1]/length(c(x,y))
}

# Performance -------------------------------------------------------------

calc_perf <- function(detect, signal){
  n_detect <- sum(detect)
  sens <- sum(detect==1&signal==1)/sum(signal) 
  FDR <- sum(detect==1&signal==0)/sum(detect)
  return(data.frame(n_detect, sens, FDR))
}

summary_perf <- function(data){
  power <- mean(data$n_detect>0, na.rm=T)
  sens <- mean(data$sens[data$n_detect>0], na.rm=T)
  FDR <- mean(data$FDR[data$n_detect>0], na.rm=T)
  AUC <- mean(data$AUC[data$n_detect>0], na.rm=T)
  return(data.frame(power, sens, FDR, AUC))
}  

run_simulation <- function(iter, margin_drug, n_AE, H0, H1, n_signal, drug_bias,
                           num_AE_bias, AE_bias, tree_str, cut=0, type='additive'){
  
  margin_AE <- as.vector(rmultinom(1,sum(margin_drug),runif(n_AE)))
  df <- data.frame(NodeID=tree_str$tree[1:n_AE], n111=numeric(n_AE), n101=numeric(n_AE),
                   n011=numeric(n_AE), n001=numeric(n_AE))
   
  signal_ind <- sample(1:n_AE,n_signal)

  df[-signal_ind,2:5] <- t(sapply(margin_AE[-signal_ind], rmultinom, n=1, prob=margin_drug*H0))
  df[signal_ind,2:5] <- t(sapply(margin_AE[signal_ind], rmultinom, n=1, prob=margin_drug*H1))
  
  signal_vec <- tree_signal(signal=signal_ind, H0=H0, H1=H1, tree_str=tree_str, margin_AE=margin_AE, type=type)
  signal <- ifelse(signal_vec>cut, 1, 0)
  
  d_n111 <- rbinom(nrow(df), df$n111, p=drug_bias[1])
  d_n101 <- rbinom(nrow(df), df$n101, p=drug_bias[2])
  d_n011 <- rbinom(nrow(df), df$n011, p=drug_bias[3])
  d_n001 <- rbinom(nrow(df), df$n001, p=drug_bias[4])
  
  #if(AE_bias!=0 & num_AE_bias!=0){
  #  AE_bias_list <- sample(1:nrow(df),num_AE_bias)
  #  p <- runif(num_AE_bias,0,AE_bias)
  #  for(j in 2:5){
  #    df[AE_bias_list,j] <- df[AE_bias_list,j] - rbinom(num_AE_bias,df[AE_bias_list,j],p=p)
  #  }
  #}
  
  for(j in signal_ind){
    p <- AE_bias
    df[j,2:5] <- df[j,2:5] - rbinom(4,unlist(df[j,2:5]),p=p)
  }

  df$n111 <- df$n111 - d_n111
  df$n101 <- df$n101 - d_n101
  df$n011 <- df$n011 - d_n011
  df$n001 <- df$n001 - d_n001
  df[df<0] <- 0
  
  omega_res <- omega_shrinkage(data=df, tree_str)
  chisq_res <- chisq_method(data=df, tree_str)
  CRR_CSS_res <- CSS(data=df, tree_str)
  CRR_res <- CRR_CSS_res[,1:5]
  CSS_res <- CRR_CSS_res[,c(1:2,6:7)]
  LLR_res <- LLR_method(data=df, tree_str)
  
  omega_res$detect <- ifelse(omega_res$omega025>0,1,0)
  chisq_res$detect <- ifelse(chisq_res$chisq>2,1,0)
  CRR_res$detect <- ifelse(CRR_res$n111>=3&CRR_res$PRR>2&CRR_res$chisq>4&CRR_res$CRR>2,1,0)
  CSS_res$detect <- ifelse(CSS_res$PRR025>1&CSS_res$CSS>1,1,0)
  LLR_res$detect <- ifelse(LLR_res$p_value<0.05,1,0)
  
  omega_perf <- chisq_perf <- CRR_perf <- CSS_perf <- 
    LLR_perf <- data.frame(iter=iter, n_detect=NA, sens=NA, FDR=NA, AUC=NA)
  
  omega_perf[,-c(1,5)] <- calc_perf(omega_res$detect, signal)
  chisq_perf[,-c(1,5)] <- calc_perf(chisq_res$detect, signal)
  CRR_perf[,-c(1,5)] <- calc_perf(CRR_res$detect, signal)
  CSS_perf[,-c(1,5)] <- calc_perf(CSS_res$detect, signal)
  LLR_perf[,-c(1,5)] <- calc_perf(LLR_res$detect, signal)

  if(any(signal_vec>cut)){
    omega_perf$AUC <- roc(signal, omega_res$omega025, quiet=T)$auc
    chisq_perf$AUC <- roc(signal, chisq_res$chisq, quiet=T, na.rm=T)$auc
    CRR_perf$AUC <- roc(signal, CRR_res$CRR, quiet=T)$auc
    CSS_perf$AUC <- roc(signal, CSS_res$CSS, quiet=T)$auc
    LLR_perf$AUC <- roc(signal, LLR_res$LLR, quiet=T)$auc
  }
  
  return(list(omega=omega_perf, chisq=chisq_perf, CRR=CRR_perf, CSS=CSS_perf, LLR=LLR_perf))
}

# Utils -------------------------------------------------------------------

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}
