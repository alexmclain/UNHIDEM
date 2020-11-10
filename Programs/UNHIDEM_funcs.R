require(RcppArmadillo)
require(Rcpp)

source("Programs/E_step_functions.R")
sourceCpp("Programs/UNHIDEM_cpp0_5.cpp")
sourceCpp("Programs/UNHIDEM_cpp0_5_w_covs.cpp")
sourceCpp("Programs/LM_by_col.cpp")
sourceCpp("Programs/LM_w_COVS_by_col.cpp")
sourceCpp("Programs/Row_sum.cpp")
sourceCpp("Programs/Col_sum.cpp")
sourceCpp("Programs/MVM.cpp")



UNHIDEM <- function(Y,Z,X=NULL,alpha=0.05,ep=1e-8,B=5,adj=3,maxit=1000,one.sided=TRUE,Y_test=NULL,Z_test=NULL,verbose=FALSE,signal=NULL,eta_i=NULL,plot_ind=FALSE){
  
  M <- dim(Z)[2]
  N <- dim(Z)[1]
  if(!is.null(X)){X <- as.matrix(X)}
  full_res <- UNHIDEM0.96(Y=Y,Z=Z,X=X,alpha,verbose=verbose,signal=signal,maxit=maxit,eta_i=eta_i,one.sided=one.sided,ep=ep,B=B,adj=adj,plot_ind=plot_ind,Y_test=Y_test,Z_test=Z_test)
  res <- list(full_res=full_res, total_iter=full_res$count)
  return(res)
}


UNHIDEM0.96 <- function(Y,Z,X=NULL,alpha,verbose=TRUE,signal,maxit=1000,eta_i=NULL,one.sided=TRUE,ep=1e-8,B=5,adj=NULL,plot_ind=FALSE,Y_test=NULL,Z_test=NULL,X_test=NULL){
  
  if(is.null(adj)){adj <- 3/(1+1*I(one.sided))}
  
  ##### Setting initial values and initializing outputs ####
  M <- dim(Z)[2]
  N <- dim(Z)[1]
  if(!is.null(X)){
    p <- dim(X)[2]
    beta_X <- rep(0,p)
    }
  delta <- beta_t <- beta_var <- beta_tilde <- beta_tilde_var <- rep(0,M)
  W_ast    <- rep(0,N)
  W_ast_var<- W_ast+1
  count <- 0
  conv_check <- 0
  cor_vec <- rep(0,M) 
  plot_dat <- NULL
  Z_2 <- Z*Z
  Xt_conv1 <- Xt_conv2 <- 0
  signal_track <- report_pred <- NULL
  while(count < maxit & conv_check < 1){
    
    W_ast_old <- W_ast
    W_ast2_old <- W_ast_var
    count <- count+1
    fact <- 1/(count+1)
    
    ### Performing the M-step.
    if(count==1){LR_update <- LR_cpp(Y,Z,X)}
    if(count>1){LR_update <- theta_update_0_5.cpp(Y,Z,X,W_ast,W_ast_var,delta,beta_tilde,beta_tilde_var,Z_2)}
    
    beta_t_new <- LR_update$coef[,2]
    beta_var_new <- LR_update$obs_SE[,2]^2
    obs_t <- LR_update$coef[,2]/LR_update$obs_SE[,2]
    
    ### Performing the Mixing step.
    if(count>1){cor_vec <- corr_func(beta_t_new,beta_t,obs_t,B=B,one.sided = one.sided)}
    beta_t <- beta_t_new*fact + beta_t*(1-fact)
    beta_var <- beta_var_new*fact^2 + beta_var*(1-fact)^2 + 2*fact*(1-fact)*cor_vec*sqrt(beta_var_new)*sqrt(beta_var)
    
    ### Performing the E-step.
    E_step <- E_step_func(beta_t,beta_var,df=N-3-p,adj=adj,alpha=alpha,one.sided=one.sided)
    beta_tilde <- E_step$beta_tilde
    beta_tilde_var <- E_step$beta_tilde_var
    delta <- E_step$delta
    
    ### Updating W and W2
    W_W2_update <- W_update_fun(Z,Z_2,beta_tilde,beta_tilde_var,delta)
    W_ast <- W_W2_update$W_ast
    W_ast_var <- W_W2_update$W_ast_var

    if(var(W_ast)>0){
      if(count>1){ ### Check Convergence
        Xt_conv1 <- cor(W_ast_old,W_ast) + cor(W_ast2_old+W_ast_old^2,W_ast_var+W_ast^2)
        if((Xt_conv1) > (2-ep)){conv_check <- conv_check+1}
      }
      ### Run calibration model
      mod <- W_regression(Y,W_ast,W_ast_var + W_ast^2,X)
      Y_pred <- mod$Y_pred  
    }
    if(var(W_ast)==0){
      cat("Warning loop completely recycled back to beta=0.\n")
      break
    }
    
    ### Getting prediction error if eta_i or test data given.
    if(!is.null(eta_i)){report_pred = mean((eta_i-Y_pred)^2)}
    if(!is.null(Y_test) & !is.null(Z_test)){
      W_W2_test <- W_update_fun(Z_test,Z_2=NULL,E_step$beta_tilde,E_step$beta_tilde_var,E_step$delta)
      
      Y_pred_test <- cbind(1,W_W2_test$W_ast,X_test)%*%mod$coef
      report_pred = mean((Y_test-Y_pred_test)^2)
    }

    ### Performing hypothesis testing on current estimates.
    MTR_res <- MTR(E_step,alpha,signal)
    if(!is.null(signal)){signal_track <- c(MTR_res$BH_sum$LFDR_sum[3],M-MTR_res$BH_sum$LFDR_sum[1],MTR_res$BH_sum$LFDR_sum[2],length(signal)-MTR_res$BH_sum$LFDR_sum[2])}

    ### Outputting results if verbose=TRUE
      if((count %in% c(seq(100,maxit,100)) | conv_check==1) & verbose){report_func(count,E_step,MTR_res,Xt_conv1,signal_track,report_pred)}    
    
    ### Storing iteration data
    plot_dat <- rbind(plot_dat,c(count,conv_check,signal_track,sum(MTR_res$BH_res$LFDR),sum(delta),cor(W_ast,Y),report_pred,mean(cor_vec),-log10(2-(Xt_conv1))))
  }

  
  
  
  ### Formatting iteration data for outputting
  if(!is.null(signal_track)){colnames(plot_dat) <- c("Iter","Conv_check","FP","TN","TP","FN","Total_Disc","Sum_delta","Corr_W_Y","Pred_err","Avg_corr","Conv")}else{
    if(!is.null(report_pred)){
      colnames(plot_dat) <- c("Iter","Conv_check","Total_Disc","Sum_delta","Corr_W_Y","Pred_err","Avg_corr","Conv")
    }else{
      colnames(plot_dat) <- c("Iter","Conv_check","Total_Disc","Sum_delta","Corr_W_Y","Avg_corr","Conv")
    }
  }
  plot_dat <- data.frame(plot_dat)
  
  ### Plotting iteration results if plot_ind=TRUE and either eta_i or test data is given
  if(plot_ind){
    if(is.null(report_pred)){cat("Warning: cannot plot without eta_i or test data.\n")}else{
      plot_function_res(plot_dat,!is.null(Y_test))
    }
  }
  
  if(conv_check==0){cat("Warning: convergence criteria not met. Set different convergence criteria or raise maximum number of iterations.\n")}
  
  return(list(E_step=E_step,Calb_mod=mod,count=count,plot_dat=plot_dat))
}


MTR <- function(E_step,alpha,signal=NULL){
  
  p_vals <- E_step$p_vals
  lfdr_val <- E_step$lfdr
  p_hat <- E_step$pi0
  M <- length(p_vals)
  alpha_hat <- alpha/p_hat
  
  T_R <- p.adjust(p_vals,method = "BY")
  R_BH <- p.adjust(p_vals,method = "BH")
  
  threshold <- 0
  lfdr_val[is.na(lfdr_val)] <- 1
  if(min(lfdr_val)<alpha){threshold <- max(sort(lfdr_val)[cumsum(sort(lfdr_val))<alpha])}
  
  BH_res <- data.frame(BY = 1*I(T_R<alpha),BH = 1*I(R_BH<alpha),LFDR = 1*I(lfdr_val <= threshold))
  
  R_BY <- p.adjust(p_vals,method = "bonferroni")
  R2_BY <- p.adjust(p_vals,method = "holm")
  
  Bonf_res <- data.frame(Holm = 1*I(R2_BY<alpha),Bonf = 1*I(R_BY<alpha))
  
  BH_sum <- NULL
  Bonf_sum <- NULL
  if(!is.null(signal)){
    n_signal <- !(1:M %in% signal)
    BH_sum   <- c(sum(BH_res$BH),sum(BH_res$BH[signal]),sum(BH_res$BH[n_signal]))
    LFDR_res <- c(sum(BH_res$LFDR),sum(BH_res$LFDR[signal]),sum(BH_res$LFDR[n_signal]))
    BY_res  <- c(sum(BH_res$BY),sum(BH_res$BY[signal]),sum(BH_res$BY[n_signal]))
    BH_sum <- data.frame(BY_sum = BY_res, LFDR_sum = LFDR_res,BH_sum = BH_sum)
    
    Bonf_sum <- c(sum(Bonf_res$Bonf),sum(Bonf_res$Bonf[signal]),sum(Bonf_res$Bonf[n_signal]))
    Holm_res2 <- c(sum(Bonf_res$Holm),sum(Bonf_res$Holm[signal]),sum(Bonf_res$Holm[n_signal]))
    
    Bonf_sum <- data.frame(Holm_sum = Holm_res2,Bonf_sum=Bonf_sum)
  }
  
  
  return(list(BH_res = BH_res,Bonf_res = Bonf_res,Bonf_sum=Bonf_sum,BH_sum = BH_sum))
}


LR_cpp <- function(Y_dat,BR_dat,X=NULL){
  
  N <- length(Y_dat)
  if(!is.null(X)){
  p <- ncol(X)
  LRcpp <- LM_w_COVS_by_col(Y_dat,BR_dat,as.matrix(X))
  }else{
    p <- 0
    LRcpp <- LM_by_col(Y_dat,BR_dat)
  }
  t_val <- LRcpp$Coefficients[,2]/LRcpp$StdErr[,2]
  p_val <- pt(t_val,df=N-3-p,lower.tail = FALSE)
  
  ret <- list(coef = LRcpp$Coefficients,obs_SE = LRcpp$StdErr,T_val = t_val, p_val = p_val,sig_vec = LRcpp$sig)
  
  return(ret)
}

theta_update_0_5.cpp <- function(Y_dat,BR_dat,X,W_ast,W_ast_var,delta,beta_vec,beta_var,BR_dat2){
  # This funcion updates all of the theta estimates the fast way.
  
  N <- length(Y_dat)
  if(!is.null(X)){
    p <- ncol(X)
    LRcpp <- UNHIDEM_cpp0_5_w_covs(Y_dat,BR_dat,W_ast,W_ast_var,delta,beta_vec,beta_var,BR_dat2,as.matrix(X))
  }else{
    p <- 0
    LRcpp <- UNHIDEM_cpp0_5(Y_dat,BR_dat,W_ast,W_ast_var,delta,beta_vec,beta_var,BR_dat2)
  }
  t_val <- LRcpp$T_statistics
  p_val <- pt(t_val[,2],df=N-3-p,lower.tail = FALSE)
  
  ret <- list(coef = LRcpp$Coefficients,obs_SE = LRcpp$StdErr,T_val = t_val, p_val = p_val,log_like=LRcpp$Log_like,sig_vec = LRcpp$Sigma)
  
  return(ret)
}


corr_func <- function(beta_t_new,beta_t,T_vals,B=5,one.sided=FALSE){
  
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0
  if(!one.sided){T_vals <- abs(T_vals)}
  xx <- data.frame(a = beta_t_new , b = beta_t,T_vals=T_vals)
  brks <- c(min(T_vals),max(T_vals))
  if(B>1){
    min_prob <- 1/2*I(one.sided)+1/B*I(!one.sided)
    jumps <- 1/(B*(1+1*I(one.sided)))
    brks <- with(xx, quantile(T_vals, probs = c(0,seq(min_prob,1,jumps)),na.rm = TRUE))
  }
  xx <- within(xx, group <- cut(T_vals, breaks = brks, labels = 1:(length(brks)-1),include.lowest = TRUE))
  result <- by(xx[,1:2], xx$group, function(x) {cor(x$a, x$b)})
  result.dataframe <- as.data.frame(as.matrix(result))
  cor_vec <- result.dataframe$V1[as.numeric(xx$group)]
  return(cor_vec)
}

norm_trun_func <- function(mu,var,side="Pos")
{
  mu <- mu*I(side=="Pos") - mu*I(side=="Neg")
  T_vals <- mu/sqrt(var)
  log_norm_ratio <- dnorm(-T_vals,log=TRUE)-pnorm(-T_vals,log.p = TRUE,lower.tail = FALSE)
  norm_ratio <- exp(log_norm_ratio)
  
  beta_tilde <- (mu+sqrt(var)*norm_ratio)
  beta_tilde <- beta_tilde*I(side=="Pos") - beta_tilde*I(side=="Neg")
  beta_tilde_var <- var*(1-T_vals*norm_ratio-(norm_ratio)^2)
  return(list(beta_tilde=beta_tilde,beta_tilde_var=beta_tilde_var))
}


W_regression <- function(Y_dat,W,W2,X){
  
  Wmat <- cbind(1,W,X)
  WWpr <- t(Wmat)%*%Wmat
  WWpr[2,2] <- sum(W2)
  N <- length(Y_dat)
  if(det(WWpr)!=0){
    WWpr_inv <- solve(WWpr)
    beta_w <- WWpr_inv%*%t(Wmat)%*%Y_dat
    Y_pred <- Wmat%*%beta_w
    hat <- diag(Wmat%*%WWpr_inv%*%t(Wmat))
    resid <- (Y_dat - Y_pred)
    RSS <- sum(resid^2)
    SST <- sum((Y_dat - mean(Y_dat))^2)
    R_squared <- 1 - RSS/SST
    VCV = RSS/N*t(WWpr_inv)%*%(t(Wmat)%*%Wmat)%*%WWpr_inv
  }
  if(det(WWpr)==0){
    final_mod <- lm(Y_dat~0+Wmat)
    beta_w <- as.numeric(final_mod$coefficients)
    beta_w[is.na(beta_w)] <- 0
    Y_pred <- as.numeric(final_mod$fitted.values)
    hat <- as.numeric(influence(final_mod)$hat)
    resid <- as.numeric(final_mod$residuals)
    RSS <- sum(resid^2)
    SST <- sum((Y_dat - mean(Y_dat))^2)
    R_squared <- 1 - RSS/SST
    VCV <- vcov(final_mod)
  }
  
  ret <- list(coef = beta_w,RSS = RSS,R_squared=R_squared,Y_pred=Y_pred,hat=hat,resid=resid,VCV=VCV)
  
  return(ret)
}

plot_function_res <- function(plot_dat,Y_test_ind){
  a = quantile(plot_dat$Pred_err,probs = 0.9)
  b = min(plot_dat$Pred_err)
  par(mar=c(4.2,4.5 , 0.3, 3.1) + 0.1,mfrow=c(1,1))
  ylab_val=expression(paste("Signal  ",MSE[t]))
  if(Y_test_ind){ylab_val=expression(paste("Test  ",MSE[t]))}
  plot(plot_dat$Iter,plot_dat$Pred_err,type = "l",ylim=c(b,a),xlim=range(plot_dat$Iter),xlab="Iteration",ylab=ylab_val,lwd=2,cex.lab=1.4,cex.axis=1.2,las=1)
  p_vec <- plot_dat$Total_Disc
  b9 <- min(p_vec)
  a9 <- quantile(p_vec,probs = 0.9)
  axis(4,las=1, at=seq(b,a,length.out = 4),labels = round(seq(b9,a9,length.out = 4),0), cex.axis=1.2,las=2)
  mtext(expression(S[t]), side=4, line=2.2, cex=1.2)
  trans_crit <- (p_vec-b9)/quantile(p_vec-b9,probs = 0.9)*(a-b)+b
  col_vec <- rep(1,length(p_vec))
  if(!Y_test_ind){
    denom <- plot_dat$Total_Disc
    denom[denom==0] <- 1
    col_vec <- ifelse(plot_dat$FP/denom<alpha,"grey60",1)
  }
  points(plot_dat$Iter,trans_crit,col=col_vec,pch=19,cex=0.5)
}

report_func <- function(count,E_step,MTR_res,Xt_conv1,signal_track,report_pred){
  CC_round <- round(-log10(2-(Xt_conv1)),2)
  if(is.null(signal_track)){
    disc <- sum(MTR_res$BH_res$LFDR)
    if(!is.null(report_pred)){
      if(report_pred<=1){report_pred <- signif(report_pred,3)}
      if(report_pred>1){report_pred <- round(report_pred,1)}
      cat("Iteration=",count,"Number of discoveries (using lfdr)=",disc,"Sum(delta)=",round(sum(E_step$delta),1)," MSE(test)=",report_pred,"Convergence Crit=",CC_round,"\n")
    }else{
      cat("Iteration=",count,"Number of discoveries (using lfdr)=",disc,"Sum(delta)=",round(sum(E_step$delta),1),"Convergence Crit=",CC_round,"\n")}
  } else {
    if(report_pred<=1){report_pred <- signif(report_pred,3)}
    if(report_pred>1){report_pred <- round(report_pred,1)}
    cat("Iteration=",count," Hyp testing (using lfdr) TP=",signal_track[3]," FP=",signal_track[1]," TN=",signal_track[2]," FN=",signal_track[4]," MSE(signal)=",report_pred," Convergence Crit=",CC_round,"\n",sep = "")
  }
}

prediction.unhidem <- function(full_res,Y,Z,X=NULL,alpha=NULL,Z_2=NULL){
  E_step <- full_res$E_step
  mod <- full_res$Calb_mod
  if(is.null(Z_2)){Z_2 <- Z*Z}
  
  W_W2_update <- W_update_fun(Z,Z_2,E_step$beta_tilde,E_step$beta_tilde_var,E_step$delta)
  W_ast <- W_W2_update$W_ast
  W_ast_var <- W_W2_update$W_ast_var
  
  if(!is.null(X)){X <- as.matrix(X)}
  mod_mat <- cbind(1,W_ast,X)
  pred_mean <- mod_mat%*%mod$coef
  Var_train <- diag(mod_mat%*%mod$VCV%*%t(mod_mat)) + (mod$VCV[2,2] + mod$coef[2]^2)*W_ast_var
  
  pred_res <- data.frame(Pred=pred_mean,Var=Var_train)
  CI_train <- PI_train <- NULL
  if(!is.null(alpha)){
    pred_res$CI <- cbind(pred_mean-qt(1-alpha/2,E_step$df)*sqrt(Var_train),pred_mean+qt(1-alpha/2,E_step$df)*sqrt(Var_train))
    pred_res$PI <- cbind(pred_mean-qt(1-alpha/2,E_step$df)*sqrt(Var_train+mod$RSS/E_step$df),pred_mean+qt(1-alpha/2,E_step$df)*sqrt(Var_train+mod$RSS/E_step$df))
  }
  
  return(pred_res)
}

