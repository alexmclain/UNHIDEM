



E_step_func <- function(beta_t,beta_var,df,adj,alpha,one.sided=FALSE){
  
  ### Get and clean up T statistics
  T_vals <- beta_t/sqrt(beta_var)
  T_vals[beta_var<=0] <- 0
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0   
  M <- length(T_vals) 
  if(!one.sided){
    p_vals <- pt(abs(T_vals),df=df,lower.tail = FALSE)*2
  }else{
    p_vals <- pt(T_vals,df=N-3,lower.tail = FALSE)
  }
  
  ### Estimate the probability of a null
  p_hat  <- pi0est(p_vals,lambda = alpha)$pi0-1/M
  ### Estimate the lfdr
  delta <- 1-lfdr_T(T=T_vals,pi0=p_hat,trunc=TRUE,monotone = TRUE,adj = adj,df_val = N-3,one.sided = one.sided)$lfdr
  if(one.sided){
    tilde_pos <- norm_trun_func(beta_t,beta_var,side="Pos")
    beta_t <- tilde_pos$beta_tilde
    beta_var <- tilde_pos$beta_tilde_var
  }
  
  ret <- list(beta_tilde=beta_t,beta_tilde_var=beta_var,delta=delta,lfdr=1-delta,pi0 = p_hat,p_vals=p_vals,T_vals=T_vals,df=df)
  return(ret)
}




W_update_fun <- function(Z,Z_2,beta_tilde,beta_tilde_var,delta){
  Z_delta <- MVM(Z,delta*beta_tilde)$Res 
  W_ast<- c(Row_sum(as.matrix(Z_delta))$Rowsum)
  W_ast[is.nan(W_ast) | is.na(W_ast)] <- 0
  W_ast_var <- NULL
  if(!is.null(Z_2)){
    Z_delta2 <- MVM(Z_2,delta*(beta_tilde_var+beta_tilde^2-delta*beta_tilde^2))$Res 
    W_ast_var <- c(Row_sum(as.matrix(Z_delta2))$Rowsum)
    W_ast_var[is.nan(W_ast_var) | is.na(W_ast_var)] <- 0
  }
  return(list(W_ast=W_ast,W_ast_var=W_ast_var))
}




pi0est <- function (p, lambda = seq(0.05, 0.95, 0.05), pi0.method = c("smoother", 
                                                                      "bootstrap"), smooth.df = 3, smooth.log.pi0 = FALSE, ...) 
{
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda)
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  }
  else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", 
                                       ll, ".", sep = ""), "If length of lambda greater than 1,", 
                       "you need at least 4 values.")))
  }
  else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  }
  else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec = lambda))[ind])/(length(p) * 
                                                                   (1 - lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      }
      else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W/(m^2 * (1 - lambda)^2)) * (1 - W/m) + (pi0 - 
                                                         minpi0)^2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    }
    else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda, 
              pi0.smooth = pi0Smooth))
}

lfdr <- function (p, pi0 = NULL, trunc = TRUE, monotone = TRUE, transf = c("probit","logit"), adj = 1.5, eps = 10^-8, ...) 
{
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  }
  else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    
    y[y<=0] <- 1e-10
    lfdr <- pi0 * dnorm(x)/y
  }
  else {
    x <- log((p + eps)/(1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x)/(1 + exp(x))^2
    lfdr <- (pi0 * dx)/y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  f <- list(x=sort(x),y=lfdr[order(x)])
  res <- list(lfdr=lfdr_out,f=f)
  return(res)
}


lfdr_T <- function (T, pi0 = NULL, trunc = TRUE, monotone = TRUE, adj = 1.5,df_val=NULL,one.sided=FALSE, ...) 
{
  lfdr_out <- T
  rm_na <- !is.na(T)
  T <- T[rm_na]
  if (is.null(pi0)) {
    stop("pi0 must be given.")
  }
  n <- length(T)
  myd <- density(T, adjust = adj)
  mys <- smooth.spline(x = myd$x, y = myd$y)
  y <- predict(mys, T)$y
  y[y<=0] <- 1e-10
  lfdr <- pi0 * dt(T,df=df_val)/y
  
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone & !one.sided) {
    T_neg <- T[T<0]
    lfdr_neg <- lfdr[T<0]
    T_pos <- T[T>=0]
    lfdr_pos <- lfdr[T>=0]
    
    o <- order(T_neg, decreasing = FALSE)
    ro <- order(o)
    lfdr[T<0] <- cummax(lfdr_neg[o])[ro]
    
    o <- order(T_pos, decreasing = FALSE)
    ro <- order(o)
    lfdr[T>=0] <- cummin(lfdr_pos[o])[ro]
  }
  if (monotone & one.sided) {
    o <- order(T, decreasing = TRUE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  f <- list(x=sort(T),y=lfdr[order(T)])
  res <- list(lfdr=lfdr_out,f=f)
  return(res)
}
