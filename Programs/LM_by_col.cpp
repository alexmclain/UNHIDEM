// [[Rcpp::depends(RcppArmadillo)]]

#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List LM_by_col(const arma::vec y, const arma::mat X) {
  
  int n = X.n_rows, d = X.n_cols;
  
  arma::vec X_mean2 = arma::trans(mean(X,0));
  
  arma::vec X_mean = mean(X,1);
  arma::vec Xp = X_mean; 
  double mean_Xp = mean(Xp);
  arma::colvec diff = Xp - mean_Xp;
  double SST = arma::as_scalar(arma::trans(diff)*diff);
  
  arma::mat coef_mat(d,3);
  arma::mat se(d,3);
  arma::mat sig(d,1);
  arma::mat R_sq(d,1);
  arma::mat X_var(d,1);
  arma::mat X1(n,1,arma::fill::ones);
  
  for(int col=0;col<d; ++col){
    
    arma::mat X2_Xp=arma::join_rows(X1,X.col(col));
    double X_mn = mean(X.col(col));
    arma::colvec X_diff = (X.col(col)-X_mn);
    double t_var = sum(arma::trans(X_diff)*X_diff);
    X_var(col) = t_var;
    
    arma::colvec coef_Xp = arma::solve(X2_Xp, Xp); 
    arma::colvec resid_Xp = Xp - X2_Xp*coef_Xp; 
    
    double SSE = arma::as_scalar(arma::trans(resid_Xp)*resid_Xp);
    double R2 = 1-SSE/SST;
    
    
    arma::mat X2=arma::join_horiz(X2_Xp,Xp);
    
    int n = X2.n_rows, k = X2.n_cols;
    
    arma::colvec coef = arma::solve(X2, y); 
    arma::colvec resid = y - X2*coef; 
    
    double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
    arma::colvec stderrest = 
      arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X2)*X2)) );
    double s = sqrt(sig2);
    
    coef_mat(col,0) = coef(0);
    coef_mat(col,1) = coef(1);
    coef_mat(col,2) = coef(2);
    se(col,0) = stderrest(0);
    se(col,1) = stderrest(1);
    se(col,2) = stderrest(2);
    sig(col,0) = s; 
    R_sq(col,0) = R2; 
  }
  
  arma::vec mul = mean(sig,0);
  arma::colvec t_val = coef_mat.col(1)/se.col(1);
  
  return List::create(Named("Coefficients") = coef_mat,Named("StdErr") = se, Named("mul") = mul,
                      Named("R_sq")=R_sq ,Named("sig")=sig ,Named("X_var") = X_var,Named("t_val")=t_val,Named("X_mean2")=X_mean2);
}








