


source("Programs/UNHIDEM_funcs.R")

Ex_dat <- read.csv("Data/Example_data.csv")
Ex_dat_tru <- read.csv("Data/Example_data_truth.csv")


Y <- Ex_dat[,1]
Z <- as.matrix(Ex_dat[,-1])
M <- dim(Z)[2]


##Truth
beta_tr <- Ex_dat_tru$x  #True beta coefficients
signal <- c(1:M)[beta_tr!=0] #Indicies of non-null beta coefficients
eta_i <- apply(t(Z)*beta_tr,2,sum) #The true signal

## Set convergence criteria
maxit <- 1e5
ep <- 1e-8
alpha <- 0.1
B=5
adj=3
one.sided=FALSE
verbose=TRUE
plot_ind=TRUE

##### First run an example imputting the true signal and indicies on non-null beta coefficients
unhidem <- UNHIDEM1.46(Y,Z,alpha,ep = ep, B=B, adj=adj,maxit=maxit,one.sided = one.sided,verbose=verbose,plot_ind = plot_ind,eta_i = eta_i,signal = signal)

full_res <- unhidem$full_res
pred_res <- prediction.unhidem(full_res,Y,Z,alpha)
head(pred_res)
MT_res <- MTR(full_res$E_step,alpha,signal)
MT_res$BH_sum



##### Second run an analysis using test data

Ex_dat_test <- read.csv("Data/Example_test_data.csv")
Y_test <- Ex_dat_test[,1]
Z_test <- as.matrix(Ex_dat_test[,-1])

unhidem <- UNHIDEM1.46(Y,Z,alpha,ep = ep, B=B, adj=adj,maxit=maxit,one.sided = one.sided,Y_test=Y_test,Z_test=Z_test,verbose=verbose,plot_ind = plot_ind)

pred_res_test <- prediction.unhidem(full_res,Y_test,Z_test,alpha)
head(pred_res_test)

## Proportion of test PIs that contain the test observation
mean(1*I(Y_test>pred_res_test$PI[,1] & Y_test<pred_res_test$PI[,2]))

## Proportion of test CIs that contain the test true signal
eta_test <- apply(t(Z_test)*beta_tr,2,sum) #The true signal for test data
mean(1*I(eta_test>pred_res_test$CI[,1] & eta_test<pred_res_test$CI[,2]))







