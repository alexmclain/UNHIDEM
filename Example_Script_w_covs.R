
source("Programs/UNHIDEM_funcs.R")

Ex_dat <- read.csv("Data/Example_data_w_covs.csv")
Ex_dat_tru <- read.csv("Data/Example_data_truth_w_covs.csv")


Y <- Ex_dat[,1]
X <- Ex_dat[,3:4]
Z <- as.matrix(Ex_dat[,-c(1:4)])
M <- dim(Z)[2]


##Truth
beta_tr <- Ex_dat_tru$x[1:M]  #True beta coefficients
beta_X_tr <- Ex_dat_tru$x[(M+1):(M+3)] #True coefficients for adjustment variables and intercept
signal <- c(1:M)[beta_tr!=0] #Indicies of non-null beta coefficients
eta_i <- apply(t(Z)*beta_tr,2,sum) #The true signal

### Set modeling parameters
ep <- 1e-8 #Convergence criteria
maxit <- 1e5 #Max iterations
alpha <- 0.1 #Type I error
B=5 # Number of groups for correlation estimation
adj=3 # Used for estimation of p-value distribution

one.sided=FALSE
verbose=TRUE
plot_ind=TRUE

##### First run an example imputting the true signal and indicies on non-null beta coefficients
unhidem <- UNHIDEM(Y,Z,X,alpha,ep = ep, B=B, adj=adj,maxit=maxit,one.sided = one.sided,verbose=verbose,plot_ind = plot_ind,eta_i = eta_i,signal = signal)

full_res <- unhidem$full_res
pred_res <- prediction.unhidem(full_res,Y,Z,X,alpha)
head(pred_res)

MT_res <- MTR(full_res$E_step,alpha,signal)
MT_res$BH_sum

