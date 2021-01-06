#############################################################################################################################
####### In this script, we only consider single change-point in the middle. Using the weakly sparse as an alternative #######
#############################################################################################################################
### sourcing functions
library("ggplot2")
library("glmnet")
library("matrixcalc")
library("fields")
library("sparsevar")
library("vars")
library("MTS")
library("mvtnorm")
library("xtable")
library("pracma")
library("lattice")
source("functions_SBDetection_Peiliang.R")

### simulation settings
reps <- 50
k <- 20
T <- 300
brk <- c(floor(T/2), T+1)
m <- length(brk)
p.t <- 1

gammas <- rep(0.25, 2)
ranks <- c(1, 3)
L_true <- vector("list", 2)
S_true <- vector("list", 2)
signals <- c(0.7, 0.8)

### generating sparse components
for(i in 1:2){
    L_true[[i]] = S_true[[i]] <- matrix(0, k, k)
}
for(i in 1:2){
    for(j in 1:(k-1)){
        S_true[[i]][j,j+1] <- (-1)^(i) * signals[i]
    }
}


### generating low-rank components by using the spectral decomposition
### instead randomly create the low-rank components

set.seed(100)
L_basis <- randortho(k)
singular_vals <- c(0.5, 0.1, 0.05)
for(i in 1:2){
    for(j in 1:ranks[i]){
        L_true[[i]] = L_true[[i]] + singular_vals[j] * (L_basis[,j] %*% t(L_basis[,j]))
    }
}

### check the rank and print the matrices of low-rank components
for(i in 1:2){
    print(qr(L_true[[i]])$rank)
    L_true[[i]] <- (signals[i] * gammas[i] / max(L_true[[i]])) * L_true[[i]]
}
L.full <- cbind(L_true[[1]], L_true[[2]])
print(plot.matrix(abs(L.full), 2))
print(norm(L_true[[1]] - L_true[[2]], "2"))

### force the stationary time series and make information ratio be embedded into the low-rank components
phi_full <- NULL
phi_true <- vector("list", 2)
max_eigen <- rep(0, 2)
Ms <- 0
for(i in 1:2){
    phi_true[[i]] <- S_true[[i]] + L_true[[i]]
    max_eigen[i] <- max(abs(eigen(phi_true[[i]])$values))
    phi_true[[i]] <- phi_true[[i]] * 0.9 / max_eigen[i]
    S_true[[i]] <- S_true[[i]] * 0.9 / max_eigen[i]
    L_true[[i]] <- L_true[[i]] * 0.9 / max_eigen[i]
    phi_full <- cbind(phi_full, phi_true[[i]])
    Ms <- max(Ms, max(S_true[[i]]))
}
print(plot.matrix(abs(phi_full), 2))
print(norm(phi_true[[1]] - phi_true[[2]], "2")) ### present the differences of two transaction matrices
print(paste("low rank components jump:", norm(L_true[[2]] - L_true[[1]], "2")))
print(paste("sparse components jump:", norm(S_true[[2]] - S_true[[1]], "2")))

############ Grid Search and cross-validation for estimation ##########
### let's apply the following method, first we select the tuning parameter for low-rank
### then we fix the low-rank tuning and then we select the tuning parameter for sparse
source("fista_LpS.R")

###### Grid search function ######
ts.separate <- function(data, npart = 5){
    n <- dim(data)[1]
    if(n < 10){
        break
        print("Warning: Too few data points!")
    }
    else{
        ntest <- floor(n / npart)
        ntrain <- n - ntest
        data_train <- data[1:ntrain,]
        data_test <- data[-(1:ntrain),]
    }
    return(list(train = data_train, test = data_test))
}

tuning.selection <- function(data, lambda.seq, mu.seq){
    ### separate data time series into training and testing set
    separate <- ts.separate(data)
    data_train <- separate$train
    data_test <- separate$test
    
    ntrain <- dim(data_train)[1]
    ntest <- dim(data_test)[1]
    
    X_train <- data_train[1:(ntrain-1),]
    Y_train <- data_train[2:ntrain,]
    
    X_test <- data_test[1:(ntest-1),]
    Y_test <- data_test[2:ntest,]
    
    ### grid search part: find the tuning parameter that minimizes the prediction error on test set
    n <- dim(data)[1]
    k <- dim(data)[2]
    grid <- matrix(0, nrow = length(lambda.seq), ncol = length(mu.seq))
    for(r in 1:length(lambda.seq)){
        for(c in 1:length(mu.seq)){
            lambda <- lambda.seq[r]
            mu <- mu.seq[c]
            fit <- fista.LpS(X_train, Y_train, lambda, mu, niter = 100, backtracking = TRUE, diag(k))
            if(qr(fit$lr.comp)$rank == 0 || qr(fit$lr.comp)$rank > k/2){
                grid[r,c] <- Inf
            }else{
                residual <- Y_test - X_test %*% (fit$sparse.comp + fit$lr.comp)
                grid[r,c] <- (1/ntest) * norm(residual, "F")^2
            }
        }
    }
    idx <- which(grid == min(grid), arr.ind = TRUE)
    final.lambda <- lambda.seq[idx[1]]
    final.mu <- mu.seq[idx[2]]
    return(list(grid = grid, lambda = final.lambda, mu = final.mu))
}


###### Change point detection functions ######
detect.LpS <- function(data, lambda1.seq, lambda2.seq, mu1.seq, mu2.seq, skip = 50){
    sse <- c()
    n <- dim(data)[1]
    lambda1.selected = lambda2.selected <- rep(0, n-2*skip)
    mu1.selected = mu2.selected <- rep(0, n-2*skip)
    
    lr_est <- vector("list", n-2*skip)
    sp_est <- vector('list', n-2*skip)
    for(t in (skip+1):(n-skip)){
        ###### segmentation of dataset ######
        seg1 <- data[1:t,]
        seg2 <- data[((t+1):n),]
        
        ###### tuning parameter selection for the left-side segment ######
        ret <- tuning.selection(seg1, lambda1.seq, mu1.seq)
        lambda1 = ret$lambda
        mu1 = ret$mu
        
        lambda1.selected[t-skip] <- lambda1
        mu1.selected[t-skip] <- mu1
        
        ### fitting the left-side dataset
        X1 <- seg1[1:(t-1),]
        Y1 <- seg1[2:t,]
        fit1 <- fista.LpS(X1, Y1, lambda1, mu1, 100, TRUE, diag(k))
        
        ###### tuning parameter selection for the right-side segment ######
        ret <- tuning.selection(seg2, lambda2.seq, mu2.seq)
        lambda2 = ret$lambda
        mu2 = ret$mu
        
        lambda2.selected[t-skip] <- lambda2
        mu2.selected[t-skip] <- mu2
        
        # fitting the right-side dataset
        X2 <- seg2[1:(n-t-1),]
        Y2 <- seg2[2:(n-t),]
        fit2 <- fista.LpS(X2, Y2, lambda2, mu2, 100, TRUE, diag(k))
        
        ###### find the residual and full SSE ######
        phi.hat1 <- t(fit1$lr.comp + fit1$sparse.comp)
        phi.hat2 <- t(fit2$lr.comp + fit2$sparse.comp)
        
        lr_est[[t-skip]] <- cbind(t(fit1$lr.comp), t(fit2$lr.comp))
        sp_est[[t-skip]] <- cbind(t(fit1$sparse.comp), t(fit2$sparse.comp))
        
        # pred.err1 <- (1/t) * norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        # pred.err2 <- (1/(n-t)) * norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        pred.err1 <- (1/n) * norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        pred.err2 <- (1/n) * norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        sse <- c(sse, pred.err1 + pred.err2)
        #print(paste("Finished time point:", t, "Lambdas:", lambda1, lambda2, "Mus:", mu1, mu2))
    }
    
    idx <- which.min(sse) + skip
    
    L_hat1 <- lr_est[[idx-skip]][,c(1:k)]; L_hat2 <- lr_est[[idx-skip]][,-(c(1:k))]
    S_hat1 <- sp_est[[idx-skip]][,c(1:k)]; S_hat2 <- sp_est[[idx-skip]][,-(c(1:k))]
    return(list(cp = idx, S_hat1 = S_hat1, S_hat2 = S_hat2, L_hat1 = L_hat1, L_hat2 = L_hat2, sse = sse))
}

detect.sparse <- function(data, skip = 50){
    sse <- c()
    n <- dim(data)[1]
    d <- dim(data)[2]
    sse <- c()
    for(t in skip:(T-skip)){
        ### separate into two datasets: [1,t] and [t+1, N]
        ## segment 1:
        seg.data1 <- data[1:t,]
        fit1 <- fitVAR(seg.data1, p=1, penalty = "ENET", alpha = 1, nfolds = 5, method = "cv", intercept = FALSE)
        sse1 <- (1/n) * sum(fit1$residuals^2)
        
        ## segment 2:
        seg.data2 <- data[(t+1):(T),]
        fit2 <- fitVAR(seg.data2, p=1, penalty = "ENET", alpha = 1, nfolds = 5, method = "cv", intercept = FALSE)
        sse2 <- (1/n) * sum(fit2$residuals^2)
        
        ## total sse
        sse <- c(sse, sse1 + sse2)
        #print(t-skip+1)
    }
    
    idx <- which.min(sse) + skip
    
    ### refit the transition matrix
    seg.data1 <- data[1:idx,]
    seg.data2 <- data[(idx+1):T,]
    fit1 <- fitVAR(seg.data1, p=1, penalty = "ENET", alpha = 1, nfolds = 5, method = "cv", intercept = FALSE)
    fit2 <- fitVAR(seg.data2, p=1, penalty = "ENET", alpha = 1, nfolds = 5, method = "cv", intercept = FALSE)
    
    seg1_mat <- fit1$A[[1]]
    seg2_mat <- fit2$A[[1]]
    
    return(list(cp = idx, S_hat1 = seg1_mat, S_hat2 = seg2_mat, sse = sse))
}

fixtuning.detect.LpS <- function(data, lambda, mu, skip = 50){
    sse <- c()
    n <- dim(data)[1]
    for(t in (skip+1):(n-skip)){
        ###### segmentation of dataset ######
        seg1 <- data[1:t,]
        seg2 <- data[((t+1):n),]
        
        ### fitting the left-side dataset
        X1 <- seg1[1:(t-1),]
        Y1 <- seg1[2:t,]
        fit1 <- fista.LpS(X1, Y1, lambda[1], mu[1], 100, TRUE, diag(k))
        
        ### fitting the right-side dataset
        X2 <- seg2[1:(n-t-1),]
        Y2 <- seg2[2:(n-t),]
        fit2 <- fista.LpS(X2, Y2, lambda[2], mu[2], 100, TRUE, diag(k))
        
        ###### find the residual and full SSE ######
        phi.hat1 <- t(fit1$lr.comp + fit1$sparse.comp)
        phi.hat2 <- t(fit2$lr.comp + fit2$sparse.comp)
        pred.err1 <- norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        pred.err2 <- norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        sse <- c(sse, pred.err1 + pred.err2)
    }
    idx <- which.min(sse) + skip
    
    ####### refit the L+S model to see the estimators and evaluate #######
    seg1 <- data[1:idx,]
    seg2 <- data[(idx+1):n,]
    
    ### fitting the left-side dataset
    X1 <- seg1[1:(idx-1),]
    Y1 <- seg1[2:idx,]
    fit1 <- fista.LpS(X1, Y1, lambda[1], mu[1], 100, TRUE, diag(k))
    
    ### fitting the right-side dataset
    X2 <- seg2[1:(n-idx-1),]
    Y2 <- seg2[2:(n-idx),]
    fit2 <- fista.LpS(X2, Y2, lambda[2], mu[2], 100, TRUE, diag(k))
    
    phi.hat1 <- t(fit1$lr.comp + fit1$sparse.comp)
    phi.hat2 <- t(fit2$lr.comp + fit2$sparse.comp)
    return(list(cp = idx, S_hat1 = t(fit1$sparse.comp), S_hat2 = t(fit2$sparse.comp), 
                L_hat1 = t(fit1$lr.comp), L_hat2 = t(fit2$lr.comp), sse = sse))
}

#########################################################################################################
### recording vectors
model.times <- rep(0, reps)
model.cps <- rep(0, reps)
model.sparse <- vector('list', reps) 
model.lowrank <- vector('list', reps)
model.curve <- vector('list', reps)

alter.times <- rep(0, reps)
alter.cps <- rep(0, reps)
alter.mat <- vector('list', reps)
alter.curve <- vector('list', reps)

#############################################################################################################
#################### Starting epoches: use minimizing MSPE to find the tuning parameters ####################
#############################################################################################################
for(epoch in 1:reps){
    ### generating data
    set.seed(epoch)
    skip <- 50
    e.sigma <- as.matrix(0.01*diag(k));
    try = var.sim.break(T, arlags = seq(1, p.t, 1), malags = NULL, phi = phi_full, sigma = e.sigma, brk = brk)
    data <- try$series
    data <- as.matrix(data)
    # MTSplot(data)
    
    ######################## Use grid search to find L+S result ########################
    lambda1.seq <- seq(0.1, 0.5, length.out = 5)
    lambda2.seq <- seq(0.25, 0.5, length.out = 5)

    mu1.seq <- seq(1.05, 2.5, length.out = 5)
    mu2.seq <- seq(0.75, 2.5, length.out = 5)

    ptm <- proc.time()
    result <- detect.LpS(data, lambda1.seq, lambda2.seq, mu1.seq, mu2.seq, skip = 50)
    model.times[epoch] <- (proc.time() - ptm)[3]

    model.cps[epoch] <- result$cp + 1
    print(paste("Estimated change point:", result$cp))

    model.lowrank[[epoch]] <- cbind(result$L_hat1, result$L_hat2)
    model.sparse[[epoch]] <- cbind(result$S_hat1, result$S_hat2)
    model.curve[[epoch]] <- result$sse
    
    ######################## Use sparse model to find result ########################
    ptm <- proc.time()
    fit <- detect.sparse(data, skip = 50)
    alter.times[epoch] <- (proc.time() - ptm)[3]

    alter.cps[epoch] <- fit$cp
    print(paste("Estimated change point: ", fit$cp))

    alter.mat[[epoch]] <- cbind(fit$S_hat1, fit$S_hat2)
    alter.curve[[epoch]] <- fit$sse
    ################################################################################################
    print(c(qr(result$L_hat1)$rank, qr(result$L_hat2)$rank))
    print(paste("==================================================", epoch))
}

#########################################################################################################
### saving the results and analyze locally
save(model.times, file = 'Model_times.RData')
save(model.cps, file = 'Model_cps.RData')
save(model.sparse, file = 'Model_sparse.RData')
save(model.lowrank, file = 'Model_lowrank.RData')
save(model.curve, file = "Model_curve.RData")

save(alter.times, file = 'Alter_times.RData')
save(alter.cps, file = 'Alter_cps.RData')
save(alter.mat, file = 'Alter_mat.RData')
save(alter.curve, file = 'Alter_curve.RData')
