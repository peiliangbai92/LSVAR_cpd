#'@param data: data stream
#'@param h: window size
#'@param start: starting point
#'@param lambda.seq: sequence of tuning parameter of sparse we need to select from
#'@param mu.seq: sequence of tuning parameter of low-rank we select from
#'@param pts: first step candidate change points
#'@param omega: tuning parameter for selecting screened change points
#'@param skip: skipping boundary time points
##########################################################################################
# ts.separate <- function(data, npart = 5){
#     n <- dim(data)[1]
#     if(n < 10 && n >= 4){
#         ntest <- 2
#         ntrain <- n - ntest
#         data_train <- data[1:ntrain,]
#         data_test <- data[-(1:ntrain),]
#     }
#     else{
#         ntest <- floor(n / npart)+1
#         ntrain <- n - ntest
#         data_train <- data[1:ntrain,]
#         data_test <- data[-(1:ntrain),]
#     }
#     return(list(train = data_train, test = data_test))
# }
# 
# tuning.selection <- function(data, lambda.seq, mu.seq, npart = 3){
#     ### separate data time series into training and testing set
#     separate <- ts.separate(data, npart = npart)
#     data_train <- separate$train
#     data_test <- separate$test
#     
#     ntrain <- dim(data_train)[1]
#     ntest <- dim(data_test)[1]
#     
#     X_train <- data_train[1:(ntrain-1),]
#     Y_train <- data_train[2:ntrain,]
#     
#     X_test <- data_test[1:(ntest-1),]
#     Y_test <- data_test[2:ntest,]
#     
#     ### grid search part: find the tuning parameter that minimizes the prediction error on test set
#     n <- dim(data)[1]
#     k <- dim(data)[2]
#     grid <- matrix(0, nrow = length(lambda.seq), ncol = length(mu.seq))
#     for(r in 1:length(lambda.seq)){
#         for(c in 1:length(mu.seq)){
#             lambda <- lambda.seq[r]
#             mu <- mu.seq[c]
#             fit <- fista.LpS(X_train, Y_train, lambda, mu, niter = 100, backtracking = TRUE, diag(k))
#             if(qr(fit$lr.comp)$rank == 0 || qr(fit$lr.comp)$rank > k/2){
#                 grid[r,c] <- Inf
#             }else{
#                 residual <- Y_test - X_test %*% (fit$sparse.comp + fit$lr.comp)
#                 grid[r,c] <- (1/ntest) * norm(residual, "F")^2
#             }
#         }
#     }
#     idx <- which(grid == min(grid), arr.ind = TRUE)
#     final.lambda <- lambda.seq[idx[1]]
#     final.mu <- mu.seq[idx[2]]
#     return(list(grid = grid, lambda = final.lambda, mu = final.mu))
# }

detect.LpS <- function(data, lambda, mu, skip = 50){
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
        
        lr_est[[t-skip]] <- cbind(t(fit1$lr.comp), t(fit2$lr.comp))
        sp_est[[t-skip]] <- cbind(t(fit1$sparse.comp), t(fit2$sparse.comp))
        
        pred.err1 <- (1/n) * norm(Y1 - X1 %*% t(phi.hat1), "F")^2
        pred.err2 <- (1/n) * norm(Y2 - X2 %*% t(phi.hat2), "F")^2
        sse <- c(sse, pred.err1 + pred.err2)
        print(paste("Finished time point:", t))
    }
    
    idx <- which.min(sse) + skip
    
    L_hat1 <- lr_est[[idx-skip]][,c(1:k)]; L_hat2 <- lr_est[[idx-skip]][,-(c(1:k))]
    S_hat1 <- sp_est[[idx-skip]][,c(1:k)]; S_hat2 <- sp_est[[idx-skip]][,-(c(1:k))]
    return(list(cp = idx, S_hat1 = S_hat1, S_hat2 = S_hat2, L_hat1 = L_hat1, L_hat2 = L_hat2, sse = sse))
}

##########################################################################################
################# First step: rolling window method to select candidates #################
##########################################################################################
first.step.detect <- function(data, h, lambda, mu, skip = 3){
    n <- dim(data)[1]
    k <- dim(data)[2]
    s <- 1
    e <- s + h
    candi_cp <- c()
    while(e <= n-1){
        interval_data <- data[s:e, ]
        
        ### use single change point detection method to find the candidates
        fit <- detect.LpS(interval_data, lambda, mu, skip = skip)
        current_cp <- fit$cp
        candi_cp <- c(candi_cp, current_cp + s)
        
        s <- s + floor(0.25 * h)
        if(s + h <= n){
            e <- s + h
        }else{
            e <- n
        }
        # print(paste("Finished interval:", s, e))
    }
    first.step.cp <- candi_cp[-1]
    first.step.cp <- first.step.cp[order(first.step.cp)]
    
    ### remove repeat candidate change points
    m <- length(first.step.cp)
    pts <- c()
    for(mm in 1:(m-1)){
        if(abs(first.step.cp[mm] - first.step.cp[mm+1]) > 2){
            pts <- c(pts, first.step.cp[mm])
        }
    }
    pts <- c(pts, first.step.cp[m])  ## final points after the first rolling window selection.
    return(pts)
}

############################################################################################
################# Second step: screening with backward selection algorithm #################
############################################################################################
break.var.lps <- function(data, pts, lambda, mu){
    n <- dim(data)[1]
    k <- dim(data)[2]
    m <- length(pts)
    L.n <- rep(0, m+1)
    if(m == 0){
        pts.temp <- c(1, n+1)
    }else{
        pts.temp <- c(1, pts, n+1)
    }
    for(mm in 1:(m+1)){
        data.temp <- data[(pts.temp[mm]):(pts.temp[mm+1]-1),]
        n.temp <- dim(data.temp)[1]
        x.temp <- data.temp[1:(n.temp-1),]
        y.temp <- data.temp[2:n.temp,]
        
        ### select tuning parameters
        # tuning.params <- tuning.selection(data.temp, lambda.seq, mu.seq, npart = npart)
        # lambda <- tuning.params$lambda
        # mu <- tuning.params$mu
        
        ### estimate the coefficient matrices and calculate the corresponding prediction error
        try <- fista.LpS(x.temp, y.temp, lambda, mu, niter = 20, backtracking = TRUE, diag(k))
        est.coef <- t(try$lr.comp) + t(try$sparse.comp)
        pred.error <- y.temp - x.temp %*% t(est.coef)
        L.n[mm] <- sum(pred.error^2)
    }
    return(list(L.n = sum(L.n)))
}

backward.selection <- function(data, pts, lambda, mu){
    n <- dim(data)[1]
    k <- dim(data)[2]
    m <- length(pts)
    L.n <- rep(0, m)
    L.n.current <- rep(0, 1)
    
    try <- break.var.lps(data, pts, lambda, mu)
    L.n.current <- try$L.n
    for(mm in 1:m){
        pts.temp <- pts[-mm]
        try <- break.var.lps(data, pts.temp, lambda, mu)
        L.n[mm] <- try$L.n
    }
    return(list(L.n = L.n, L.n.current = L.n.current))
}

second.step.detect <- function(data, pts, omega, lambda, mu){
    m <- length(pts)
    if(m == 0){
        break
    }
    mm <- 0; ic <- 0
    while(mm < m){
        # print(mm)
        mm <- mm + 1
        try <- backward.selection(data, pts, lambda, mu)
        L.n <- try$L.n; L.n.curr <- try$L.n.current
        if(min(L.n) + (m-1)*omega >= L.n.curr + m*omega){
            ic <- L.n.curr + (m - mm + 1) * omega
            break
        }else{
            pts <- pts[-which(L.n == min(L.n))]
            print(pts)
        }
    }
    return(list(pts = pts, ic = ic))
}


