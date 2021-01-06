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
source("fista_LpS.R")
########################### dynamical programming model ##############################
L.n <- function(data, s, e, lambda, mu){
    temp_data <- data[s:e, ]
    n.temp <- dim(temp_data)[1]
    p <- dim(temp_data)[2]
    
    X <- temp_data[1:(n.temp-1),]
    y <- temp_data[2:n.temp, ]
    try <- fista.LpS(X, y, lambda, mu, niter = 50, backtracking = TRUE, diag(p))
    est.coef <- t(try$lr.comp) + t(try$sparse.comp)
    pred.error <- y - X %*% t(est.coef)
    
    res <- sum(pred.error^2)
    return(res)
}
dp.detect <- function(data, s, t, flag, lambda, mu){
    ret <- c()
    n <- dim(data)[1]
    p <- dim(data)[2]
    
    ### penalty tuning parameter
    gamma <- log(n)*p*0.05
    
    ### use dynamical programming method
    while(s < n-1){
        s <- s+1
        while(t < n && flag == 0){
            t <- t+1
            L.temp <- c()
            L.total <- ifelse(t-s >= gamma, L.n(data, s, t, lambda, mu), 0)
            for(l in (s+2):(t-2)){
                ### separately consider the left, right segments for [s,l] and [l+1,t], and total [s,t]
                L.left <- ifelse(t-s >= gamma, L.n(data, s, l, lambda, mu), 0)
                L.right <- ifelse(t-s >= gamma, L.n(data, l, t, lambda, mu), 0)
                L.temp <- c(L.temp, L.left + L.right + gamma)
            }
            if(min(L.temp) < L.total){
                local_min <- which.min(L.temp) + s
                ret <- c(ret, local_min)
                s <- local_min
                flag <- 1
            }
            print(paste("subinterval time point:", t, sep = " "))
        }
        flag <- 0
        print(paste("current time point:", s, sep = " "))
        
    }
    return(list(final.cps = ret))
}


