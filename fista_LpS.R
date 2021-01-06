###### FISTA low-rank plus sparse estimating functions ######
fista.LpS <- function(A, b, lambda, mu, niter, backtracking = TRUE, x.true){
    tnew = t <- 1
    p <- dim(A)[2]
    x1 <- matrix(0, nrow = p, ncol = p)
    xnew1 = xnew2 <- x1
    y1 = y2 <- x1
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b
    
    if(backtracking == TRUE){
        L <- norm(A, "F")^2 / 5
        gamma <- 2
    }else{
        L <- norm(A, "F")^2
    }
    
    obj.val = rel.err <- c()
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            found <- FALSE
            while(found == FALSE){
                y <- y1 + y2
                prox1 <- prox.sparse.func(y1, y, A, b, 2*L.bar, lambda, AtA, Atb)
                prox2 <- prox.nuclear.func(y2, y, A, b, 2*L.bar, mu, AtA, Atb)
                
                ### Restricted solution space
                for(j in 1:p){
                    for(k in 1:p){
                        if(abs(prox2[j,k]) > 0.25){
                            prox2[j,k] <- 0.25 * sign(prox2[j,k])
                        }
                    }
                }
                
                prox <- prox1 + prox2
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    found <- TRUE
                }else{
                    L.bar <- L.bar * gamma
                }
            }
            L <- L.bar
        }
        x1 <- xnew1 
        x2 <- xnew2
        xnew1 <- prox1
        xnew2 <- prox2
        t = tnew
        tnew <- (1 + sqrt(1 + 4*t^2))/2
        y1 <- xnew1 + (t - 1) / tnew * (xnew1 - x1)
        y2 <- xnew2 + (t - 1) / tnew * (xnew2 - x2)
        xnew <- xnew1 + xnew2
        
        obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew1, lambda) + nuclear.pen(xnew2, mu))
        rel.err <- c(rel.err, norm(xnew - x.true, "F") / norm(x.true, "F"))
    }
    return(list(sparse.comp = xnew1, lr.comp = xnew2, obj.val = obj.val, rel.err = rel.err))
}
shrinkage <- function(y, tau){
    z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for(i in 1:nrow(y)){
        for(j in 1:ncol(y)){
            z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
        }
    }
    return(z)
}
shrinkage.lr <- function(y, tau){
    z <- rep(0, length(y))
    for(i in 1:length(y)){
        z[i] <- sign(y[i]) * max(0, abs(y[i]) - tau)
    }
    return(z)
}
gradf.func <- function(x, AtA, Atb){
    return(AtA %*% x - Atb)
}
nuclear.pen <- function(x, lambda){
    d <- svd(x)$d
    return(lambda * sum(d))
}
sparse.pen <- function(x, lambda){
    return(lambda*sum(x))
}
f.func <- function(x, A, b){
    return(0.5 * norm(A %*% x - b, "F")^2)
}
Q.func <- function(x, y, A, b, L, AtA, Atb){
    return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + 0.5 * L * norm(x - y, "F")^2)
}
prox.nuclear.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    d <- shrinkage.lr(svd(Y)$d, 2*lambda / L)
    return(svd(Y)$u %*% diag(d) %*% t(svd(Y)$v))
}
prox.sparse.func <- function(w1, y, A, b, L, lambda, AtA, Atb){
    Y <- w1 - (1 / L) * gradf.func(y, AtA, Atb)
    return(shrinkage(Y, 2*lambda / L))
}
obj.func <- function(x.lr, x.sparse, A, b, lambda, mu){
    ### x.sparse is a list
    m <- length(x.sparse)
    loss <- 0
    for(i in 1:m){
        loss <- loss + f.func((x.lr[[i]] + x.sparse[[i]]), A, b) + sparse.pen(x.sparse[[i]], lambda) + nuclear.pen(x.lr[[i]], mu)
    }
    return(loss)
}
