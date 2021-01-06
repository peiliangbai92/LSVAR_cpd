fista.sparse <- function(A, b, lambda, d, niter, backtracking = TRUE, phi.true){
  tnew = t <- 1
  x <- matrix(0, d, d)
  xnew <- x
  y <- x
  AtA <- t(A) %*% A
  Atb <- t(A) %*% b
  
  obj.val = rel.err <- c()
  if(backtracking == TRUE){
    L <- norm(A, "2")^2 / 5
    eta <- 2
  }else{
    L <- norm(A, "2")^2
  }
  for(i in 1:niter){
    if(backtracking == TRUE){
      L.bar <- L
      flag <- FALSE
      while(flag == FALSE){
        prox <- prox.func(y, A, b, L.bar, lambda, AtA, Atb)
        if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
          flag <- TRUE
        }else{
          L.bar <- L.bar * eta
        }
      }
      L <- L.bar
    }
    x <- xnew
    xnew <- prox
    t <- tnew
    tnew <- (1 + sqrt(1 + 4*t^2)) / 2
    y <- xnew + ((t - 1) / tnew) * (xnew - x)
    
    obj.val <- c(obj.val, f.func(xnew, A, b) + g.func(xnew, lambda))
    rel.err <- c(rel.err, norm(xnew - phi.true, "F") / norm(phi.true, "F"))
  }
  return(list(phi.hat = xnew, obj.vals = obj.val, rel.err = rel.err))
}

shrinkage <- function(y, tau){
  z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  for(i in 1:nrow(y)){
    for(j in 1:ncol(y)){
      if(abs(y[i,j]) > tau){
        z[i,j] <- sign(y[i,j]) * (y[i,j] - tau)
      }else{
        z[i,j] <- sign(y[i,j]) * (tau/1.5)
      }
    }
  }
  return(z)
}
f.func <- function(x, A, b){
  return(0.5 * norm(A %*% x - b, "F")^2)
}
gradf.func <- function(x, AtA, Atb){
  return(AtA %*% x - Atb)
}
g.func <- function(x, lambda){
  return(lambda*sum(x))
}
Q.func <- function(x, y, A, b, L, AtA, Atb){
  return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + 0.5 * L * norm(x - y, "F")^2)
}
prox.func <- function(y, A, b, L, lambda, AtA, Atb){
  Y <- y - (1 / L) * gradf.func(y, AtA, Atb)
  return(shrinkage(Y, 2*lambda / L))
}
