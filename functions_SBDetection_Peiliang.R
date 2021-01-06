

######## FUNCTIONS #################


var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k = nrow(sigma)
  m <- length(brk)
  nT = nobs + skip
  at = rmvnorm(nT, rep(0, k), sigma)
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  nma = length(malags)
  if (nma > 0) {
    malags = sort(malags)
    q = malags[nma]
  }
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }
  
  if (m > 1){
    
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }
  
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}

pred <- function(Y,phi,p,T,k,h){
  concat.Y <- matrix(0,k,p+h); concat.Y[,1:p] <- Y[,(T-p+1):T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,p+j-i];}
    concat.Y[,p+j] <- temp; 
  }
  return(as.matrix(concat.Y[,p+h]))
}

pred.block <- function(Y,phi,p,T,k,h){
  concat.Y <- matrix(0,k,p+h); concat.Y[,1:p] <- Y[,(T-p+1):T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,p+j-i];}
    concat.Y[,p+j] <- temp; 
  }
  return(as.matrix(concat.Y[,(p+1):(p+h)]))
}

pred.block.new <- function(Y,phi.1,phi.2,p,T,k,h,sep){
  concat.Y <- matrix(0,k,T+h); concat.Y[,1:T] <- Y[,1:T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    if( j + T <= sep ){
      for (i in 1:p){temp <- temp +  phi.1[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    }
    if( j + T > sep ){
      for (i in 1:p){temp <- temp +  phi.2[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    }
    # for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    concat.Y[,T+j] <- temp; 
  }
  return(as.matrix(concat.Y[,(T+1):(T+h)]))
}

pred.block.new.local <- function(Y,phi.1,phi.2,p,T,k,h,sep){
  concat.Y <- matrix(0,k,p+h); concat.Y[,1:p] <- Y[,(T-p+1):T];
  for ( j in 1:h){
    temp <- matrix(0,k,1);
    if( j + T <= sep ){
      for (i in 1:p){temp <- temp +  phi.1[,((i-1)*k+1):(i*k)]%*%concat.Y[,p+j-i];}
    }
    if( j + T > sep ){
      for (i in 1:p){temp <- temp +  phi.2[,((i-1)*k+1):(i*k)]%*%concat.Y[,p+j-i];}
    }
    # for (i in 1:p){temp <- temp +  phi[,((i-1)*k+1):(i*k)]%*%concat.Y[,T+j-i];}
    concat.Y[,p+j] <- temp; 
  }
  return(as.matrix(concat.Y[,(p+1):(p+h)]))
}

soft <- function(L,weight,lambda){
  for (i in 1:length(L[1,])){
    lambda <- lambda*(1+weight[i])
    if ( L[i] > lambda){L[i] <- L[i] - lambda}
    if ( L[i] < -lambda){L[i] <- L[i] + lambda}
    if ( abs(L[i]) <= lambda){L[i] <- 0}
  }
  return(L)
}

soft.full <- function(L,lambda,k,p,n){
  
  # for(kk in 1:n){
  #   temp <- L[,((kk-1)*k*p+1):(kk*k*p)];
  #   nrm <- sum(abs(temp))
  #   if ( nrm <= lambda){ L[,((kk-1)*k*p+1):(kk*k*p)] <- matrix(0,k,k*p)  }
  #   if ( nrm > lambda) { L[,((kk-1)*k*p+1):(kk*k*p)] <- L[,((kk-1)*k*p+1):(kk*k*p)] - matrix((lambda/(p*k^2)),k,k*p); }
  # }
  # nrm <- sum(abs(L))
  # if ( nrm <= lambda){ L <- matrix(0,k,k*p)  }
  # if ( nrm > lambda) { L <- L - matrix((lambda/(p*k^2)),k,k*p);  }
  
  
  for (i in 1:length(L[1,])){
    for(j in 1:length(L[,1])){
      if ( L[j,i] > lambda){L[j,i] <- L[j,i] - lambda}
      if ( L[j,i] < -lambda){L[j,i] <- L[j,i] + lambda}
      if ( abs(L[j,i]) <= lambda){L[j,i] <- 0}
    }
  }
  return(L)
}

soft.group <- function(L,weight,group.lag,lambda){
  L <- L/(1+weight);
  for (i in 1:length(group.lag[,1])){
    temp <- group.lag[i,]; temp <- temp[which(temp!=0)];
    L[temp] <- (max(0,1-lambda/(sqrt(sum((L[temp])^2)))))*L[temp]
  }
  return(L*(1+weight))
}

var.lasso.brk <- function(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4)){

  k <- length(data[1,]); T <- length(data[,1]); T.1 <- T;
  iter <- matrix(0,k,length(lambda));
  phi.hat <- matrix(0,k,k*p); phi.hat.fista <- matrix(0,max.iteration,k*p);
  pred.error <- rep(0,length(lambda)); phi.hat.temp <- matrix(0,k,k*p*length(lambda));
  Y <- as.matrix(t(data)); 
  Y <- as.matrix(Y[,(p+1):T.1]); 
  # Y <- Y%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  Z <- matrix(0,k*p,T.1-p); 
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for ( i in 1:(T.1-p)){
    for ( j in 1:p){
      Z[((j-1)*k+1):(j*k),i] <- t(data[i+p-j,])
    }
  }
  step.size <- 1/( max(svd(Z)$d)  )^2; 


    for (ll in 1:length(lambda)){
      
      
      # temp <- foreach(ii=1:k, .combine=rbind, .export=c("soft")) %dopar% {
      #   l <- 2;
      #   while( l < max.iteration){
      #     l <- l+1;
      #     phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
      #     phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
      #     phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
      #     if ( max(abs(phi.new - phi.temp)) < tol) {break;}
      #     if (max(abs(phi.new - phi.temp)) > tol ) {
      #       phi.hat.fista[l,] <- phi.new;
      #       # print(l);
      #     }
      #   }
      #   phi.new
      # }
      # phi.hat.temp <- temp;
      
      
      for ( ii in 1:k){
        l <- 2;
        while( l < max.iteration){
          l <- l+1;
          phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
          phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
          phi.new <- soft(phi.new,rep(0,k*p),lambda[ll]);
          if ( max(abs(phi.new - phi.temp)) < tol) {break;}
          if (max(abs(phi.new - phi.temp)) > tol ) {
            phi.hat.fista[l,] <- phi.new;
            # print(l);
            }
        }

        # print("l="); print(l)
        iter[ii,ll] <- l;
        phi.hat.temp[ii,((ll-1)*k*p+1):(ll*k*p)] <- phi.new;
      }
          # forecast <- sapply(c(1:(T.2-T.1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,T.1,k,jjj)  )
      forecast <- matrix(0,k,T.1-p);
      forecast <- sapply(c((p):(T.1-1)), function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)],p,jjj,k,1)  )
      pred.error[ll] <- sum((t(data[(p+1):(T.1),])-forecast)^2) + 0*lambda[ll]*sum(abs(phi.hat.temp[,((ll-1)*k*p+1):(ll*k*p)]));
    }
    ll.final <- 1
    phi.hat <- phi.hat.temp[,((ll.final-1)*k*p+1):(ll.final*k*p)];
  
  


  
  
  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
}

var.break.fit <- function(method,data, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3),initial.phi = NULL){
  method.full <- c("LASSO");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]);
  iter <- matrix(0,k,length(lambda));
  n <- T - p;
  Y <- matrix(0,k*p,n);
  for( i in p:(T-1)){
    Y[,(i-p+1)] <- sapply(c(1:p), function(jjj) data[i-jjj+1,]  )
  }
  
  C <- vector("list",n);
  # for(jjj in 1:n){C[[jjj]] <- as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(data[jjj+p,],1,k)); }
  C <- lapply(c(1:n), function(jjj) as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(data[jjj+p,],1,k))    );
  C.sum <- matrix(0,k*p*n,k);
  C.sum[1:(k*p),] <- C[[1]];
  for(i in 2:n){C.sum[((i-1)*k*p+1):(i*k*p),] <- C.sum[((i-2)*k*p+1):((i-1)*k*p),] + C[[i]] }
  C.sum.new <- matrix(0,k*p*n,k);
  C.sum.new[1:(k*p),] <- C.sum[((n-1)*k*p+1):(n*k*p),];
  for(i in 2:n){C.sum.new[((i-1)*k*p+1):(i*k*p),] <- C.sum[((n-1)*k*p+1):(n*k*p),] - C.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  
  D <- vector("list",n);
  D <- lapply(c(1:n), function(jjj) as.matrix(Y[,jjj],k*p,1)%*%t(as.matrix(Y[,jjj],k*p,1))    );
  D.sum <- matrix(0,k*p*n,k*p);
  D.sum[1:(k*p),] <- D[[1]];
  for(i in 2:n){D.sum[((i-1)*k*p+1):(i*k*p),] <- D.sum[((i-2)*k*p+1):((i-1)*k*p),] + D[[i]] }
  D.sum.new <- matrix(0,k*p*n,k*p);
  D.sum.new[1:(k*p),] <- D.sum[((n-1)*k*p+1):(n*k*p),];
  # D.sum.new[(k*p+1):(n*k*p),] <- sapply(c(2:n), function(jjj) D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((jjj-2)*k*p+1):((jjj-1)*k*p),]    )
  for(i in 2:n){D.sum.new[((i-1)*k*p+1):(i*k*p),] <- D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  D.sum.new.inv <- matrix(0,k*p*n,k*p);
  for(jjj in 1:n){D.sum.new.inv[((jjj-1)*k*p+1):(jjj*k*p),] <- solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),] +  (tol)*diag(k*p) );   }
  # D.sum.new.inv <- sapply(c(1:n), function(jjj)  solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),])   )
  
  D.new <- matrix(0,n,p*k);
  for(w in 1:n){D.new[w,] <- as.matrix((as.matrix(Y[,w],k*p,1))^2,p*k,1);  }
  
  
  # D.new <- rep(0,n);
  # D.new <- sapply(c(1:n), function(jjj)   sum(  (as.matrix(Y[,jjj],k*p,1))^2  )  );
  # 
  
  phi.hat <- matrix(0,k,k*p*n);
  if (!is.null(initial.phi)){phi.hat <- initial.phi;}
  active <- rep(0,n);
  active <- sapply(c(1:n), function(jjj) if ( sum((phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])^2) != 0  ) {jjj} else {0}   )
  active <- active[which(active!=0)]
  phi.new <- matrix(0,k,k*p*n);
  # phi.temp <- phi.hat;

  
  
  # X <- matrix(0,n,n*k*p);
  # for ( i in 1:n){
  #   for(j in 1:i){
  #     X[i,((j-1)*p*k+1):(j*p*k)] <-  t(as.matrix(Y[,i],k*p,1));
  #   }
  # }
  # step.size <- 1/( max(svd(X)$d)  )^2;
  # tol <- 0.5*step.size;
  # step.size <- (0.25)*tol;
  # step.size <- 0.01;
  
  # print(step.size)
  
  
  
  

  
   if (method == 'LASSO'){
    for (ll in 1:length(lambda)){
        l <- 2;

        while( l < max.iteration){
          l <- l+1; 
          phi.compare <- phi.hat;
          # cl <- makePSOCKcluster(2);
          # registerDoParallel(cl);
          # foreach (ii=1:n) %dopar%  {
          for (ii in 1:n){
            
              E <- vector("list",n);
              E <- lapply(c(1:n), function(jjj) D.sum.new[((max(jjj,ii)-1)*k*p+1):(max(jjj,ii)*k*p),]%*%t(phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])   )
              E <- Reduce("+",E);
              E <- E - (1)*D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%t(phi.hat[,((ii-1)*k*p+1):(ii*k*p)]);
            
              S <- C.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]  - E;
              
              
            # B <- t(S)%*%D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%S
            # B <- sqrt(sum(B^2));
            # if( lambda[ll] >= B) { phi.temp <- matrix(0,k,k*p);   }
            # if( lambda[ll] < B) { phi.temp <- (1 - lambda[ll]/B)*D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%S   }
              
              S <- soft.full(S,lambda[ll],k,p,n)
              # phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
              phi.temp <- matrix(0,k*p,k);
              for(ww in 1:(k*p)){phi.temp[ww,] <- S[ww,]/(sum(D.new[(ii):(n),ww]))    }
              
              
              # phi.temp <- ( S )/sum(as.vector(D.new[(ii):(n)]))
              phi.temp <- t(phi.temp);
              
              phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
              phi.new[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
            
            
            # phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
            # phi.temp <- t(phi.temp);
            # phi.temp <- soft.full(phi.temp,lambda[ll],k,p,n)
            # phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
            # phi.new[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
          }

          
          if ( max(abs(phi.new - phi.compare)) < tol) {
            print(max(abs(phi.new - phi.compare)));
            break;} 
          if (max(abs(phi.new - phi.compare)) > tol ) {
              phi.hat <- phi.new;    
              print(max(abs(phi.new - phi.compare)))
          }
          if ( max(abs(phi.new - phi.compare)) > 10^5) {
            print("NOT CONVERGED");
            break;}
        }
        # print("l="); print(l);

    }

  }
  
  # stopCluster(cl)
  
  # for(i in 1:(n*p*k)){
  #   for( j in 1:k){
  #     if ( abs (phi.hat[j,i] <= tol) ){phi.hat[j,i] <- 0;}
  #   }
  # }
  return(list(phi.hat = phi.hat, iter = iter))
}


var.break.fit.block <- function(method,data, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3),initial.phi = NULL, blocks, cv.index){
  method.full <- c("LASSO");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]); n.new <- length(blocks) - 1;
  iter <- matrix(0,k,length(lambda));
  n <- T - p;
  
  Y.b <- vector("list",n.new); y.b <- vector("list",n.new);
  y.b <- lapply(c(1:n.new), function(jjj)   data[(blocks[jjj]+1):(blocks[jjj+1]),]  );
  y.b[[1]] <- data[(p+1):(blocks[2]),];
  
  Y <- matrix(0,k*p,n);
  for( i in p:(T-1)){
    Y[,(i-p+1)] <- sapply(c(1:p), function(jjj) data[i-jjj+1,]  )
  }
  
  Y.b[[1]] <- Y[,(1):(blocks[2]-p)];
  Y.b[2:n.new] <- lapply(c(2:(n.new)), function(jjj)   Y[,(blocks[jjj]+1-p):(blocks[jjj+1]-p)]  );
  
  cv.l <- length(cv.index);
  if( cv.l > 0){
    for(t.1 in 1:cv.l){
      tt <- length(y.b[[cv.index[t.1]]][,1]);
      y.b[[cv.index[t.1]]] <- y.b[[cv.index[t.1]]][-tt,];
      Y.b[[cv.index[t.1]]] <- Y.b[[cv.index[t.1]]][,-tt];
    }
  }
  # Y.b <- lapply(c(1:(n.new-1)), function(jjj)   Y[,(blocks[jjj]+1):(blocks[jjj+1])]  );
  # Y.b[[n.new]] <- Y[,(blocks[n.new]+1):(blocks[n.new+1])];
  
  C <- vector("list",n.new);
  C <- lapply(c(1:n.new), function(jjj) as.matrix( Y.b[[jjj]]%*%y.b[[jjj]]  )    );
  C.sum <- matrix(0,k*p*n.new,k);
  C.sum[1:(k*p),] <- C[[1]];
  for(i in 2:n.new){C.sum[((i-1)*k*p+1):(i*k*p),] <- C.sum[((i-2)*k*p+1):((i-1)*k*p),] + C[[i]] }
  C.sum.new <- matrix(0,k*p*n.new,k);
  C.sum.new[1:(k*p),] <- C.sum[((n.new-1)*k*p+1):(n.new*k*p),];
  for(i in 2:n.new){C.sum.new[((i-1)*k*p+1):(i*k*p),] <- C.sum[((n.new-1)*k*p+1):(n.new*k*p),] - C.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  
  D <- vector("list",n.new);
  D <- lapply(c(1:n.new), function(jjj) as.matrix( Y.b[[jjj]]%*%t(Y.b[[jjj]])  )     );
  D.sum <- matrix(0,k*p*n.new,k*p);
  D.sum[1:(k*p),] <- D[[1]];
  for(i in 2:n.new){D.sum[((i-1)*k*p+1):(i*k*p),] <- D.sum[((i-2)*k*p+1):((i-1)*k*p),] + D[[i]] }
  D.sum.new <- matrix(0,k*p*n.new,k*p);
  D.sum.new[1:(k*p),] <- D.sum[((n.new-1)*k*p+1):(n.new*k*p),];
  # D.sum.new[(k*p+1):(n*k*p),] <- sapply(c(2:n), function(jjj) D.sum[((n-1)*k*p+1):(n*k*p),] - D.sum[((jjj-2)*k*p+1):((jjj-1)*k*p),]    )
  for(i in 2:n.new){D.sum.new[((i-1)*k*p+1):(i*k*p),] <- D.sum[((n.new-1)*k*p+1):(n.new*k*p),] - D.sum[((i-2)*k*p+1):((i-1)*k*p),] }
  D.sum.new.inv <- matrix(0,k*p*n.new,k*p);
  for(jjj in 1:n.new){D.sum.new.inv[((jjj-1)*k*p+1):(jjj*k*p),] <- solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),] +  1*(1*10^(-6))*diag(k*p) );   }
  # D.sum.new.inv <- sapply(c(1:n), function(jjj)  solve(D.sum.new[((jjj-1)*k*p+1):(jjj*k*p),])   )
  
  # D.new <- vector("list",n.new);
  # D.new <- rep(0,n.new);
  # D.new <- sapply(c(1:n.new), function(jjj) sum(diag(Y.b[[jjj]]%*%t(Y.b[[jjj]])))     );
  D.new <- matrix(0,n.new,p*k);
  for(w in 1:n.new){D.new[w,] <- as.matrix(diag(Y.b[[w]]%*%t(Y.b[[w]])),1,p*k);  }
  # print(D.new)
  # D.new <- sapply(c(1:n.new), function(jjj) as.matrix(diag(Y.b[[jjj]]%*%t(Y.b[[jjj]])),p*k,1)     );
  # print(cumsum(D.new))
  
  
  phi.hat <- 0.00+matrix(0,k,k*p*n.new);
  if (!is.null(initial.phi)){phi.hat <- initial.phi;}
  active <- rep(0,n.new);
  active <- sapply(c(1:n.new), function(jjj) if ( sum((phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])^2) != 0  ) {jjj} else {0}   )
  active <- active[which(active!=0)]
  phi.new <- matrix(0,k,k*p*n.new);
  # phi.temp <- phi.hat;
  
  flag <- 0;
  if (method == 'LASSO'){
    for (ll in 1:length(lambda)){
      # if (flag == 1){ ll <- max(ll-1,1);}
      l <- 2;
      
      while( l < max.iteration){
        if(l == floor(0.5*max.iteration)){tol <- (3/2)*tol;}
        if(l == floor(0.75*max.iteration)){tol <- (4/3)*tol;}
        l <- l+1; 
        phi.compare <- phi.hat;

        # phi.new <- foreach(ii=1:n.new, .combine=cbind, .export=c("soft.full")) %dopar% {
        #   E <- vector("list",n.new);
        #   E <- lapply(c(1:n.new), function(jjj) D.sum.new[((max(jjj,ii)-1)*k*p+1):(max(jjj,ii)*k*p),]%*%t(phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])   )
        #   E <- Reduce("+",E);
        #   E <- E - D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%t(phi.hat[,((ii-1)*k*p+1):(ii*k*p)]);
        # 
        #   S <- C.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]  - E;
        # 
        #   S <- soft.full(S,lambda[ll],k,p,n.new)
        #   # phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
        #   # phi.temp <- t(phi.temp);
        #   # phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
        #   # phi.temp
        #   t(D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S ))
        # }
        # 
        # phi.hat <- phi.new;
        

        for (ii in 1:n.new){

          E <- vector("list",n.new);
          E <- lapply(c(1:n.new), function(jjj) D.sum.new[((max(jjj,ii)-1)*k*p+1):(max(jjj,ii)*k*p),]%*%t(phi.hat[,((jjj-1)*k*p+1):(jjj*k*p)])   )
          E <- Reduce("+",E);
          E <- E - D.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%t(phi.hat[,((ii-1)*k*p+1):(ii*k*p)]);

          S <- C.sum.new[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]  - E;
          
          # if(ii == 1) {print(max(abs(S)))}
          S <- soft.full(S,lambda[ll],k,p,n.new)
          # if(ii == 1) {print("SOFT"); print(max(abs(S)))}
          phi.temp <- D.sum.new.inv[((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p),]%*%( S )
          # phi.temp <- ( S )/sum(as.vector(D.new[(ii):(n.new)]))
          prod <- matrix(0,k*p,k*p);
          for(ww in 1:(k*p)){prod[ww,ww] <- 1/(sum(D.new[(ii):(n.new),ww]))    }
          
          # if(lambda<100){
          # print("D.new"); print(min(diag(prod)));
          # }
          
          # phi.temp <- matrix(0,k*p,k);
          # for(ww in 1:(k*p)){phi.temp[ww,] <- S[ww,]/(sum(D.new[(ii):(n.new),ww]))    }
          # if(ii == 1 && l == 3){print(diag(prod))}
          # phi.temp <- S;
          # phi.temp <- ( S )/sum(as.vector(D.new[(ii):(n.new)]))
          
          # phii <- as.vector(abs(phi.temp));
          # tune.lasso <- quantile(phii[which(phii!=0)],0.1);
          # if(lambda > 100 && lambda <115){
          # print("tune.lasso"); print(tune.lasso);
          # }
          

          # phi.temp <- soft.full(phi.temp,0.05,k,p,n.new)
          
          
          
          
          phi.temp <- t(phi.temp);

          

          phi.hat[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;
          phi.new[,((max(ii,ii)-1)*k*p+1):(max(ii,ii)*k*p)] <- phi.temp;

        }
        
        # if ( max(abs(phi.new - phi.compare)) > 10^3){phi.new <- soft.full(phi.new,0.5,k,p,n.new);}
        phi.new <- soft.full(phi.new,0.05,k,p,n.new);
        
        
        
        # max(abs(phi.new))
        # max(abs(phi.compare))
        
        if ( max(abs(phi.new - phi.compare)) < tol) {
          # flag <- 0;
          # print(max(abs(phi.new - phi.compare)));
          # print(max(abs(phi.new)));
          break;} 
        if (max(abs(phi.new - phi.compare)) > tol ) {
          phi.hat <- phi.new;    
          # print(max(abs(phi.new - phi.compare)))
        }
        if ( max(abs(phi.new - phi.compare)) > 10^5) {
          # print(max(abs(phi.new)));
          print("NOT CONVERGED");
          # lambda <- 2*lambda;
          flag <- 1;
          break;}
      }
      # print("l="); print(l);
      
    }
    
  }
  
  # stopCluster(cl)
  
  # for(i in 1:(n*p*k)){
  #   for( j in 1:k){
  #     if ( abs (phi.hat[j,i] <= tol) ){phi.hat[j,i] <- 0;}
  #   }
  # }
  return(list(phi.hat = phi.hat, iter = iter, flag = flag))
}



break.var <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts ){
  k <- length(data[1,]); T <- length(data[,1]); m <- length(pts); L.n <- rep(0,m+1);
  if ( m == 0) { pts.temp <- c(1,T+1);}
  if ( m > 0){  pts.temp <- rep(0,m+2); pts.temp[1] <- 1; pts.temp[m+2] <- T+1; pts.temp[(2):(m+1)] <- pts;}
  for(mm in 1:( m+1 )){
    # print("mmmmm"); print(mm)
    data.temp <- data[(pts.temp[mm]):(pts.temp[(mm+1)]-1),];
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = 10^(-4))
    L.n[mm] <- try$pred.error;
  }
  return(list(L.n = sum(L.n)))
}

backward.sel <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts){
  k <- length(data[1,]); T <- length(data[,1]);
  m <- length(pts); L.n <- rep(0,m); L.n.curr <- rep(0,1);
  
  ###### SHVAR FUN FULL PTS ##################
  try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts )
  L.n.curr <- try$L.n;
  
  for( mm in 1:m){
    pts.temp <- pts[-mm];
    
    ### SHVAR FIT FUNCTION ###################
    try <- break.var(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.temp );
    L.n[mm] <- try$L.n;
    
    
  }
  return(list(L.n = L.n, L.n.curr = L.n.curr  ))
  
} 

second.step <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega){
  m <- length(pts); if( m == 0){break;}
  mm <- 0; ic <- 0;
  while(mm < m){
    mm <- mm + 1;
    try <- backward.sel(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts);
    L.n <- try$L.n; L.n.curr <- try$L.n.curr;
    if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
      ic <- L.n.curr + (m - mm + 1)*omega;
      break;}
    if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
      pts <- pts[-which(L.n == min(L.n))];
      print(pts);
      }
    
  }
  
  
  return(list(pts = pts, ic = ic ))
}



second.step.middle <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts){
  m <- length(pts); if( m == 0){break;}
  mm <- 0; ic <- 0;
  while(mm < (m-1)){
    mm <- mm + 1;
    try <- backward.sel(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts);
    L.n <- try$L.n; L.n.curr <- try$L.n.curr;
    pts <- pts[-which(L.n == min(L.n))];
    # print(pts);
    # if(mm == (m-1)){pts.final <- pts;}
    # if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
    #   ic <- L.n.curr + (m - mm + 1)*omega;
    #   break;}
    # if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
    #   pts <- pts[-which(L.n == min(L.n))];
    #   print(pts);
    # }
    
  }
  
  
  return(list(pts = pts, ic = ic ))
}



break.var.local <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega.loc, ind ){
  k <- length(data[1,]); T <- length(data[,1]); m <- length(pts); 
  bounds <- vector("list",2*m);
  for(i in 1:m){
    bounds[[(2*i-1)]] <- c(pts[i] - omega.loc, pts[i] - 1 );
    bounds[[(2*i)]] <- c(pts[i], pts[i] + omega.loc );
  }
  # bounds.temp <- bounds;
  if(ind != 0){
    bounds[[(2*ind-1)]] <- c(pts[ind] - omega.loc, pts[ind] + omega.loc );
    bounds <- bounds[-(2*ind)];
  }
  
  m.new <- length(bounds);
  L.n <- rep(0,m.new);
  phi.local <- vector("list",m.new);
  
  if ( m == 0) { pts.temp <- c(1,T+1);}
  if ( m > 0){  pts.temp <- rep(0,m+2); pts.temp[1] <- 1; pts.temp[m+2] <- T+1; pts.temp[(2):(m+1)] <- pts;}
  r <- foreach(mm=1:m.new, .inorder = FALSE, .export = c("var.lasso.brk","soft","pred")) %dopar% {
    temp <- vector("list",2);
    data.temp <- data[(bounds[[mm]][1]):(bounds[[mm]][2]),];
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = tol)
    temp[[1]] <- try$pred.error;
    temp[[2]] <- try$phi.hat;
    temp
  }
  L.n <- as.vector(sapply(c(1:m.new), function(jj) r[[jj]][[1]]));
  phi.local <- lapply(c(1:m.new), function(jj) r[[jj]][[2]]);
  
  # for(mm in 1:( m.new )){
  #   # print("mmmmm"); print(mm)
  #   data.temp <- data[(bounds[[mm]][1]):(bounds[[mm]][2]),];
  #   try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = tol)
  #   L.n[mm] <- try$pred.error;
  #   phi.local[[mm]] <- try$phi.hat;
  # }
  
  return(list(L.n = sum(L.n), phi.local = phi.local))
}


break.var.local.new <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega.loc){
  k <- length(data[1,]); T <- length(data[,1]); m <- length(pts); 
  bounds.1 <- vector("list",2*m); bounds.2 <- vector("list",m);
  for(i in 1:m){
    bounds.1[[(2*i-1)]] <- c(pts[i] - omega.loc, pts[i] - 1 );
    bounds.1[[(2*i)]] <- c(pts[i], pts[i] + omega.loc );
    bounds.2[[(i)]] <- c(pts[i] - omega.loc, pts[i] + omega.loc );
  }

  L.n.1 <- rep(0,2*m);
  phi.local.1 <- vector("list",2*m);
  
  L.n.2 <- rep(0,m);
  phi.local.2 <- vector("list",m);
  
  if ( m == 0) { pts.temp <- c(1,T+1);}
  if ( m > 0){  pts.temp <- rep(0,m+2); pts.temp[1] <- 1; pts.temp[m+2] <- T+1; pts.temp[(2):(m+1)] <- pts;}
  
  r.1 <- foreach(mm=1:(2*m), .inorder = FALSE, .export = c("var.lasso.brk","soft","pred")) %dopar% {
    temp <- vector("list",2);
    data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]),];
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = tol)
    temp[[1]] <- try$pred.error;
    temp[[2]] <- try$phi.hat;
    temp
  }
  L.n.1 <- as.vector(sapply(c(1:(2*m)), function(jj) r.1[[jj]][[1]]));
  phi.local.1 <- lapply(c(1:(2*m)), function(jj) r.1[[jj]][[2]]);
  
  r.2 <- foreach(mm=1:(m), .inorder = FALSE, .export = c("var.lasso.brk","soft","pred")) %dopar% {
    temp <- vector("list",2);
    data.temp <- data[(bounds.2[[mm]][1]):(bounds.2[[mm]][2]),];
    try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = tol)
    temp[[1]] <- try$pred.error;
    temp[[2]] <- try$phi.hat;
    temp
  }
  L.n.2 <- as.vector(sapply(c(1:(m)), function(jj) r.2[[jj]][[1]]));
  phi.local.2 <- lapply(c(1:(m)), function(jj) r.2[[jj]][[2]]);
  
  # for(mm in 1:( m.new )){
  #   # print("mmmmm"); print(mm)
  #   data.temp <- data[(bounds[[mm]][1]):(bounds[[mm]][2]),];
  #   try <- var.lasso.brk(data = data.temp, weight = NULL, lambda, p, max.iteration = 1000, tol = tol)
  #   L.n[mm] <- try$pred.error;
  #   phi.local[[mm]] <- try$phi.hat;
  # }
  
  return(list(L.n.1 = L.n.1, L.n.2 = L.n.2, phi.local.1 = phi.local.1, phi.local.2 = phi.local.2))
}



backward.sel.local <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega.loc){
  k <- length(data[1,]); T <- length(data[,1]);
  m <- length(pts); L.n <- rep(0,m); L.n.curr <- rep(0,1);
  
  ###### SHVAR FUN FULL PTS ##################
  try <- break.var.local(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc, ind = 0 )
  L.n.curr <- try$L.n; phi.local <- try$phi.local;
  
  for( mm in 1:m){
    pts.temp <- pts[-mm];
    
    ### SHVAR FIT FUNCTION ###################
    try <- break.var.local(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc, ind = mm );
    L.n[mm] <- try$L.n;
    
    
  }
  return(list(L.n = L.n, L.n.curr = L.n.curr, phi.local = phi.local ))
  
} 

second.step.local <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega, omega.loc){
  m <- length(pts); if( m == 0){break;}
  mm <- 0; ic <- 0;
  while(mm < m){
    mm <- mm + 1;
    try <- backward.sel.local(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc);
    L.n <- try$L.n; L.n.curr <- try$L.n.curr; phi.local <- try$phi.local;
    if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
      ic <- L.n.curr + (m - mm + 1)*omega;
      break;}
    if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
      pts <- pts[-which(L.n == min(L.n))];
      print(pts);
    }
    
  }
  
  
  return(list(pts = pts, ic = ic, phi.local = phi.local ))
}

second.step.local.new <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega, omega.loc, al){
  m <- length(pts); if( m == 0){break;}
  try <- break.var.local.new(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc);
  # print("DONE")
  L.n.1 = try$L.n.1; L.n.2 = try$L.n.2; phi.local.1 = try$phi.local.1; phi.local.2 = try$phi.local.2; 
  omega.new <- abs((sum(L.n.2)-sum(L.n.1)));
  # print("OMEGA"); print(omega.new);
  # print("L.n.1"); print(L.n.1);
  # print("L.n.2"); print(L.n.2);
  L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2; L.n.plot <- rep(0,m+1); L.n.plot[1] <- sum(L.n.1) + 0*(m)*omega; 
  # L.n.plot[m+2] <- sum(L.n.2);
  mm <- 0; ic <- 0; add.temp <- 0; pts.full <- vector("list",m+1); pts.full[[1]] <- pts; ind.pts <- rep(0,m);
  while(mm < m){
    mm <- mm + 1;
    L.n.temp <- rep(0,length(pts));
    for(i in 1:length(pts)){
      L.n.temp[i] <- sum(L.n.1.temp) - L.n.1.temp[(2*i-1)] - L.n.1.temp[(2*i)] + L.n.2.temp[i] + 1*add.temp;
    }
    ll <- min(which.min(L.n.temp)); ind.pts[mm] <- ll;
    pts <- pts[-ll]; L.n.1.temp <- L.n.1.temp[-c(2*ll-1,2*ll)]; add.temp <- add.temp + 1*L.n.2.temp[ll]; 
    L.n.2.temp <- L.n.2.temp[-ll]; 
    L.n.plot[mm+1] <- L.n.temp[ll] + 0*(m - mm)*omega + 0*add.temp;
    pts.full[[mm+1]] <- pts;
    # print(pts);
    
    
    # try <- backward.sel.local(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc);
    # L.n <- try$L.n; L.n.curr <- try$L.n.curr; phi.local <- try$phi.local;
    # if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
    #   ic <- L.n.curr + (m - mm + 1)*omega;
    #   break;}
    # if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
    #   pts <- pts[-which(L.n == min(L.n))];
    #   print(pts);
    # }
    
  }
  ind <- 0;
  a <- as.vector(diff(L.n.plot));
  omega.new <- 1;
  for(i in 1:m){
    if(L.n.plot[m+1] - L.n.plot[i] > 0){omega.new <- max(omega.new,L.n.plot[m+1] - L.n.plot[i])}
  }
  # print(a);
  b <- quantile(a[which(a > 0)],1-al);
  # b <- (1-al)*omega.new;
  for(j in 1:m){
    if(max(a[1:j]) <= b){ind <- ind + 1;}
  }
  
  if(ind >=1){
    for(i in 1:ind){
      phi.local.1 <- phi.local.1[-c(2*ind.pts[i]-1,2*ind.pts[i])];
      
    }
  }
  
  aa <- a[a>0];
  plot( seq(1,length(aa),1), aa, type = 'o', col = "blue", main=c("Screening"))
  plot( seq(1,length(a),1), a, type = 'o', col = "blue", main=c("Screening"))
    if (ind != 0){abline(v=  ind  );}
  
  
  return(list(pts = pts.full[[ind+1]], ic = ic, phi.local = phi.local.1, L.n.plot = L.n.plot, L.n.1 = L.n.1, L.n.2 = L.n.2, omega = b  ))
}

second.step.local.cluster <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega, omega.loc, al){
    m <- length(pts); if( m == 0){break;}
    try <- break.var.local.new(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc);
    # print("DONE")
    L.n.1 = try$L.n.1; L.n.2 = try$L.n.2; phi.local.1 = try$phi.local.1; phi.local.2 = try$phi.local.2; 
    omega.new <- abs((sum(L.n.2)-sum(L.n.1)));
    pts.deleted <- rep(0,m);
    # print("OMEGA"); print(omega.new);
    # print("L.n.1"); print(L.n.1);
    # print("L.n.2"); print(L.n.2);
    L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2; L.n.plot <- rep(0,m+1); L.n.plot[1] <- sum(L.n.1) + 0*(m)*omega; 
    # L.n.plot[m+2] <- sum(L.n.2);
    mm <- 0; ic <- 0; add.temp <- 0; pts.full <- vector("list",m+1); pts.full[[1]] <- pts; ind.pts <- rep(0,m);
    while(mm < m){
        mm <- mm + 1;
        L.n.temp <- rep(0,length(pts));
        for(i in 1:length(pts)){
            L.n.temp[i] <- sum(L.n.1.temp) - L.n.1.temp[(2*i-1)] - L.n.1.temp[(2*i)] + L.n.2.temp[i] + 1*add.temp;
        }
        ll <- min(which.min(L.n.temp)); ind.pts[mm] <- ll;
        pts.deleted[mm] <- pts[ll];
        pts <- pts[-ll]; L.n.1.temp <- L.n.1.temp[-c(2*ll-1,2*ll)]; add.temp <- add.temp + 1*L.n.2.temp[ll]; 
        L.n.2.temp <- L.n.2.temp[-ll]; 
        L.n.plot[mm+1] <- L.n.temp[ll] + 0*(m - mm)*omega + 0*add.temp;
        pts.full[[mm+1]] <- pts;
        # print(pts);
        
        
        # try <- backward.sel.local(data, lambda, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts, omega.loc);
        # L.n <- try$L.n; L.n.curr <- try$L.n.curr; phi.local <- try$phi.local;
        # if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
        #   ic <- L.n.curr + (m - mm + 1)*omega;
        #   break;}
        # if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
        #   pts <- pts[-which(L.n == min(L.n))];
        #   print(pts);
        # }
        
    }
    
    # ind <- 0;
    # a <- as.vector(diff(L.n.plot));
    # omega.new <- 1;
    # for(i in 1:m){
    #   if(L.n.plot[m+1] - L.n.plot[i] > 0){omega.new <- max(omega.new,L.n.plot[m+1] - L.n.plot[i])}
    # }
    # # print(a);
    # b <- quantile(a[which(a > 0)],1-al);
    # # b <- (1-al)*omega.new;
    # for(j in 1:m){
    #   if(max(a[1:j]) <= b){ind <- ind + 1;}
    # }
    # 
    # if(ind >=1){
    #   for(i in 1:ind){
    #     phi.local.1 <- phi.local.1[-c(2*ind.pts[i]-1,2*ind.pts[i])];
    #     
    #   }
    # }
    # 
    # aa <- a[a>0];
    # plot( seq(1,length(aa),1), aa, type = 'o', col = "blue", main=c("Screening"))
    # plot( seq(1,length(a),1), a, type = 'o', col = "blue", main=c("Screening"))
    # if (ind != 0){abline(v=  ind  );}
    
    
    ind <- 0;
    b <- 1;
    print(pts.deleted)
    print(L.n.plot)
    a <- as.vector(diff(L.n.plot));
    # print(a);
    plot(sort(a),type="o");
    omega.new <- 1;
    pts.sel <- c();
    if( m == 1){
        if( abs(a/(L.n.plot[1]+0.01)) <= 0.3){pts.sel <- c();}
        if( abs(a/(L.n.plot[1]+0.01)) > 0.3){pts.sel <- pts;}
    }
    if( m > 1){
        clus.1 <- kmeans(a, centers = 1); fit.1 <- clus.1$betweenss/clus.1$totss; print(fit.1);
        if(fit.1 > 0.8){pts.sel <- c();}
        if( fit.1 <= 0.8 ){
            clus.2 <- kmeans(a, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; print(fit.2);
            loc <- clus.2$cluster;
            if( clus.2$centers[1] > clus.2$centers[2]  ){pts.sel <- pts.deleted[which(loc==1)]; b <- a[min(which(loc==1))];}
            if( clus.2$centers[1] < clus.2$centers[2]  ){pts.sel <- pts.deleted[which(loc==2)]; b <- a[min(which(loc==2))];}
        }
    }
    
    pts.sel <- sort(pts.sel);
    print("FINAL POINTS:");
    print(pts.sel);
    
    return(list(pts = pts.sel, ic = ic, phi.local = phi.local.1, L.n.plot = L.n.plot, L.n.1 = L.n.1, L.n.2 = L.n.2, omega = b  ))
}


second.step.block <- function(data, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts, omega, blocks){
  m <- length(pts); if( m == 0){break;}
  n.new <- length(blocks) - 1;
  data.new <- data[blocks[(2):(n.new+1)],];
  pts.ind <- rep(0,m);
  for(i in 1:m){
    pts.temp <- pts[i];
    for(j in 1:(n.new+1)){
      if(pts.temp == blocks[j]){pts.ind[i] <- j-1}
    }
  }
  both <- matrix(0,2,m); both[1,] <- pts; both[2,] <- pts.ind;
  mm <- 0; ic <- 0;
  while(mm < m){
    mm <- mm + 1;
    try <- backward.sel(data.new, lambda, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.ind);
    L.n <- try$L.n; L.n.curr <- try$L.n.curr;
    if ( min(L.n) + (m-1)*omega >= L.n.curr + m*omega   ) {
      ic <- L.n.curr + (m - mm + 1)*omega;
      break;}
    if ( min(L.n) + (m-1)*omega < L.n.curr + m*omega   ){
      pts <- pts[-which(L.n == min(L.n))];
      print(pts);
    }
    
  }
  
  
  return(list(pts = pts, ic = ic ))
}

plot.new <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    # if (k < 4) {
    #   par(mfcol = c(k, 1))
    #   for (j in 1:k) {
    #     plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
    #          type = "l")
    #   }
    # }

    if (k >= 1) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05;
      # plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
      #      ylim = yl)
      plot(tdx, data[, 1], xlab = "seconds", ylab = " ", type = "l", 
           ylim = yl, xaxt="n")
      axis(1, c(1,500,1000,1500,2000), c(0,50,100,150,200))
      axis(2)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}


plot.new.high <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    # if (k < 4) {
    #   par(mfcol = c(k, 1))
    #   for (j in 1:k) {
    #     plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
    #          type = "l")
    #   }
    # }
    
    if (k >= 1) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05;
      # plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
      #      ylim = yl)
      plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
           ylim = yl, xaxt="n")
      axis(1, c(1,20,40,60,80), c(0,20,40,60,80))
      axis(2)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}


MTSplot.new <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    if (k < 0) {
      par(mfcol = c(k, 1))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if (k == 0) {
      par(mfcol = c(2, 2))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if ((k > 0) && (k < 1)) {
      par(mfcol = c(3, 2), mai = c(0.3, 0.3, 0.3, 0.3))
      k1 = 6
      jcnt = 0
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l", cex.axis = 0.8)
        jcnt = jcnt + 1
        if ((jcnt == k1) && (k > 6)) {
          jcnt = 0
          cat("Hit return for more plots: ", "\n")
          readline()
        }
      }
    }
    if (k > 0) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05
      plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
           ylim = yl)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}


selection.plot <- function (data, caltime = NULL){
  if (!is.matrix(data)) 
    data = as.matrix(data)
  if (is.ts(data)) {
    plot(data)
  }
  else {
    nT = dim(data)[1]
    tdx = c(1:nT)
    if (length(caltime) > 1) 
      tdx = caltime
    k = dim(data)[2]
    if (k < 0) {
      par(mfcol = c(k, 1))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if (k == 0) {
      par(mfcol = c(2, 2))
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l")
      }
    }
    if ((k > 0) && (k < 1)) {
      par(mfcol = c(3, 2), mai = c(0.3, 0.3, 0.3, 0.3))
      k1 = 6
      jcnt = 0
      for (j in 1:k) {
        plot(tdx, data[, j], xlab = "time", ylab = colnames(data)[j], 
             type = "l", cex.axis = 0.8)
        jcnt = jcnt + 1
        if ((jcnt == k1) && (k > 6)) {
          jcnt = 0
          cat("Hit return for more plots: ", "\n")
          readline()
        }
      }
    }
    if (k > 0) {
      par(mfcol = c(1, 1))
      yl = range(data) * 1.05
      plot(tdx, data[, 1], xlab = "time", ylab = " ", type = "l", 
           ylim = yl)
      for (j in 2:k) {
        lines(tdx, data[, j], lty = j, col = j)
      }
    }
  }
  par(mfcol = c(1, 1))
}

ar.est <- function(method, data, weight = NULL, lambda, p, break.pts, r.n, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3)){
  method.full <- c("LASSO");
  if ( !(method %in% method.full) ){print("ERROR"); break; }
  
  k <- length(data[1,]); T <- length(data[,1]); T.1 <- T; m.hat <- length(break.pts) + 1;
  ind.remain <- rep(0,(2+2*length(break.pts))); ind.remain[1] <- p; ind.remain[(2+2*length(break.pts))] <- T.1;
  if(length(break.pts) >= 1){
    for(i in 1:length(break.pts)){ind.remain[(2*i)] <- break.pts[i] - r.n - 1; ind.remain[(2*i+1)] <- break.pts[i] + r.n + 1;  }
  }
  iter <- matrix(0,k,length(lambda));
  phi.hat <- matrix(0,k,k*m.hat*p); phi.hat.fista <- matrix(0,max.iteration,k*m.hat*p);
  pred.error <- rep(0,length(lambda)); phi.hat.temp <- matrix(0,k,k*m.hat*p*length(lambda)); std.res <- rep(0,length(lambda));
  
  Y <- as.matrix(t(data)); Y <- Y[,-seq(1,p,1)];
  Z <- matrix(0,k*p,T.1-p); 
  # Z <- Z%*%( diag(T-p) - matrix(1/T,T-p,T-p)  );
  for ( i in 1:(T.1-p)){
    for ( j in 1:p){
      Z[((j-1)*k+1):(j*k),i] <- t(data[i+p-j,])
    }
  }
  
  if(length(break.pts) >= 1){
    for (i in 1:length(break.pts)){
      Y <- Y[,-seq(break.pts[i]-r.n,break.pts[i]+r.n,1)];
      # Z <- Z[,-seq(break.pts[i]-r.n,break.pts[i]+r.n,1)];
    }
  }

  
  n <- length(Y[1,]);
  Z.new <- matrix(0,T.1-p,k*m.hat*p);
  del.ind <- c();
  if(length(break.pts) >= 1){
    Z.new[(1:(break.pts[1]-r.n-1)),1:(k*p)] <- t(Z[1:(k*p),(1:(break.pts[1]-r.n-1))]);
    if( m.hat > 2 ){
      for(i in 1:(m.hat-2)){
        # ind <- break.pts[i]-r.n-1;
        Z.new[((break.pts[i]+r.n+1):(break.pts[i+1]-r.n-1)),((i)*k*p+1):((i+1)*k*p)] <- t(Z[1:(k*p),((break.pts[i]+r.n+1):(break.pts[i+1]-r.n-1))]);
      }
    }
    
    Z.new[((break.pts[m.hat-1]+r.n+1):(T.1-p)),((m.hat-1)*k*p+1):((m.hat)*k*p)] <- t(Z[1:(k*p),((break.pts[m.hat-1]+r.n+1):(T.1-p))])
    del.ind <- c();
    for(i in 1:(T.1-p)){
      if (sum(Z.new[i,]^2) == 0){del.ind <- c(del.ind,i)}
    }
    Z.new <- Z.new[-del.ind,]
    Z <- t(Z.new);
  }
  

  
  

  step.size <- 1/( max(svd(Z)$d)  )^2;
  # print(step.size)
  # step.size <- 10^(-3);
  
  
  for (ll in 1:length(lambda)){
    
    
    temp <- foreach(ii=1:k, .combine=rbind, .export=c("soft")) %dopar% {
      if(ll > 1){
        phi.hat.fista[1,] <- phi.hat.temp[ii,((ll-2)*k*m.hat*p+1):((ll-1)*k*m.hat*p)];
        phi.hat.fista[2,] <- phi.hat.temp[ii,((ll-2)*k*m.hat*p+1):((ll-1)*k*m.hat*p)];
      }
      l <- 2;
      while( l < max.iteration){
        l <- l+1; 
        phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
        phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
        phi.new <- soft(phi.new,rep(0,k*m.hat*p),lambda[ll]);
        if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
        if (max(abs(phi.new - phi.temp)) > tol ) {
          phi.hat.fista[l,] <- phi.new; 
          # print(l);  
        }
      }
      phi.new
    }
    
    phi.hat.temp[,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)] <- temp;
    
    # for ( ii in 1:k){
    #   if(ll > 1){
    #     phi.hat.fista[1,] <- phi.hat.temp[ii,((ll-2)*k*m.hat*p+1):((ll-1)*k*m.hat*p)];
    #     phi.hat.fista[2,] <- phi.hat.temp[ii,((ll-2)*k*m.hat*p+1):((ll-1)*k*m.hat*p)];
    #   }
    #   l <- 2;
    #   while( l < max.iteration){
    #     l <- l+1; 
    #     phi.temp <- phi.hat.fista[l-1,] + ((l-2)/(l+1))*(phi.hat.fista[l-1,] - phi.hat.fista[l-2,])
    #     phi.new <- phi.temp + step.size*t(Z%*%(t(Y[ii,]-phi.temp%*%Z)));
    #     phi.new <- soft(phi.new,rep(0,k*m.hat*p),lambda[ll]);
    #     if ( max(abs(phi.new - phi.temp)) < tol) {break;} 
    #     if (max(abs(phi.new - phi.temp)) > tol ) {
    #       phi.hat.fista[l,] <- phi.new; 
    #       # print(l);  
    #     }
    #   }
    #   
    #   # print("l="); print(l)
    #   iter[ii,ll] <- l;
    #   phi.hat.temp[ii,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)] <- phi.new;
    # }
    
    forecast <- matrix(0,k,T.1);
    for(i in 1:m.hat){
      len <- ind.remain[(2*i)] - ind.remain[(2*i-1)];
      delay <- floor(2*len/3);
      delay <- 0;
      l.b <-  ind.remain[(2*i-1)] + delay +  1; u.b <- ind.remain[(2*i)];
      forecast[,(l.b):(u.b)] <- sapply(c((l.b-1):(u.b-1)), 
                                       function(jjj) pred(t(data),phi.hat.temp[,((ll-1)*k*m.hat*p+(i-1)*k*p+1):((ll-1)*k*m.hat*p+(i)*k*p)],p,jjj,k,1)  )
    }
    
    if( ll == 1){
      del.ind <- c();
      for(i in 1:(T.1)){
        if (sum(forecast[,i]^2) == 0){del.ind <- c(del.ind,i)}
      }
    }
    
    # foo <- t(forecast);
    # plot(c(1:300),foo[,1],type="l",col="red")
    # lines(c(1:300),data[,1],col="green")
    
    residual <- t(data) - forecast;
    BIC.temp <- 0;
    for(i in 1:m.hat){
      l.b <-  ind.remain[(2*i-1)] + 0 +  1; u.b <- ind.remain[(2*i)];
      temp <- AIC.BIC(residual[,(l.b):(u.b)],phi.hat.temp[,((ll-1)*k*m.hat*p+(i-1)*k*p+1):((ll-1)*k*m.hat*p+(i)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    
    
    residual <- residual[,-del.ind]; nn <- length(residual[1,]); 
    # check <- AIC.BIC(residual,phi.hat.temp[,((ll-1)*k*m.hat*p+1):(ll*k*m.hat*p)]);
    # pred.error[ll] <- sum(residual^2)/(k*nn); std.res[ll] <- sd(residual[,]^2);
    # pred.error[ll] <- check$BIC;
    pred.error[ll] <- BIC.temp;
    # pred.error[ll] <- check$AIC;
  print(ll)
  }
  ll.final <- which(pred.error==min(pred.error)); ll.final <- min(ll.final);
  # sd.final <- 1*median(std.res)/(sqrt(k*nn));
  # # sd.final <- 1*max(std.res)/(sqrt(k*nn));
  # # sd.final <- 3*sd(pred.error);
  # # sd.final <- quantile(pred.error[1:ll.final],0.975) + pred.error[ll.final];
  # # sd.final <- 1*median(std.res)/(sqrt(nn));
  # p.error <- abs(pred.error[ll.final:length(lambda)] - (pred.error[ll.final] + sd.final));
  # ll.final <- which(p.error==min(p.error)); ll.final <- min(ll.final);
  phi.hat <- phi.hat.temp[,((ll.final-1)*k*m.hat*p+1):(ll.final*k*m.hat*p)];
  
  # for(i in 1:k){
  #   for(j in 1:(k*m.hat*p)){
  #     if(abs(phi.hat[i,j]) < lambda[ll.final]  ){phi.hat[i,j] <- 0;}
  #   }
  # }
  
  return(list(phi.hat = phi.hat, iter = iter, pred.error = pred.error, tune.final = lambda[ll.final]))
  
  
  
}

mspe.plot <- function(method,pred.error,lambda,tune.final,jj){
  plot( lambda, pred.error, type = 'o', col = "blue", main=c(jj,method))
  abline(v=tune.final)
}

mspe.plot.new <- function(method,pred.error,lambda,tune.final,jj){
  plot( seq(1,length(lambda),1), pred.error, type = 'o', col = "blue", main=c(jj,method))
  abline(v=  which(lambda == tune.final)  )
  # abline(v=tune.final)
}

estimation.check <- function(phi,phi.hat){
  k <- length(phi[,1]); p <- (length(phi[1,]))/k; N <- (length(phi.hat[1,]))/(k*p); L <- matrix(0,k,k);
  l2.error <- rep(0,N); true.zero <- rep(0,N); true.non.zero <- rep(0,N); L.hat <- matrix(0,k,k*N); L.hat.final <- rep(0,N);
  false.zero <- rep(0,N); false.non.zero <- rep(0,N); 
  
  count.non <- 0; count.zero <- 0;
  for (i in 1:k){
    for (j in 1:(k*p)){
      if ( phi[i,j] != 0  ){count.non <- count.non + 1;}
      if ( phi[i,j] == 0  ){count.zero <- count.zero + 1;}
    }
  }
  for (i in 1:k){
    for (j in 1:k){
      temp <- 0; 
      for (l in 1:p){if ( phi[i,((l-1)*k+j)] !=0   ) {temp <- temp+1;}}
      L[i,j] <- temp;
    }
  }
  for (jj in 1:N){
    phi.temp <- phi.hat[,((jj-1)*k*p+1):(jj*k*p)];
    count.false.zero <- 0; count.false.non.zero <- 0; count.true.non.zero <- 0; count.true.zero <- 0;
    for (i in 1:k){
      for (j in 1:(k*p)){
        if ( phi[i,j] != 0 && phi.hat[i,((jj-1)*k*p+j)] == 0   ){count.false.zero <- count.false.zero + 1;}
        if ( phi[i,j] == 0 && phi.hat[i,((jj-1)*k*p+j)] != 0   ){count.false.non.zero <- count.false.non.zero + 1;}
        if ( phi[i,j] == 0 && phi.hat[i,((jj-1)*k*p+j)] == 0   ){count.true.zero <- count.true.zero + 1;}
        if ( phi[i,j] != 0 && phi.hat[i,((jj-1)*k*p+j)] != 0   ){count.true.non.zero <- count.true.non.zero + 1;}
        if ( phi[i,j] == 0 ){phi.temp[i,j] <- 0;}
      }
    }
    l2.error[jj] <- sum((phi.temp-phi)^2);
    true.zero[jj] <- count.true.zero; true.non.zero[jj] <- count.true.non.zero;
    false.zero[jj] <- count.false.zero; false.non.zero[jj] <- count.false.non.zero;
  }
  for (jj in 1:N){
    for (i in 1:k){
      for (j in 1:k){
        temp <- 0; 
        for (l in 1:p){if ( phi.hat[i,(((jj-1)*k*p)+((l-1)*k)+j)] !=0   ) {temp <- temp+1;}}
        L.hat[i,(((jj-1)*k)+j)] <- temp;
      }
    }
    
  }
  L.hat.final <- sapply(c(1:N), function(jj) (sum((L.hat[,((jj-1)*k+1):(jj*k)]-L)^2))/(sum(L))  )
  
  return( list(l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
               true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
               false.non.zero.median = median(false.non.zero), lag.error.mean = mean(L.hat.final),
               lag.error.sd = sd(L.hat.final)) )
}

estimation.check.new <- function(phi,phi.final,k,p,m.hat,pts.final.full){
  N <- length(phi.final); 
  l2.error <- rep(0,N); true.zero <- rep(0,N); true.non.zero <- rep(0,N); 
  false.zero <- rep(0,N); false.non.zero <- rep(0,N); 
  
  count.non <- 0; count.zero <- 0;
  for (i in 1:k){
    for (j in 1:(k*m.hat*p)){
      if ( phi[i,j] != 0  ){count.non <- count.non + 1;}
      if ( phi[i,j] == 0  ){count.zero <- count.zero + 1;}
    }
  }

  for (jj in 1:N){
    phi.temp <- phi.final[[jj]]; pt.temp <- pts.final.full[[jj]];
    if( length(phi.temp[1,]) < k*m.hat*p  ){
      phi.temp.new <- matrix(0,k,k*m.hat*p); phi.temp.new[,(1):(k*p)] <- phi.temp[,(1):(k*p)]; ind <- 0;
      for(i in 2:m.hat){
        if( pt.temp[i-1] == 0  ){ind <- ind+1; phi.temp.new[,((i-1)*k*p+1):(i*k*p)] <- phi.temp.new[,((i-2)*k*p+1):((i-1)*k*p)];    }
        if( pt.temp[i-1] != 0  ){phi.temp.new[,((i-1)*k*p+1):(i*k*p)] <- phi.temp[,((i-1-ind)*k*p+1):((i-ind)*k*p)];    }
      }
      phi.temp <- phi.temp.new;
    }
    
    count.false.zero <- 0; count.false.non.zero <- 0; count.true.non.zero <- 0; count.true.zero <- 0;
    for (i in 1:k){
      for (j in 1:(k*m.hat*p)){
        if ( phi[i,j] != 0 && phi.temp[i,j] == 0   ){count.false.zero <- count.false.zero + 1;}
        if ( phi[i,j] == 0 && phi.temp[i,j] != 0   ){count.false.non.zero <- count.false.non.zero + 1;}
        if ( phi[i,j] == 0 && phi.temp[i,j] == 0   ){count.true.zero <- count.true.zero + 1;}
        if ( phi[i,j] != 0 && phi.temp[i,j] != 0   ){count.true.non.zero <- count.true.non.zero + 1;}
        if ( phi[i,j] == 0 ){phi.temp[i,j] <- 0;}
      }
    }
    l2.error[jj] <- sqrt(sum((phi.temp-phi)^2)/sum(phi^2));
    true.zero[jj] <- count.true.zero; true.non.zero[jj] <- count.true.non.zero;
    false.zero[jj] <- count.false.zero; false.non.zero[jj] <- count.false.non.zero;
  }

  
  return( list(l2.error.mean = mean(l2.error), l2.error.sd = sd(l2.error), true.zero.median = median(true.zero),
               true.non.zero.median = median(true.non.zero), false.zero.median = median(false.zero),
               false.non.zero.median = median(false.non.zero)) )
}

AIC.BIC <- function(residual,phi){
  k <- length(phi[,1]); k.lam <- length(phi[1,]); T.new <- length(residual[1,]); count <- 0;
  for (i in 1:k){for (j in 1:k.lam){if(phi[i,j] != 0){count <- count + 1;}}}
  # for ( i in 1:T.new){residual[,i] <- residual[,i] - as.matrix(rowMeans(residual))}
  # sigma.hat <- (1/(T.new-1))*(residual%*%t(residual)) + 10^(-10)*diag(k);
  sigma.hat <- 0*diag(k);
  for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[,i]%*%t(residual[,i]);  }
  sigma.hat <- (1/(T.new))*sigma.hat;
  ee.temp <- min(eigen(sigma.hat)$values); 
  if(ee.temp <= 0){sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp)+10^(-4))*diag(k);}
  # sigma.hat <- (1/(T.new-1))*(t(residual)%*%(residual)) + 0*10^(-10)*diag(T.new);
  log.det <- log(det(sigma.hat)+ 0*10^(-10)); 
  # print(det(sigma.hat));
  return(list(AIC = log.det + 2*count/T.new, BIC = log.det + log(T.new)*count/T.new))
}

AIC.BIC.CV <- function(residual,phi){
  k <- length(phi[,1]); k.lam <- length(phi[1,]); T.new <- length(residual[1,]); count <- 0;
  for (i in 1:k){for (j in 1:k.lam){if(phi[i,j] != 0){count <- count + 1;}}}
  for ( i in 1:T.new){residual[,i] <- residual[,i] - as.matrix(rowMeans(residual))}
  sigma.hat <- (1/(T.new-1))*(residual%*%t(residual)) 
  log.det <- log(det(sigma.hat)); if (abs(det(sigma.hat)) < 10^(-8) ){log.det <- -8^(1)}
  return(list(AIC = log.det + 2*count/T.new, BIC = log.det + log(T.new)*count/T.new))
}


plot.matrix <- function (phi,p,name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  }
  else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p - 
                                                    s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
  at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
  }
  else {
    at2 = c()
  }
  at <- c(at, at2)
  se2 = seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)), col.regions = rgb.palette, 
                  colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = name, 
                                                                         cex = 1), panel = function(...) {
                                                                           panel.levelplot(...)
                                                                           panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + 
                                                                                                                   p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                                                                                                                                                   k + m * s), lwd = 0.5)
                                                                           bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                           b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                           b1 <- c(bl1, b23)
                                                                           panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                           panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                         }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                                   cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0, 
                                                                                                                                              tck = c(0, 0))))
  return(L2)
}



first.step <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4)){
  
  test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  phi.hat.full <- test$phi.hat;
  ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  temp.lam <- ll;
  # temp.lam <- quantile(phi.hat.full,ll)
  ll <- c(0);
  brk.points.list <- vector("list",length(ll));
  
  for(j in 1:length(ll)){
    
    phi.hat <- phi.hat.full;
    # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
    
    n <- T - p;
    m.hat <- 0; brk.points <- rep(0,n);
    
    for (i in 2:n)
    {
      if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
        m.hat <- m.hat + 1; brk.points[m.hat] <- i;
      }
    }
    
    loc <- rep(0,m.hat);
    brk.points <- brk.points[1:m.hat];
    brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
    m.hat <- length(brk.points);
    if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
    loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
    brk.points.list[[j]] <- brk.points;
  }
  
  
  return(list(brk.points = brk.points.list, phi.hat = phi.hat.full))
}


first.step.cv <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4)){
  
  
  kk <- length(lambda);
  cv <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  
  for (i in 1:kk) {
    test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n);
      
      for (iii in 2:n)
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- iii;
        }
      }
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
      m.hat <- length(brk.points);
      if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
      loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    
    
    phi.full.all <- vector("list",T-p);
    phi.temp.cv <- matrix(0,k,k*p);
    forecast <- matrix(0,k,T-p);
    for(j in (p+1):T){
      phi.temp.cv <- phi.temp.cv + phi.hat.full[,((j-p-1)*k*p+1):((j-p)*k*p)];
      phi.full.all[[(j-p)]] <- phi.temp.cv;
      forecast[,(j-p)] <- pred(t(data),phi.temp.cv,p,j-1,k,1)
    }
    
    phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    if( m.hat > 1){
      for(jj in 1:m.hat){
        phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
      }
    }
    phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    
    residual <- t(data[(p+1):T,]) - forecast;
    BIC.temp <- 0;
    if (m.hat == 0){
      temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    if(m.hat >=1){
      temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
      if ( m.hat >= 2){
        for(ii in 2:m.hat){
          l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
          temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
          BIC.temp <- BIC.temp + temp$BIC;
        }
      }
      temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
      BIC.temp <- BIC.temp + temp$BIC;
    }
    
    
    cv[i] <- sum( (forecast - t(data[(p+1):T,])  )^2 );
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  phi.hat.full <- phi.final[[lll]];
  
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}



first.step.cv.new <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k.org <- length(data.temp[1,]);
  data.temp <- data.temp[-cv.index,];
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  
  for (i in 1:kk) {
    if ( i == 1){
      test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi)
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n);
      
      for (iii in 2:n)
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- iii;
        }
      }
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
      m.hat <- length(brk.points);
      if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
      loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    
    
    phi.full.all <- vector("list",T-p);
    phi.temp.cv <- matrix(0,k,k*p);
    forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
    for(j in (p+1):T){
      phi.temp.cv <- phi.temp.cv + phi.hat.full[,((j-p-1)*k*p+1):((j-p)*k*p)];
      phi.full.all[[(j-p)]] <- phi.temp.cv;
      forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
    }
    
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
    }
    
    forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
    residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
      }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    # if( m.hat > 1){
    #   for(jj in 1:m.hat){
    #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
    #   }
    # }
    # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    # 
    # residual <- t(data.temp[(p+1):T,]) - forecast;
    # BIC.temp <- 0;
    # if (m.hat == 0){
    #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    # if(m.hat >=1){
    #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    #   if ( m.hat >= 2){
    #     for(ii in 2:m.hat){
    #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
    #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
    #       BIC.temp <- BIC.temp + temp$BIC;
    #     }
    #   }
    #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    
    
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + cv.var[lll] )   ){ind.new <- ind.new + 1;}
    }
  }
  lll <- lll + ind.new;
  phi.hat.full <- phi.final[[lll]];
  print("CV"); print(cv);
  print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll]))
}



first.step.cv.new.blocks <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index, blocks){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  # if ( cv.l > 0 )  {data.temp <- data.temp[-cv.index,];}
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  flag.full <- rep(0,kk);
  
  for (i in 1:kk) {
    if ( i == 1){
      # test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = NULL, blocks = blocks, cv.index = cv.index)
      flag.full[i] <- test$flag;
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi, blocks = blocks, cv.index = cv.index)
      flag.full[i] <- test$flag;
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # print("TILL HERE!")
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      # print("BREAKPOINTS"); print(brk.points);
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      # if(i == kk){print("INITIAL"); print(brk.points);}
      brk.points <- brk.points[which(brk.points > (2+1)*min(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*min(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(p,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
              
          }
        }
      }

      
      
      
      
      
      
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (1+0)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   # for(i.1 in 2:(length(loc))){
      #   #   if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #   #   # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   # }
      #   
      #   loc <- loc[which(ind==0)]
      #   # brk.points <- brk.points[loc];
      #   
      # }
      
      ######## AGAIN ########################
      # m.hat <- length(brk.points);
      # loc <- rep(0,m.hat);
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) == (1+1)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   for(i.1 in 2:(length(loc))){
      #     if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #     # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   }
      #   
      #   loc <- loc[which(ind==0)]
      #   brk.points <- brk.points[loc];
      #   
      # }

      brk.points.list[[j]] <- brk.points;
    }
    
    
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    # print("TILL HERE!")
    # print(brk.points)
    
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,k,T);
    phi.full.all[[1]] <- phi.hat[,(1):(k*p)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*k*p+1):(i.1*k*p)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data.org),phi.full.all[[i.1-1]],p,blocks[i.1],k,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,k,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],p,blocks[cv.index[j]+1]-1,k,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    residual <- (forecast - t(data.org))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org)])  );
    }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    
    if(FALSE) {
    
    phi.full.all <- vector("list",T-p);
    phi.temp.cv <- matrix(0,k,k*p);
    forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
    for(j in (p+1):T){
      phi.hat.full.temp <- matrix(0,k,k*p);
      ind <- 1;
      for(j.3 in 1:n.new){if( j > blocks[j.3+1]  ) {ind <- ind + 1;}  }
      for(j.4 in 1:ind){phi.hat.full.temp <- phi.hat.full.temp + phi.hat.full[,((j.4-1)*k*p+1):((j.4)*k*p)];}
      
      
      # phi.temp.cv <- phi.temp.cv + phi.hat.full[,((ind-1)*k*p+1):((ind)*k*p)];
      phi.temp.cv <- phi.hat.full.temp;
      phi.full.all[[(j-p)]] <- phi.temp.cv;
      forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
    }
    
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
    }
    
    forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
    residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
    }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
    # if( m.hat > 1){
    #   for(jj in 1:m.hat){
    #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
    #   }
    # }
    # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
    # 
    # residual <- t(data.temp[(p+1):T,]) - forecast;
    # BIC.temp <- 0;
    # if (m.hat == 0){
    #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    # if(m.hat >=1){
    #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    #   if ( m.hat >= 2){
    #     for(ii in 2:m.hat){
    #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
    #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
    #       BIC.temp <- BIC.temp + temp$BIC;
    #     }
    #   }
    #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
    #   BIC.temp <- BIC.temp + temp$BIC;
    # }
    
    
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
    # cv[i] <- BIC.temp;
    # temp <- AIC.BIC.CV(residual,phi.hat.temp);
    # cv[i] <- temp$BIC;
    }
    # print("====================================")
  }
  
  # lll <- min(which(cv[which(flag.full==0)] == min(cv[which(flag.full==0)])));
  lll <- min(which(cv == min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + 1*sqrt(n.new/T)*cv.var[lll] ) && flag.full[i.3] != 1   ){ind.new <- ind.new + 1;}
    }
  }
  # lll <- lll + ind.new;
  lll <- lll + 0;
  # lll <- min(lll+1,kk);
  phi.hat.full <- phi.final[[lll]];
  # print("CV"); print(cv);
  # print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll], phi.hat = phi.final[[lll]]))
}

first.step.cv.new.blocks.rank <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index, blocks, initial=NULL){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  # if ( cv.l > 0 )  {data.temp <- data.temp[-cv.index,];}
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  flag.full <- rep(0,kk);
  
  for (i in 1:kk) {
    if ( i == 1){
      if(  is.null(initial) ){
        # test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
        test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = NULL, blocks = blocks, cv.index = cv.index)
        flag.full[i] <- test$flag;
      }
      if(  !is.null(initial) ){
        # test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
        test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial, blocks = blocks, cv.index = cv.index)
        flag.full[i] <- test$flag;
      }
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi, blocks = blocks, cv.index = cv.index)
      flag.full[i] <- test$flag;
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # print("TILL HERE!")
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      # print("BREAKPOINTS"); print(brk.points);
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      # if(i == kk){print("INITIAL"); print(brk.points);}
      brk.points <- brk.points[which(brk.points > (2+1)*min(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*min(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(p,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
            
          }
        }
      }
      
      
      
      
      
      
      
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (1+0)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   # for(i.1 in 2:(length(loc))){
      #   #   if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #   #   # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   # }
      #   
      #   loc <- loc[which(ind==0)]
      #   # brk.points <- brk.points[loc];
      #   
      # }
      
      ######## AGAIN ########################
      # m.hat <- length(brk.points);
      # loc <- rep(0,m.hat);
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) == (1+1)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   for(i.1 in 2:(length(loc))){
      #     if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #     # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   }
      #   
      #   loc <- loc[which(ind==0)]
      #   brk.points <- brk.points[loc];
      #   
      # }
      
      brk.points.list[[j]] <- brk.points;
    }
    
    
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    # print("TILL HERE!")
    # print(brk.points)
    
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,k,T);
    phi.full.all[[1]] <- phi.hat[,(1):(k*p)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*k*p+1):(i.1*k*p)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data.org),phi.full.all[[i.1-1]],p,blocks[i.1],k,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,k,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],p,blocks[cv.index[j]+1]-1,k,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    residual <- (forecast - t(data.org))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org)])  );
    }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    
    if(FALSE) {
      
      phi.full.all <- vector("list",T-p);
      phi.temp.cv <- matrix(0,k,k*p);
      forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
      for(j in (p+1):T){
        phi.hat.full.temp <- matrix(0,k,k*p);
        ind <- 1;
        for(j.3 in 1:n.new){if( j > blocks[j.3+1]  ) {ind <- ind + 1;}  }
        for(j.4 in 1:ind){phi.hat.full.temp <- phi.hat.full.temp + phi.hat.full[,((j.4-1)*k*p+1):((j.4)*k*p)];}
        
        
        # phi.temp.cv <- phi.temp.cv + phi.hat.full[,((ind-1)*k*p+1):((ind)*k*p)];
        phi.temp.cv <- phi.hat.full.temp;
        phi.full.all[[(j-p)]] <- phi.temp.cv;
        forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
      }
      
      for(j in (1):cv.l){
        forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
      }
      
      forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
      residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
      var.matrix <-  matrix(0,k,m.hat+1);
      if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
      if ( m.hat >= 1){
        var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
        var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
      }
      if ( m.hat >=2 ){
        for(mm in 2:m.hat){
          var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
        }
      }
      
      
      if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
      if ( m.hat >= 1){
        for (i.1 in 1:cv.l){
          ind <- 0;
          for (i.2 in 1:m.hat){
            if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
          }
          cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
        }
      }
      
      cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
      
      # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
      # if( m.hat > 1){
      #   for(jj in 1:m.hat){
      #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
      #   }
      # }
      # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
      # 
      # residual <- t(data.temp[(p+1):T,]) - forecast;
      # BIC.temp <- 0;
      # if (m.hat == 0){
      #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      # if(m.hat >=1){
      #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      #   if ( m.hat >= 2){
      #     for(ii in 2:m.hat){
      #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
      #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
      #       BIC.temp <- BIC.temp + temp$BIC;
      #     }
      #   }
      #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      
      
      cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
      # cv[i] <- BIC.temp;
      # temp <- AIC.BIC.CV(residual,phi.hat.temp);
      # cv[i] <- temp$BIC;
    }
    # print("====================================")
  }
  
  # lll <- min(which(cv[which(flag.full==0)] == min(cv[which(flag.full==0)])));
  lll <- min(which(cv == min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + 1*sqrt(n.new/T)*cv.var[lll] ) && flag.full[i.3] != 1   ){ind.new <- ind.new + 1;}
    }
  }
  lll <- lll + ind.new;
  # lll <- lll + 0;
  # lll <- min(lll+1,kk);
  phi.hat.full <- phi.final[[lll]];
  # print("CV"); print(cv);
  # print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll], phi.hat = phi.final[[lll]]))
}

first.step.cv.new.blocks.isb <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index, blocks){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k.org <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  # if ( cv.l > 0 )  {data.temp <- data.temp[-cv.index,];}
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  flag.full <- rep(0,kk);
  
  
  for (i in 1:kk) {
    if ( i == 1){
      # test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = NULL, blocks = blocks, cv.index = cv.index)
      flag.full[i] <- test$flag;
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      if (flag.full[i-1] == 1){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi, blocks = blocks, cv.index = cv.index)
      flag.full[i] <- test$flag;
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # print("TILL HERE!")
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      # print("BREAKPOINTS"); print(brk.points);
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      # if(i == kk){print("INITIAL"); print(brk.points);}
      brk.points <- brk.points[which(brk.points > (2+1)*min(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*min(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(p,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
            
          }
        }
      }
      
      
      
      
      
      
      
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (1+0)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   # for(i.1 in 2:(length(loc))){
      #   #   if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #   #   # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   # }
      #   
      #   loc <- loc[which(ind==0)]
      #   # brk.points <- brk.points[loc];
      #   
      # }
      
      ######## AGAIN ########################
      # m.hat <- length(brk.points);
      # loc <- rep(0,m.hat);
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) == (1+1)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   for(i.1 in 2:(length(loc))){
      #     if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #     # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   }
      #   
      #   loc <- loc[which(ind==0)]
      #   brk.points <- brk.points[loc];
      #   
      # }
      
      brk.points.list[[j]] <- brk.points;
    }
    
    
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    # print("TILL HERE!")
    # print(brk.points)
    
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,k,T);
    phi.full.all[[1]] <- phi.hat[,(1):(k*p)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*k*p+1):(i.1*k*p)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data),phi.full.all[[i.1-1]],p,blocks[i.1],k,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,k,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],p,blocks[cv.index[j]+1]-1,k,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    residual <- (forecast - t(data.org))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org)])  );
    }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    
    if(FALSE) {
      
      phi.full.all <- vector("list",T-p);
      phi.temp.cv <- matrix(0,k,k*p);
      forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
      for(j in (p+1):T){
        phi.hat.full.temp <- matrix(0,k,k*p);
        ind <- 1;
        for(j.3 in 1:n.new){if( j > blocks[j.3+1]  ) {ind <- ind + 1;}  }
        for(j.4 in 1:ind){phi.hat.full.temp <- phi.hat.full.temp + phi.hat.full[,((j.4-1)*k*p+1):((j.4)*k*p)];}
        
        
        # phi.temp.cv <- phi.temp.cv + phi.hat.full[,((ind-1)*k*p+1):((ind)*k*p)];
        phi.temp.cv <- phi.hat.full.temp;
        phi.full.all[[(j-p)]] <- phi.temp.cv;
        forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
      }
      
      for(j in (1):cv.l){
        forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
      }
      
      forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
      residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
      var.matrix <-  matrix(0,k,m.hat+1);
      if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
      if ( m.hat >= 1){
        var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
        var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
      }
      if ( m.hat >=2 ){
        for(mm in 2:m.hat){
          var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
        }
      }
      
      
      if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
      if ( m.hat >= 1){
        for (i.1 in 1:cv.l){
          ind <- 0;
          for (i.2 in 1:m.hat){
            if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
          }
          cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
        }
      }
      
      cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
      
      # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
      # if( m.hat > 1){
      #   for(jj in 1:m.hat){
      #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
      #   }
      # }
      # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
      # 
      # residual <- t(data.temp[(p+1):T,]) - forecast;
      # BIC.temp <- 0;
      # if (m.hat == 0){
      #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      # if(m.hat >=1){
      #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      #   if ( m.hat >= 2){
      #     for(ii in 2:m.hat){
      #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
      #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
      #       BIC.temp <- BIC.temp + temp$BIC;
      #     }
      #   }
      #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      
      
      cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
      # cv[i] <- BIC.temp;
      # temp <- AIC.BIC.CV(residual,phi.hat.temp);
      # cv[i] <- temp$BIC;
    }
    # print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + 1*sqrt(n.new/T)*cv.var[lll] )   ){ind.new <- ind.new + 1;}
    }
  }
  lll <- lll + ind.new;
  # lll <- min(lll+1,kk);
  phi.hat.full <- phi.final[[lll]];
  # print("CV"); print(cv);
  # print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll], phi.hat = phi.final[[lll]]))
}

first.step.cv.new.blocks.data <- function(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4),cv.index, blocks){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); k.org <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  # if ( cv.l > 0 )  {data.temp <- data.temp[-cv.index,];}
  kk <- length(lambda);
  cv <- rep(0,kk); cv.var <- rep(0,kk);
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); k <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  
  for (i in 1:kk) {
    if ( i == 1){
      # test <- var.break.fit(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol)
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = NULL, blocks = blocks, cv.index = cv.index)
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[(i-1)]]
      test <- var.break.fit.block(method,data.temp, weight = NULL, lambda[i], p, max.iteration = max.iteration, tol = tol, initial.phi = initial.phi, blocks = blocks, cv.index = cv.index)
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    # print("TILL HERE!")
    # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
    # temp.lam <- ll;
    # temp.lam <- quantile(phi.hat.full,ll)
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      
      phi.hat <- phi.hat.full;
      # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
      
      n <- T - p;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*k*p+1):(iii*k*p)] )^2 ) != 0   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      # print("BREAKPOINTS"); print(brk.points);
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      # if(i == kk){print("INITIAL"); print(brk.points);}
      brk.points <- brk.points[which(brk.points > (2+1)*min(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*min(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(p,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
            
          }
        }
      }
      
      
      
      
      
      
      
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (1+0)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   # for(i.1 in 2:(length(loc))){
      #   #   if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #   #   # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   # }
      #   
      #   loc <- loc[which(ind==0)]
      #   # brk.points <- brk.points[loc];
      #   
      # }
      
      ######## AGAIN ########################
      # m.hat <- length(brk.points);
      # loc <- rep(0,m.hat);
      # # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)*min(blocks.size)  ) {loc[mm] <- mm;}  }}
      # # loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
      # if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) == (1+1)*max(p,min(blocks.size))  ) {loc[mm] <- mm;}  }}
      # loc <- loc[which(loc!=0)]; 
      # if( length(loc) > 1){
      #   ind <- rep(0,length(loc));
      #   for(i.1 in 2:length(loc)){
      #     if( (loc[i.1]-loc[i.1-1]) <= 1){ind[i.1] <- 1;}
      #   }
      #   for(i.1 in 2:(length(loc))){
      #     if( ind[i.1] == 1 && ind[i.1-1] == 1 ){ind[i.1] <- 0;}
      #     # if( ind[i.1] == 1 && ind[i.1+1] == 1 ){ind[i.1] <- 0;}
      #   }
      #   
      #   loc <- loc[which(ind==0)]
      #   brk.points <- brk.points[loc];
      #   
      # }
      
      brk.points.list[[j]] <- brk.points;
    }
    
    
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    # print("TILL HERE!")
    # print(brk.points)
    
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,k,T);
    phi.full.all[[1]] <- phi.hat[,(1):(k*p)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*k*p+1):(i.1*k*p)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data),phi.full.all[[i.1-1]],p,blocks[i.1],k,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,k,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],p,blocks[cv.index[j]+1]-1,k,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    residual <- (forecast - t(data.org))^2;
    var.matrix <-  matrix(0,k,m.hat+1);
    if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
    if ( m.hat >= 1){
      var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
      var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org)])  );
    }
    if ( m.hat >=2 ){
      for(mm in 2:m.hat){
        var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
      }
    }
    
    
    if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
    if ( m.hat >= 1){
      for (i.1 in 1:cv.l){
        ind <- 0;
        for (i.2 in 1:m.hat){
          if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
        }
        cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
      }
    }
    
    cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
    
    
    if(FALSE) {
      
      phi.full.all <- vector("list",T-p);
      phi.temp.cv <- matrix(0,k,k*p);
      forecast <- matrix(0,k,T-p); forecast.new <- matrix(0,k,cv.l);
      for(j in (p+1):T){
        phi.hat.full.temp <- matrix(0,k,k*p);
        ind <- 1;
        for(j.3 in 1:n.new){if( j > blocks[j.3+1]  ) {ind <- ind + 1;}  }
        for(j.4 in 1:ind){phi.hat.full.temp <- phi.hat.full.temp + phi.hat.full[,((j.4-1)*k*p+1):((j.4)*k*p)];}
        
        
        # phi.temp.cv <- phi.temp.cv + phi.hat.full[,((ind-1)*k*p+1):((ind)*k*p)];
        phi.temp.cv <- phi.hat.full.temp;
        phi.full.all[[(j-p)]] <- phi.temp.cv;
        forecast[,(j-p)] <- pred(t(data.temp),phi.temp.cv,p,j-1,k,1)
      }
      
      for(j in (1):cv.l){
        forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j]-1-p-j+1)]],p,cv.index[j]-1,k,1)
      }
      
      forecast.all <- matrix(0,k,T.org-p); forecast.all[,cv.index] <- forecast.new; forecast.all[,-cv.index] <- forecast;
      residual <- (forecast.all - t(data.org[((p+1):T.org),]))^2;
      var.matrix <-  matrix(0,k,m.hat+1);
      if ( m.hat == 0){var.matrix <- sapply(c(1:k), function(jjj)  var(residual[jjj,])  );}
      if ( m.hat >= 1){
        var.matrix[,1] <- t(as.vector(sapply(c(1:k), function(jjj)  var(residual[jjj,(1):(brk.points[1]-1)])  )));
        var.matrix[,(m.hat+1)] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(m.hat)]):(T.org-p)])  );
      }
      if ( m.hat >=2 ){
        for(mm in 2:m.hat){
          var.matrix[,mm] <- sapply(c(1:k), function(jjj)  var(residual[jjj,(brk.points[(mm-1)]):(brk.points[mm]-1)])  );
        }
      }
      
      
      if ( m.hat == 0) {cv.var[i] <- sum(var.matrix);}
      if ( m.hat >= 1){
        for (i.1 in 1:cv.l){
          ind <- 0;
          for (i.2 in 1:m.hat){
            if (cv.index[i.1] > brk.points[i.2]){ind <- ind + 1;}
          }
          cv.var[i] <- cv.var[i] + sum(var.matrix[,(ind+1)]);
        }
      }
      
      cv.var[i] <- (1/(k*cv.l))*sqrt(cv.var[i]);
      
      # phi.hat.temp <- matrix(0,k,(k*p*(m.hat+1)));
      # if( m.hat > 1){
      #   for(jj in 1:m.hat){
      #     phi.hat.temp[,((jj-1)*k*p+1):((jj)*k*p)] <- phi.full.all[[(brk.points[jj]-p-1)]];
      #   }
      # }
      # phi.hat.temp[,((m.hat)*k*p+1):((m.hat+1)*k*p)] <- phi.full.all[[(T-p)]];
      # 
      # residual <- t(data.temp[(p+1):T,]) - forecast;
      # BIC.temp <- 0;
      # if (m.hat == 0){
      #   temp <- AIC.BIC.CV(residual[,(1):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      # if(m.hat >=1){
      #   temp <- AIC.BIC.CV(residual[,(1):(brk.points[1]-1)],phi.hat.temp[,((1-1)*k*m.hat*p+(1-1)*k*p+1):((1-1)*k*m.hat*p+(1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      #   if ( m.hat >= 2){
      #     for(ii in 2:m.hat){
      #       l.b <-  brk.points[(ii-1)]; u.b <- brk.points[ii]-1;
      #       temp <- AIC.BIC.CV(residual[,(l.b):(u.b)],phi.hat.temp[,((1-1)*k*m.hat*p+(ii-1)*k*p+1):((1-1)*k*m.hat*p+(ii)*k*p)]);
      #       BIC.temp <- BIC.temp + temp$BIC;
      #     }
      #   }
      #   temp <- AIC.BIC.CV(residual[,(brk.points[m.hat]):(T-p)],phi.hat.temp[,((1-1)*k*m.hat*p+(m.hat)*k*p+1):((1-1)*k*m.hat*p+(m.hat+1)*k*p)]);
      #   BIC.temp <- BIC.temp + temp$BIC;
      # }
      
      
      cv[i] <- (1/(k*cv.l))*sum( (forecast.new - t(data.org[cv.index,])  )^2 );
      # cv[i] <- BIC.temp;
      # temp <- AIC.BIC.CV(residual,phi.hat.temp);
      # cv[i] <- temp$BIC;
    }
    print("====================================")
  }
  
  lll <- min(which(cv==min(cv)));
  ind.new <- 0;
  if (lll < kk){
    for(i.3 in (lll+1):(kk)){
      if ( cv[i.3] < (cv[lll] + 1*sqrt(n.new/T)*cv.var[lll] )   ){ind.new <- ind.new + 1;}
    }
  }
  ind.new <- 0;
  lll <- lll + ind.new;
  # lll <- min(lll+1,kk);
  phi.hat.full <- phi.final[[lll]];
  # print("CV"); print(cv);
  # print("CV.VAR"); print(cv.var);
  
  # test <- var.break.fit(method,data.temp, weight = NULL, lambda, p, max.iteration = max.iteration, tol = tol)
  # phi.hat.full <- test$phi.hat;
  # ll <- seq(0.05,max(abs(phi.hat.full)),0.05);
  # temp.lam <- ll;
  # # temp.lam <- quantile(phi.hat.full,ll)
  # ll <- c(0);
  # brk.points.list <- vector("list",length(ll));
  # 
  # for(j in 1:length(ll)){
  #   
  #   phi.hat <- phi.hat.full;
  #   # phi.hat <- soft.full(phi.hat.full,temp.lam[j],1,1,1);
  #   
  #   n <- T - p;
  #   m.hat <- 0; brk.points <- rep(0,n);
  #   
  #   for (i in 2:n)
  #   {
  #     if ( sum((phi.hat[,((i-1)*k*p+1):(i*k*p)] )^2 ) != 0   ){
  #       m.hat <- m.hat + 1; brk.points[m.hat] <- i;
  #     }
  #   }
  #   
  #   loc <- rep(0,m.hat);
  #   brk.points <- brk.points[1:m.hat];
  #   brk.points <- brk.points[which(brk.points > (p+3))]; brk.points <- brk.points[which(brk.points < (n))];
  #   m.hat <- length(brk.points);
  #   if (m.hat > 1){for(mm in 1:(m.hat-1)){ if( abs(brk.points[mm+1] - brk.points[mm]) <= (p+1)  ) {loc[mm] <- mm;}  }}
  #   loc <- loc[which(loc!=0)]; if( length(loc) > 0){brk.points <- brk.points[-loc];}
  #   brk.points.list[[j]] <- brk.points;
  # }
  
  
  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv.final = lambda[lll], phi.hat = phi.final[[lll]]))
}

generate.phi.new <- function(k,p.t,l.bound,max.iter,seed){
  
  for ( try in 1:max.iter){
    phi <- matrix(0,k,k*p.t)
    set.seed(try+seed)
    base <- matrix(2*(runif((k^2)*p.t)-1/2),k,k*p.t);
    for(i in 1:k){
      for(j in 1:(k*p.t)){
        d <- (abs(abs(i-j)-1)+1);
        if ( abs(base[i,j]/d) > l.bound   ){phi[i,j] <- base[i,j]/d; }
        if ( abs(base[i,j]/d) <= l.bound   ){phi[i,j] <- 0; }
      }
    }
    
    companion.phi <- matrix(0,k*p.t,k*p.t);
    companion.phi[1:k,] <- phi;
    if( p.t > 1){
      for(i in 1:(p.t-1)){companion.phi[(i*k+1):((i+1)*k),((i-1)*k+1):((i)*k)] <- diag(k)       }
    }
    aaa <- eigen(companion.phi)$values; aaa <- Mod(aaa);
    # print("TRY=="); print(try)
    if(max(aaa) < 1){break;}
  }
  
  return(list(phi = phi, try = try))
}

detection.refinement <- function(data,omega.2,break.pts,phi.hat,p){
  T <- length(data[,1]); k <- length(data[1,]); l <- length(break.pts);
  if(l >= 1){
    for(i in 1:l){
      lb <- break.pts[i] - omega.2; ub <- break.pts[i] + omega.2;
      phi.1 <- phi.hat[,((i-1)*k*p+1):((i)*k*p)];
      phi.2 <- phi.hat[,((i)*k*p+1):((i+1)*k*p)];
      break.pts[i] <- loc.finder(data, lb, ub, phi.1, phi.2, p);
  
    }
    
  }
  return(break.pts)
}

detection.refinement.local <- function(data,omega.2,break.pts,phi.hat,p,phi.hat.actual=NULL){
  T <- length(data[,1]); k <- length(data[1,]); l <- length(break.pts);
  if(l >= 1){
    for(i in 1:l){
      # print(i)
      lb <- break.pts[i] - omega.2; ub <- break.pts[i] + omega.2;
      phi.1 <- phi.hat[[(2*i-1)]];
      phi.2 <- phi.hat[[(2*i)]];
      # phi.1 <- phi.hat.actual[,((i-1)*k*p):((i)*k*p)];
      # phi.2 <- phi.hat.actual[,((i)*k*p):((i+1)*k*p)];
      break.pts[i] <- loc.finder.par(data, lb, ub, phi.1, phi.2, p);
      
    }
    
  }
  return(break.pts)
}


loc.finder <- function(data, lb, ub, phi.1, phi.2, p){
  T <- length(data[,1]); k <- length(data[1,]);
  temp <- rep(0,ub-lb+1);
  for(i in lb:ub){
    
    forecast <- matrix(0,k,ub-lb+1);
    forecast <- pred.block.new.local(t(data),phi.1,phi.2,p,lb-1,k,ub-lb+1,i)
    
    
    temp[i-lb+1] <- sum((forecast -  t(data[lb:ub,]))^2)
    
  }
  lll <- min(which(temp==min(temp)));
  # print(lll);
  return(lll+lb-1)
}

loc.finder.par <- function(data, lb, ub, phi.1, phi.2, p){
  T <- length(data[,1]); k <- length(data[1,]);
  # temp <- rep(0,ub-lb+1);
  h <- ub-lb+1; sep <- i; T.1 <- lb-1;
  forecast.1 <- foreach(j=1:h, .combine=cbind) %dopar% {
    temp.2 <- matrix(0,k,1);
      for (i.1 in 1:p){temp.2 <- temp.2 +  phi.1[,((i.1-1)*k+1):(i.1*k)]%*%t(data)[,j+T.1-i.1];}
    temp.2
  }
  forecast.2 <- foreach(j=1:h, .combine=cbind) %dopar% {
    temp.2 <- matrix(0,k,1);
    for (i.1 in 1:p){temp.2 <- temp.2 +  phi.2[,((i.1-1)*k+1):(i.1*k)]%*%t(data)[,j+T.1-i.1];}
    temp.2
  }
  temp <- foreach(i=lb:ub, .combine=cbind) %dopar% {
    # print("UP to here!")
    forecast <- matrix(0,k,ub-lb+1);
    if(i==lb){forecast <- forecast.1;}
    if(i==ub){forecast <- forecast.2;}
    if(i>lb && i<ub){
      forecast[,(1):(i-lb+1)] <- forecast.1[,(1):(i-lb+1)];
      forecast[,(i-lb+1+1):(ub-lb+1)] <- forecast.2[,(i-lb+1+1):(ub-lb+1)];
    }
    sum((forecast -  t(data[lb:ub,]))^2)
  }
  
  temp <- as.numeric(temp)
  # print(dim(temp))
  # print(temp)
  
  # temp <- foreach(i=lb:ub, .combine=cbind) %dopar% {
  #   # print("UP to here!")
  #   forecast <- matrix(0,k,ub-lb+1);
  #   h <- ub-lb+1; sep <- i; T.1 <- lb-1;
  #   
  #   library(doParallel)
  #   cl <- makeCluster(4)
  #   registerDoParallel(cl)
  #   getDoParWorkers()
  #   
  #   temp.1 <- foreach(j=1:h, .combine=cbind) %dopar% {
  #     temp.2 <- matrix(0,k,1);
  #     if( j + T.1 <= sep ){
  #       for (i.1 in 1:p){temp.2 <- temp.2 +  phi.1[,((i.1-1)*k+1):(i.1*k)]%*%t(data)[,j+T.1-i.1];}
  #     }
  #     if( j + T.1 > sep ){
  #       for (i.1 in 1:p){temp.2 <- temp.2 +  phi.2[,((i.1-1)*k+1):(i.1*k)]%*%t(data)[,j+T.1-i.1];}
  #     }
  #     temp.2
  #   }
  #   forecast <- temp.1;
  #   
  #   # forecast <- pred.block.new.local(t(data),phi.1,phi.2,p,lb-1,k,ub-lb+1,i);
  #   sum((forecast -  t(data[lb:ub,]))^2)
  # }

  lll <- min(which(temp==min(temp)));
  # print(lll);
  return(lll+lb-1)
}

detection.check <- function(pts.final.block, brk, N){
  m <- length(brk); len <- rep(0,N);
  for(i in 1:N){len[i] <- length(pts.final.block[[i]]);}
  freq <- as.matrix(table(len)/N);
  pts.final.full <- vector("list",N);
  for(i in 1:N){
    if ( length(pts.final.block[[i]]) > (m-1)   ){pts.final.block[[i]] <- remove.extra.pts(pts.final.block[[i]], brk);}
    if ( length(pts.final.block[[i]]) == (m-1)   ){pts.final.full[[i]] <- pts.final.block[[i]];}
    if ( length(pts.final.block[[i]]) == 0 ) {  pts.final.full[[i]] <- rep(0,m-1);  }
    if ( length(pts.final.block[[i]]) > 0 && length(pts.final.block[[i]]) < (m-1) ){
      ll <- length(pts.final.block[[i]]); pts.final.full[[i]] <- rep(0,m-1);
      for(j in 1:ll){
        if (  pts.final.block[[i]][j] < (brk[1] + (1/2)*(brk[2] - brk[1]) )   ) {pts.final.full[[i]][1] <- pts.final.block[[i]][j];}
        for(kk in 2:(m-1)){
          if (  pts.final.block[[i]][j] >=  (brk[(kk-1)] + (1/2)*(brk[kk] - brk[(kk-1)]) ) && pts.final.block[[i]][j] <  (brk[(kk)] + (1/2)*(brk[kk+1] - brk[(kk)]) )   ){
            pts.final.full[[i]][kk] <- pts.final.block[[i]][j];
          }
        }
      }
    }
  }
  
  detection <- matrix(0,m,5)
  
  detection[1,1] <- c("break points"); detection[1,2] <- c("truth"); detection[1,3] <- c("mean");
  detection[1,4] <- c("std"); detection[1,5] <- c("selection rate");
  for(i in 1:(m-1)){
    detection[(i+1),1] <- c(i); detection[(i+1),2] <- c(brk[i]/T); loc <- rep(0,N);
    for(j in 1:N){
      temp <- pts.final.full[[j]]; l <- length(temp); loc[j] <- temp[i];
    }
    loc <- loc[which(loc!=0)]; T.new <- length(loc); detection[(i+1),3] <- mean(loc/T); detection[(i+1),4] <- sd(loc/T); detection[(i+1),5] <- T.new/N;
  }
  
  for(i in 2:(m)){
    for(j in 2:5){
      detection[i,j] <- round(as.numeric(detection[i,j]),digits = 4)
    }
  }
  return(list(pts.final.full = pts.final.full, detection = detection, freq = freq))
}

remove.extra.pts <- function(pts, brk){
  m.hat <- length(brk)-1;
  if(length(pts) <= m.hat){break;}
  pts.temp <- rep(0,m.hat);
  for(i in 1:m.hat){
    origin <- brk[i];
    dis <- rep(0,length(pts));
    for(j in 1:length(pts)){
      dis[j] <- abs(origin - pts[j]);
    }
    ll <- min(which.min(dis));
    pts.temp[i] <- pts[ll];
  }
  
  
  
  pts <- pts.temp;
  return(pts)
}



generate.phi <- function(k,p.t,l.bound,max.iter,seed){
  
  for ( try in 1:max.iter){
    phi <- matrix(0,k,k*p.t)
    set.seed(try+seed)
    base <- matrix(2*(runif(k^2)-1/2),k,k);
    for ( i in 1:p.t){
      for( j in 1:k){
        for ( l in 1:k){
          if ( abs(base[j,l]/(i*(1+abs(j-l)))) > l.bound  ){phi[j,((i-1)*k+l)] <- base[j,l]/(i*(1+abs(j-l)));}
          else {phi[j,((i-1)*k+l)] <- 0;}
          
        }
      }
    }
    companion.phi <- matrix(0,k*p.t,k*p.t);
    companion.phi[1:k,] <- phi;
    if( p.t > 1){
      for(i in 1:(p.t-1)){companion.phi[(i*k+1):((i+1)*k),((i-1)*k+1):((i)*k)] <- diag(k)       }
    }
    aaa <- eigen(companion.phi)$values; aaa <- Mod(aaa);
    # print("TRY=="); print(try)
    if(max(aaa) < 1){break;}
  }
  
  return(list(phi = phi, try = try))
}

generate.phi.new <- function(k,p.t,l.bound,max.iter,seed){
  
  for ( try in 1:max.iter){
    phi <- matrix(0,k,k*p.t)
    set.seed(try+seed)
    base <- matrix(2*(runif((k^2)*p.t)-1/2),k,k*p.t);
    for(i in 1:k){
      for(j in 1:(k*p.t)){
        d <- (abs(abs(i-j)-1)+1);
        if ( abs(base[i,j]/d) > l.bound   ){phi[i,j] <- base[i,j]/d; }
        if ( abs(base[i,j]/d) <= l.bound   ){phi[i,j] <- 0; }
      }
    }
    
    companion.phi <- matrix(0,k*p.t,k*p.t);
    companion.phi[1:k,] <- phi;
    if( p.t > 1){
      for(i in 1:(p.t-1)){companion.phi[(i*k+1):((i+1)*k),((i-1)*k+1):((i)*k)] <- diag(k)       }
    }
    aaa <- eigen(companion.phi)$values; aaa <- Mod(aaa);
    # print("TRY=="); print(try)
    if(max(aaa) < 1){break;}
  }
  
  return(list(phi = phi, try = try))
}


block.finder <- function(pts,omega.local){
  nn <- length(pts);
   
  if( nn == 1){b <- pts;}
  if( nn > 1){
    b <- vector("list",nn);
    i.ind <- 1;
    jj <- 0;
    while (i.ind < nn) {
      ct <- 1;
      jj <- jj + 1;
      for (j in (i.ind+1):nn) {
        if( abs(pts[i.ind] - pts[j]  ) <= omega.local   ){ct <- ct + 1;}
      }
      b[[jj]] <- pts[(i.ind):(i.ind+ct-1)];
      i.ind <- i.ind + ct;
      # print(b[[jj]])
    }
    l <- length(b[[jj]]);
    if(b[[jj]][l] != pts[nn]  ){
      jj <- jj + 1;
      b[[(jj)]] <- c(pts[nn])   
    }
    b <- b[(1):(jj)];
  }
  
  return(b = b)
}


second.step.final <- function(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local){
  N <- length(data[,1]); k <- length(data[1,]);
  n <- length(pts.list);
  final.pts <- rep(0,n);
  for(i in 1:n){
    pts.temp <- pts.list[[i]];
    if( length(pts.temp) <= 1  ) {final.pts[i] <- pts.temp;}
    if( length(pts.temp) > 1  ){
      
      lb <- max(1,min(pts.temp) - 0.5*omega.local);
      ub <- min(N-1,max(pts.temp) + 0.5*omega.local);
      lam.temp <- (1/1)*(log(ub-lb)*log(k))/(ub-lb);
      temp <- second.step.middle(data[(lb):(ub),], lambda = lam.temp, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.temp-lb+1)
      final.pts[i] <- temp$pts + lb - 1;
    }
    
  }
  
 return(pts = final.pts) 
}

second.step.final.cluster <- function(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local, al){
  N <- length(data[,1]); k <- length(data[1,]);     
  n <- length(pts.list);
  final.pts <- c();
  for(i in 1:n){
    pts.temp <- pts.list[[i]];
    
    if( length(pts.temp) <= 1  ) {final.pts <- c(final.pts,pts.temp);}
    if( length(pts.temp) > 1  ){
      
      lb <- floor(max(1,min(pts.temp) - 0.95*omega.local))+1;
      ub <- floor(min(N-1,max(pts.temp) + 0.95*omega.local));
      lam.temp <- (1/1)*(log(ub-lb)*log(k))/(ub-lb);
      omega.local <- floor(0.95*min(abs(diff(c(p+1,pts.temp-lb+1)))));
      # temp <- second.step.middle(data[(lb):(ub),], lambda = lam.temp, p, max.iteration = 1000, tol = 10^(-4), step.size = 10^(-3), pts.temp-lb+1)
      temp <- second.step.local.new(data[(lb):(ub),], lambda = lam.temp, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.temp-lb+1, 1, omega.local, al)
      final.pts <- c(final.pts,temp$pts + lb - 1);
    }
    
  }
  
  return(pts = final.pts) 
}

##############################################################################
##############################################################################
##############################################################################


bss.block.select <- function(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), al){
  T <- length(data[,1]); k <- length(data[1,]);
  # print(k);
  
  m1 <- 0;
  while (m1 < 14){
    m1 <- m1 + 1;
    
    if( T^(0.25 + (m1-1)*0.05)/5 > 2   ){break;}
    
  }
  
  block.full <- vector("list",11-m1+1);
  NN <- length(block.full);
  pts.final <- vector("list",NN);
  
  for (i in 1:NN ){
    block.full[[i]] <- seq(0,T,floor( T/( floor(T^(0.25+0.05*(m1+i-2))    ) ) ) );
  }
  
  ind <- NN;
  mm <- 0;
  while(mm < NN){
    mm <- mm + 1
    print(lambda.1.cv)
    temp <- bss(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), block.full[[mm]], al)
    if ( length(temp$pts.3) > 0 ){
        pts.final[[mm]] <- temp$pts.3
    }
    omega.local <- temp$omega.local
    # if ( length(temp$pts.3) == 0  ){pts.final[[mm]] <- c();}
    
    
    # lambda.1.cv <- lambda.1.cv * 0.9
    if ( mm > 2 && length(pts.final[[mm]]) == length(pts.final[[mm-1]])    ){      ### mm > 3 for simulations E2
      ind <- mm;
      print("mmmm"); print(mm);
      break;
    }
    
    # if(mm > 3){
    #     ind <- mm
    #     print('mmmm'); print(mm)
    #     break;
    # }

    
  }
  
  len = rep(0, NN);
  for(i in 1:NN){
    len[i] <- length(pts.final[[i]]);
  }

  
  plot(seq(1,NN,1), len,  type = 'o', col = "blue", main=c("Length"))
  return(list( block.size = block.full[[ind]][2], pts.1 = pts.final, pts.3 = pts.final[[ind]], len = len, block.full = block.full, omega.local = omega.local ) ) 
}



bss <- function(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), blocks, al){
  
  T <- length(data[,1]); k <- length(data[1,]); final.brk.points <- c(); pts.final <- c();
  # print(k);
  ######################################################################
  ######## FIRST STEP : INITIAL BRK POINT SELECTION ####################
  ######################################################################
  
  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  bbb <- floor(n.new/5);
  aaa <- sample(1:5, 1);
  # print(aaa)
  cv.index <- seq(aaa,n.new,floor(n.new/bbb)); # cv index
  temp.first <- first.step.cv.new.blocks(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), cv.index, blocks=blocks)
  mspe.plot.new(method,temp.first$cv,lambda.1.cv,temp.first$cv.final, "1")
  fisrt.brk.points <- temp.first$brk.points;
  print(fisrt.brk.points);
  
  
  ################## OMEGA LOCAL #######################################
  n <- T - p;
  # omega.local <- floor(1*mean(blocks.size))+1*floor((1/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1;# the penalty term in the information criterion -- the higher omega, the smaller number of break points selected.
  # omega.local <- floor(1*mean(blocks.size))+floor(( log( (T^0.5)/(mean(blocks.size)) + 2)  )*(((log(n))^1.5)*(log(k))^(3/2)))+1;# the penalty term in the information criterion.
  
  # change here!! this is the original setting!!
  # omega.local <- floor(1*mean(blocks.size))+floor(( log( (T^0.5)/(mean(blocks.size)) + 2)  )*(((log(n))^1.5)*(log(k))^(3/2)))+1;# the penalty term in the information criterion.
  # omega.local <- floor((1.25)*omega.local);
  
  # omega.local <- 3800 + floor(1*mean(blocks.size))     #### change for application case
  omega.local <- 75 + floor(1*mean(blocks.size))       #### change for simulation case: E1 is 200; E2 is 500 or 1000
  
  # omega.local <- min(omega.local,0.5*(min(fisrt.brk.points)-1-p),0.5*(n - max(fisrt.brk.points)-1));
  
  
  remove.ind <- c();
  if(length(fisrt.brk.points) != 0){
    for(i in 1:length(fisrt.brk.points)){
      if ( fisrt.brk.points[i] < (omega.local-1-p)   ){remove.ind <- c(remove.ind,i);}
      if ( (T-fisrt.brk.points[i]) < (omega.local-1-p)   ){remove.ind <- c(remove.ind,i);}
    }
  }
  
  if( length(remove.ind) > 0  ){fisrt.brk.points <- fisrt.brk.points[-remove.ind];}
  
  # omega.local <- floor(min(omega.local,0.975*(min(fisrt.brk.points)-1-p),0.975*(n - max(fisrt.brk.points)-1)));
  print("omega.local"); print(omega.local);
  
  ######################################################################
  ######## SECOND STEP : LOCAL SCREENING            ####################
  ###################################################################### 
  if( length(fisrt.brk.points) != 0){

    # print(n)
    n.local <- 2*(length(fisrt.brk.points))*omega.local;
    omega <- floor((15/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1; # for blocks=100
    lambda.2 <- (1/1)*(log(2*omega.local)*log(k))/(2*omega.local) # the second tuning parameter. This default number seems to be working for many simulation and real data examples!
    # temp <- second.step.local.new(data, lambda = lambda.2, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), fisrt.brk.points, omega, omega.local, al)
    temp <- second.step.local.cluster(data, lambda = lambda.2, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), fisrt.brk.points, omega, omega.local, al)
    final.brk.points <- temp$pts;
    print(final.brk.points)
    pts.final <- final.brk.points;
    
    while( min(abs(diff(pts.final))) <  2*omega.local  ){
      if( length(pts.final) != 0){
        pts.list <- block.finder(pts.final,2*omega.local)
        # print(pts.list);
        pts.final <- second.step.final(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local)
      }
    }
    
    print(pts.final)
    
    # if( length(final.brk.points) != 0){
    #   pts.list <- block.finder(final.brk.points,2*omega.local)
    #   # print(pts.list);
    #   pts.final <- second.step.final(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local)
    #   # print(pts.final)
    # }
    # 
    # if( length(pts.final) != 0){
    #   pts.list <- block.finder(pts.final,2*omega.local)
    #   # print(pts.list);
    #   pts.final <- second.step.final(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local)
    #   print(pts.final)
    # }
  }

  
  
  
  
 return(list(pts.1 = fisrt.brk.points, pts.2 = final.brk.points, pts.3 = pts.final, omega.local = omega.local)) 
}
























bss.cluster <- function(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), blocks, al){
  
  T <- length(data[,1]); k <- length(data[1,]); final.brk.points <- c(); pts.final <- c();
  
  ######################################################################
  ######## FIRST STEP : INITIAL BRK POINT SELECTION ####################
  ######################################################################
  
  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  bbb <- floor(n.new/5);
  aaa <- sample(1:5, 1)
  cv.index <- seq(aaa,n.new,floor(n.new/bbb)); # cv index
  temp.first <- first.step.cv.new.blocks(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), cv.index, blocks=blocks)
  mspe.plot(method,temp.first$cv,lambda.1.cv,temp.first$cv.final, "1")
  fisrt.brk.points <- temp.first$brk.points;
  print(fisrt.brk.points);
  
  ######################################################################
  ######## SECOND STEP : LOCAL SCREENING            ####################
  ###################################################################### 
  if( length(fisrt.brk.points) != 0){
    n <- T - p;
    omega.local <- floor(1*mean(blocks.size))+1*floor((1/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1;# the penalty term in the information criterion -- the higher omega, the smaller number of break points selected.
    omega.local <- 2*omega.local;
    omega.local <- min(omega.local,0.5*(min(fisrt.brk.points)-1-p),0.5*(n - max(fisrt.brk.points)-1));
    omega.local <- floor(min(omega.local,0.95*(min(fisrt.brk.points)-1-p),0.95*(n - max(fisrt.brk.points)-1)));
    # print(omega.local)
    # print(n)
    n.local <- 2*(length(fisrt.brk.points))*omega.local;
    omega <- floor((15/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1; # for blocks=100
    lambda.2 <- (1/1)*(log(2*omega.local)*log(k))/(2*omega.local) # the second tuning parameter. This default number seems to be working for many simulation and real data examples!
    temp <- second.step.local.new(data, lambda = lambda.2, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), fisrt.brk.points, omega, omega.local, al)
    final.brk.points <- temp$pts;
    print(final.brk.points)
    
    if( length(final.brk.points) != 0){
      # pts.list <- block.finder(final.brk.points,omega.local)
      pts.list <- block.finder(final.brk.points,3*floor(min(blocks.size)))
      pts.final <- second.step.final(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local)
      print(pts.final)
    }
  }
  
  
  
  
  
  return(list(pts.1 = fisrt.brk.points, pts.2 = final.brk.points, pts.3 = pts.final)) 
}


bss.test <- function(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), blocks, al){
  
  n <- length(data[,1]); k <- length(data[1,]); final.brk.points <- c(); pts.final <- c();
  n.new <- length(blocks) - 1;
  ######################################################################
  ######## FIRST STEP : INITIAL BRK POINT SELECTION ####################
  ######################################################################
  

  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  # bbb <- floor(n.new/5);
  # aaa <- sample(1:5, 1)
  # cv.index <- seq(aaa,n.new,floor(n.new/bbb)); # cv index
  # temp.first <- first.step.cv.new.blocks(method, data, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, tol = tol, step.size = (1/2)*10^(-4), cv.index, blocks=blocks)
  # mspe.plot(method,temp.first$cv,lambda.1.cv,temp.first$cv.final, "1")
  # fisrt.brk.points <- temp.first$brk.points;
  # print(fisrt.brk.points);
  
  bnd.1 <- seq(1,floor(T^{0.25}),1);
  bnd.2 <- seq(n.new+1, n.new+1 - floor(T^{0.25}), -1)
  fisrt.brk.points <- blocks[-c(bnd.1,bnd.2)];
  # print(fisrt.brk.points);
  ######################################################################
  ######## SECOND STEP : LOCAL SCREENING            ####################
  ###################################################################### 
  if( length(fisrt.brk.points) != 0){
    n <- T - p;
    omega.local <- floor(1*mean(blocks.size))+1*floor((1/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1;# the penalty term in the information criterion -- the higher omega, the smaller number of break points selected.
    omega.local <- 2*omega.local;
    omega.local <- min(omega.local,0.5*(min(fisrt.brk.points)-1-p),0.5*(n - max(fisrt.brk.points)-1));
    omega.local <- floor(min(omega.local,0.95*(min(fisrt.brk.points)-1-p),0.95*(n - max(fisrt.brk.points)-1)));
    print(omega.local)
    n.local <- 2*(length(fisrt.brk.points))*omega.local;
    omega <- floor((15/1)*(((log(n))^1.5)*(log(k))^(3/2)))+1; # for blocks=100
    lambda.2 <- (1/1)*(log(2*omega.local)*log(k))/(2*omega.local) # the second tuning parameter. This default number seems to be working for many simulation and real data examples!
    temp <- second.step.local.new(data, lambda = lambda.2, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), fisrt.brk.points, omega, omega.local, al)
    final.brk.points <- temp$pts;
    print(final.brk.points)
    
    if( length(final.brk.points) != 0){
      pts.list <- block.finder(final.brk.points,omega.local)
      pts.final <- second.step.final(data, p, max.iteration = 1000, tol = tol, step.size = 10^(-3), pts.list, omega.local)
      print(pts.final)
    }
  }
  
  
  
  
  
  return(list(pts.1 = fisrt.brk.points, pts.2 = final.brk.points, pts.3 = pts.final)) 
}



################# HAOREN FUNCTIONS #########################################

get.thr <- function(z, B=200, p=NULL, scales, do.parallel=4, method=c('ar', 'sb')[1], alpha=1){
  n <- dim(z)[1]; len <- dim(z)[2]
  if(method=='ar'){
    if(is.null(p)) p <- round(log(len))
    coef <- matrix(0, nrow=n, ncol=p)
    res <- matrix(0, nrow=n, ncol=len-p)
    for(i in 1:n){
      ar.fit <- ar(z[i, ], order.max=p, method='yw')
      if(length(ar.fit$ar) > 0) coef[i, 1:length(ar.fit$ar)] <- ar.fit$ar
      res[i, ] <- ar.fit$resid[-(1:p)]
    }
  }
  if(method=='sb'){
    if(is.null(p) || !(p>0 & p<1)){
      p <- apply(z, 1, function(w){g <- get.gg(w); min(.5, 1/(((g[2]/g[1])^2)^(1/3)*len^(1/5)))})
    } else p <- rep(p, n)
  }
  if(do.parallel){
    cl <- parallel::makeCluster(do.parallel); doParallel::registerDoParallel(cl)
  }
  null.stat <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo', 'factorcpt')) %dopar% {
    if(method=='ar'){
      bz <- res[, sample(dim(res)[2], len+p, replace=TRUE)]
      for(t in 1:len) bz[, p+t] <- bz[, p+t] + apply(coef*bz[, t+(p-1):0], 1, sum)
      bz <- bz[, -(1:p)]
    }
    if(method=='sb'){
      bz <- z*0
      for(i in 1:n){
        ind <- c()
        while(length(ind)<len){
          L <- rgeom(1, p[i]); I <- sample(len, 1)
          ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
        }
        ind <- ind[1:len]
        bz[i, ] <- z[i, ind]
      }
    }
    eval <- c()
    for(sc in scales){
      cc <- factorcpt::func_coef(bz, sc)
      sgn <- sign(cc%*%t(cc))
      by <- t(factorcpt::func_input(cc, sgn))^2
      cs <- factorcpt::func_dc(by)$acs
      eval <- c(eval, apply(cs, 1, max))
    }
    eval
  }
  if(do.parallel) parallel::stopCluster(cl)
  
  # apply(null.stat, 1, function(w){quantile(w, 1)})
  apply(null.stat, 1, function(w){quantile(w, 1-alpha)})
}

sbs.alg <- function(x, thr.seq=NULL, dw=NULL, B=200, p=NULL, method=c('ar', 'sb')[1], scales=NULL, do.parallel=4){
  n <- nrow(x); T <- ncol(x); d <- n*(n+1)/2
  if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))
  if(is.null(dw)) dw <- round(max(.5*log(T)^2, 2^(-scales)))
  
  if(is.null(thr.seq)) thr.seq <- get.thr(x, B=B, p=p, scales=scales, do.parallel=do.parallel, method=method)
  
  est.cpts <- c()  
  for(i in 1:length(scales)){
    sc <- scales[i]
    cc <- factorcpt::func_coef(x, sc)
    sgn <- sign(cc%*%t(cc))
    y <- t(factorcpt::func_input(cc, sgn))^2
    len <- dim(y)[2]
    thr <- thr.seq[(i-1)*d+1:d]
    
    if(length(est.cpts)==0){
      mt <- make.tree(y, thr, dw)
      if(length(mt$est.cpts) > 0) est.cpts <- sort(mt$est.cpts)
    } else{
      brks <- c(0, est.cpts, len)
      for(b in 1:(length(brks)-1)){
        int <- (brks[b]+1):brks[b+1]
        mt <- make.tree(y[, int], thr, dw)
        est.cpts <- sort(c(est.cpts, brks[b]+mt$est.cpts))
      }
    }
  }
  est.cpts
  
}

search.chp <- function(stat, dw){
  b <- NULL; halt <- FALSE; test.stat <- 0
  len <- length(stat)
  stat[c(1:dw, (len-dw+1):len)] <- 0
  while(sum(stat) > 0){
    test.stat <- max(stat)
    b <- min(which(stat==test.stat))
    int <- max(1, b-dw):min(len, b+dw)
    if(sum(stat[int]==0)==0) break
    stat[int] <- 0
  }
  if(sum(stat)==0) halt<-TRUE
  return(list(halt=halt, test.stat=test.stat, b=b))
}

make.tree <- function(y, thr, dw){
  len <- dim(y)[2]
  tree <- list(matrix(0, 5, 1))
  est.cpts <- c()
  
  acs <- factorcpt::func_dc(y)$acs
  stat <- apply(acs*(acs > thr), 2, mean)
  sc <- search.chp(stat, dw)
  if(sc$halt) return(list(tree=tree, est.cpts=est.cpts))
  tree[[1]][1, 1] <- 1
  tree[[1]][2, 1] <- 1
  tree[[1]][3, 1] <- sc$b
  tree[[1]][4, 1] <- len
  tree[[1]][5, 1] <- sc$test.stat
  est.cpts <- c(est.cpts, tree[[1]][3, 1])
  
  j <- 1
  while(length(tree)==j){
    npc <- dim(tree[[j]])[2]
    if(sum(tree[[j]][4, ]-tree[[j]][2, ]-rep(4*dw, npc)>0)){
      ncc <- 0; i <- 1
      while(i <= npc){
        if(tree[[j]][3, i]-tree[[j]][2, i]+1>4*dw){
          s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
          acs <- factorcpt::func_dc(y[, s:e])$acs
          stat <- apply(acs*(acs > thr), 2, mean)
          sc <- search.chp(stat, dw)
          if(!sc$halt){
            if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]-1
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- s+sc$b-1
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- sc$test.stat
            est.cpts <- c(est.cpts, tree[[j+1]][3, ncc])
          }
        }
        if(tree[[j]][4, i]-tree[[j]][3, i]>4*dw){
          s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
          acs <- factorcpt::func_dc(y[, s:e])$acs
          stat <- apply(acs*(acs > thr), 2, mean)
          sc <- search.chp(stat, dw)
          if(!sc$halt){
            if(length(tree)==j) tree <- c(tree, list(matrix(0, 5, 0)))
            ncc <- ncc+1
            tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 5, 1)), 5, ncc)
            tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]
            tree[[j+1]][2, ncc] <- s
            tree[[j+1]][3, ncc] <- s+sc$b-1
            tree[[j+1]][4, ncc] <- e
            tree[[j+1]][5, ncc] <- sc$test.stat
            est.cpts <- c(est.cpts, tree[[j+1]][3, ncc])
          }
        }
        i <- i+1
      }
      j <- j+1
    } else break
  }
  list(tree=tree, est.cpts=est.cpts)
}

##

tri.kern <- function(h){
  filter <- rep(0, h+1)
  i <- 0
  while (i <= h) {
    u <- i/h
    if (u < 1/2)
      filter[i+1] <- 1
    if (u >= 1/2 & u < 1)
      filter[i+1] <- 2 * (1 - u)
    if (u > 1)
      break
    i <- i + 1
  }
  filter
}

get.gg <- function(z, M=NULL, C=2, max.K=5){
  len <- length(z)
  max.K <- max(max.K, sqrt(log(len)))
  acv <- acf(z, type="covariance", lag.max=len-1, plot=FALSE)$acf[,,1]
  if(is.null(M)){
    l <- 1; ind <- 0
    while(l < sqrt(len)){
      if(abs(acv[l+1])/acv[1] < C*sqrt(log(len)/len)){
        ind <- ind+1
      } else{
        if(ind>0) ind <- 0
      }
      if(ind==max.K) break
      l <- l+1
    }
    lam <- max(1/2, l-max.K); M <- 2*lam
  }
  k <- tri.kern(M)
  c(acv[1]+2*sum(k[-1]*acv[2:(M+1)]), 2*sum(k[-1]*(1:M)*acv[2:(M+1)]))
}


























