lasso.ss <- 
function(y, X, ss, a.lst, S.v, C.v, col.plot, verbose=TRUE, plot=TRUE) {
   check.ss.lasso(y, X, ss, a.lst, S.v, C.v, col.plot)   

   p <- dim(X)[2]
   i <- 0
   ind.v <- NULL
   col.ind <- rep(1:p)
   S.v2 <- S.v
   aic.v <- NULL
   names.X <- colnames(X)
   flag.stop <- FALSE
   flag.lasso <- FALSE
   stack.ss <- NULL
   stack.ss2 <- NULL

   A <- do.call(adiag, a.lst)

   X <- as.matrix(X)
   X.out <- as.matrix(X)
   X.in <- matrix(0, nrow(X.out),1)
   X.cand <- X%*%A

   # Step 1: Initialize all regression coefficient estimates equal to 0, and let r=y-ybar.
   if (verbose) {
      cat("Step 1: Set r=y and initialize all betas to 0.", "\n")
   }
   r <- y # Note: ybar=0 since all variables were standardized to have mean=0 and SD=1.
   beta.old <- array(0, dim=c(1,p), dimnames=list(NULL, names.X))
   row.names(beta.old) <- i
   
   # Step 2: Find the predictor among the p possible predictors that has the greatest absolute correlation with the residuals.
   i <- i + 1
   c.hat.v <- t(X.cand) %*%r
   C.hat <- max(abs(c.hat.v))
   index.maxabscor <- which.max(abs(c.hat.v))
   ind.v <- col.ind[index.maxabscor]   
   
   if (verbose) {
      cat("Step 2: Find the predictor that has the greatest absolute correlation with r.", "\n")
      cat("\t", "Note that LARS picks the maximal absolute current correlation.", "\n", "\n")
   }

   a.lst.ind <- NULL
   for (j in 1:length(a.lst)) {
      if (length(which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[index.maxabscor]))==0) {
         a.lst.ind[[j]] <- -1
      } else {
         a.lst.ind[[j]] <- which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[index.maxabscor])
      }
   }

   index.lst <- which(a.lst.ind != -1)
   index.col <- unlist(a.lst.ind[a.lst.ind != -1])

   col.ind.v <- NULL
   sub <- NULL

   for (j in 1:S.v[index.lst]) {
      col.ind.v[j] <- which(colnames(X)==dimnames(a.lst[[index.lst]])[2][[1]][j])
      sub[j] <- which(col.ind==col.ind.v[j])
   }

   col.ind <- col.ind[-sub]
   C.v[index.lst] <- index.col
   a.lst[[index.lst]] <- a.lst[[index.lst]][,index.col, drop=FALSE]
 
   if (S.v[index.lst]==1) {
      stack.ss <- 1
      stack.ss2 <- 1
   } else {
      stack.ss <- which(pmatch(ss, dimnames(a.lst[[index.lst]])[2])!="NA")
      stack.ss2 <- which(pmatch(ss, dimnames(a.lst[[index.lst]])[2])!="NA")
   }

   A <- do.call(adiag, a.lst)
   X.cand <- X%*%A

   if (verbose) {
      cat("LARS Iteration ", i, sep="", "\n")
      cat("\t", "Variable ", tail(ind.v, n=1), " added", sep="", "\n")
   }

   # Step 3: Compute c.hat.v and C.hat.
   if (verbose) {
      cat("\t", " Step 3: Compute c.hat.v and C.hat.", "\n")
   }
   c.hat.v <- t(X.cand) %*%r
   print(c.hat.v)
   C.hat <- max(abs(c.hat.v))

   # Step 4: Set s.j=sign{c.hat.v}, and calculate X.in, A.in, u.in, and a.
   if (verbose) {
      cat("\t", " Step 4: Set s.j=sign{c.hat.v}, and calculate X.in, A.in, u.in, and a.", "\n")
   }

   X.in <- X[,ind.v, drop=FALSE]

   gamma.ind.v <- NULL
   for (j in 1:dim(X.in)[2]) {
      gamma.ind.v[j] <- which(colnames(X.cand)==colnames(X.in)[j])
   } 

   X.in <- t(t(X.in)*sign(c.hat.v[gamma.ind.v]))

   a.lst.out <- a.lst
   for (j in 1:length(C.v)) {
      if (C.v[j]!=0) {
         a.lst.out[[j]] <- diag(NA, S.v[j])
      }
   }

   A.out <- do.call(adiag, a.lst.out)
   colnames(A.out) <- colnames(X)
 
   X.out <- (X%*%A.out)[,apply((X%*%A.out),2,function(x){!all(is.na(x))})] 

   g.in <- t(X.in) %*%X.in
   inv.g.in <- solve(t(X.in) %*%X.in)
   A.in <- as.vector((t(rep(1,dim(X.in)[2])) %*% inv.g.in %*% rep(1,dim(X.in)[2]))^(-1/2))
   w.in <- A.in * inv.g.in %*% rep(1,dim(X.in)[2]) 
   u.in <- X.in %*% w.in
   a <- t(X.cand) %*% u.in
 
   # Step 5: Calculate gamma.hat and d.v.  
   if (verbose) {
      cat("\t", " Step 5: Calculate gamma.hat and d.v.", "\n")
   }

   gamma.hat.v <- rep(0,p)

   for (j in col.ind) {
      comp1 <- (C.hat-(A%*%c.hat.v)[j])/(A.in-(A%*%a)[j])
      comp2 <- (C.hat+(A%*%c.hat.v)[j])/(A.in+(A%*%a)[j])
      comp <- c(comp1, comp2)
      gamma.hat.v[j] <- min(comp[comp>0])
   }

   gamma.hat <- min(gamma.hat.v[gamma.hat.v>0])
   cat("maximal abs current cor", C.hat-gamma.hat*A.in, "\n")
   gamma.hat.index <- which.min(gamma.hat.v[gamma.hat.v>0])

   if (verbose) {
      cat("\t", "\t", " Gamma values:", gamma.hat.v, "\n")
      cat("\t", "\t", " Gamma hat:", gamma.hat, "\n")
   }

   sign.v <- sign(c.hat.v)
   d.v <- rep(0,p)
   d.v[ind.v] <- sign.v[gamma.ind.v] * w.in

   # Step 6: Update beta.i <- beta.i-1 + gamma.hat*d.v, and update r.i <- y-X*beta.i.
   beta.new <- beta.old[i,] + gamma.hat*d.v

   if (verbose) {
      cat("\t", " Step 6: Update betas and r.", "\n")
      cat("\t", "\t", " Beta.new:", beta.new, "\n", "\n")  
   }

   r <- y - X.cand%*%t(beta.new%*%A)
   beta.old <- rbind(beta.old, beta.new)
   beta.old
   row.names(beta.old)[i+1] <- i

   X.select <- X[,sort(ind.v)]
   OLS <- lm(y~as.matrix(X.select)-1)
   aic.v[[i]] <- AIC(OLS)

   while(flag.stop!=TRUE) {
      i <- i + 1
      # Repeat steps 3-6 until flag.stop=TRUE.
      v1 <- pickvar.lasso.ss(r, y, X, X.cand, X.in, X.out, ss, a.lst, A, S.v, S.v2, C.v, beta.old, ind.v, col.ind, gamma.hat.index, gamma.tilda.index, flag.lasso, i, p, stack.ss, stack.ss2, verbose)
      r <- v1$r
      X.in <- v1$X.in
      X.out <- v1$X.out
      X.cand <- v1$X.cand
      a.lst <- v1$a.lst
      A <- v1$A
      S.v <- v1$S.v
      S.v2 <- v1$S.v2
      C.v <- v1$C.v

      beta.new <- v1$beta.new
      beta.old <- rbind(beta.old, beta.new)
      row.names(beta.old)[i+1] <- i
      ind.v <- v1$ind.v
      col.ind <- v1$col.ind

      X.select <- X[,sort(ind.v)]
      OLS <- lm(y~as.matrix(X.select)-1)
      aic.v[[i]] <- AIC(OLS)

      gamma.hat.index <- v1$gamma.hat.index
      gamma.tilda.index <- v1$gamma.tilda.index
      stack.ss <- v1$stack.ss
      stack.ss2 <- v1$stack.ss2
      flag.lasso <- v1$flag.lasso
      flag.stop <- v1$flag.stop
  
      if (flag.stop==TRUE) {
         beta.nonzero <- beta.old%*%A
         aic.sol <- which.min(aic.v) 
         beta.nonzero.aic <- beta.nonzero[(aic.sol+1),which(!beta.nonzero[(aic.sol+1),]==0)]
         print(beta.nonzero[(i+1), ])
         print(paste("No. of vars in model = ", dim(beta.nonzero)[2], sep=""))
         print(paste("Total no. of vars considered = ", p, sep=""))
         if (plot==TRUE) {
            ind.stack <- cbind(ind.v, stack.ss2)
            stack.ss3 <- sortrows(ind.stack)[,2]
            path.plot.lasso.ss(beta.nonzero, ind.v, i, stack.ss3, aic.sol, col.plot)
         }
         break
      }
   }
return(list(beta=beta.nonzero, beta.aic=beta.nonzero.aic, ind.v=ind.v, aic.v=aic.v, stack.ss=stack.ss2))
}