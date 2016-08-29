### Spatial scale forward stepwise regression ###

# Helper function to check user inputs for SS stepwise
check.ss.step <- function(y, X.3D, y.name, ss, epsilon) {
   if (is_numeric_vector(y)==FALSE) stop('y must be a numeric vector')
   if (round(mean(y), 3) != 0 | round(sd(y), 3) != 1) stop('y must be standardized to have mean=0 and sd=1')

   if (is.numeric(X.3D)==FALSE | is.array(X.3D)==FALSE) stop('X.3D must be a numeric array')
   if (dim(X.3D)[3] != length(ss)) stop('X.3D must have same number of stacks as number of ss')
   check.sd <- NULL
   for (level in 1:dim(X.3D)[3]) {
     check.sd[level] <- any(round(apply(X.3D[,,level], 2, sd), 3) != 1, na.rm=TRUE)
   }
   if (any(round(colMeans(X.3D), 3) != 0, na.rm=TRUE) | any(check.sd != FALSE)) stop('X.3D must be standardized to have mean=0 and sd=1')

   if (is.character(y.name)==FALSE) stop('y.name must be character')

   if (is_string_vector(ss)==FALSE) stop('ss must be a character vector')
   if (length(ss) < 1 | length(ss) > 100 ) stop('number of ss must be between 1-100')

   if (epsilon <= 0) stop('epsilon must be greater than 0')
}


# Description: Computes the resultant regression model, the associated AIC, and the residuals.
get.results.step.ss <- function(y, X.in, X.out, y.name, seq.v, ss, verbose, i) {
   out <- lm(y ~ .-1, data=X.in, x=TRUE, y=TRUE) # omitting intercept
   beta <- out$coefficients
   aic.v <- extractAIC(out)
   r <- out$residuals

   if (verbose) {
      print(i)
      X.in.names <- names(X.in)
      print(paste("Step: AIC = ", round(aic.v[2], 4), sep="" ))

      if (i==1) {
         print(paste(y.name, " ~ ", X.in.names, sep=""))
      } else if (i>1) {
         cat(y.name, paste("~"), X.in.names[1], paste("+", X.in.names[-1], sep = " "), "\n")
      }
   }

   return(list(out=out, beta=beta, aic.v=aic.v, r=r, X.in=X.in, X.out=X.out, seq.v=seq.v))
}

 

# Description: Finds the predictor that has the largest absolute correlation with the given outcome and adds it to the working design matrix.
pickvar.step.ss <- function(outcome, y, X.in, X.out, y.name, seq.v, ss, stack.ss, verbose, i, k) {
   flag.stack <- NULL
   maxcor <- NULL   
   col <- NULL

   if (dim(X.out)[2]!=1) {    
      for (j in seq.v) {
         flag.stack[j] <- ifelse(sum(is.na(X.out[,,j]))==dim(X.out[,,j])[1]* dim(X.out[,,j])[2], 1, 0)
      }
   } else {
      for (j in seq.v) {
         flag.stack[j] <- ifelse(sum(is.na(X.out[,,j]))==dim(X.in)[1], 1, 0)
      }

   }

   flag.stack <- na.omit(flag.stack)

   if (sum(flag.stack)!=0) {
       seq.v <- seq.v[-(which(flag.stack==1))]
   }
  
   for(m in seq.v){
      maxcor[m] <- max(abs(cor(X.out[,,m],outcome)), na.rm=TRUE)
      col[m] <- which.max(abs(cor(X.out[,,m],outcome)))
   }

   index <- cbind(maxcor, col)
   index.final <- c(index[which.max(index[,1]), 2], which.max(index[,1]))
   stack.ss[i] <- which.max(index[,1])


   if (i==1) {
      X.in <- as.data.frame(X.out[,index.final[1], index.final[2], drop=FALSE])
      colnames(X.in) <- paste(ss[index.final[2]], colnames(X.in), sep="_")
   } else if (i>1) {
      X.in.update <- as.data.frame(X.out[,index.final[1], index.final[2], drop=FALSE])
      colnames(X.in.update) <- paste(ss[index.final[2]], colnames(X.in.update), sep="_")
      X.in <- cbind(X.in, X.in.update)
   }

   X.out <- X.out[,-index.final[1],,drop=FALSE]   

   results <- get.results.step.ss(y, X.in, X.out, y.name, seq.v, ss, verbose, i)
   return(list(out=results$out, beta=results$beta, aic.v=results$aic.v, r=results$r, X.in=results$X.in, X.out=results$X.out, seq.v=results$seq.v, stack.ss=stack.ss))
}







### Spatial scale forward stagewise regression ###

# Helper function to check user inputs for SS stagewise
check.ss.stage <- function(y, X, X.3D, ss, increment, tolerance, col.plot) {
   if (is_numeric_vector(y)==FALSE) stop('y must be a numeric vector')
   if (round(mean(y), 3) != 0 | round(sd(y), 3) != 1) stop('y must be standardized to have mean=0 and sd=1')

   if (is_numeric_dataframe(X)==FALSE) stop('X must be a numeric data frame')
   if (any(round(colMeans(X), 3) != 0) | any(round(apply(X, 2, sd), 3) != 1)) stop('X must be standardized to have mean=0 and sd=1')

   if (is.numeric(X.3D)==FALSE | is.array(X.3D)==FALSE) stop('X.3D must be a numeric array')
   if (dim(X.3D)[3] != length(ss)) stop('X.3D must have same number of stacks as number of ss')
   check.sd <- NULL
   for (level in 1:dim(X.3D)[3]) {
     check.sd[level] <- any(round(apply(X.3D[,,level], 2, sd), 3) != 1, na.rm=TRUE)
   }
   if (any(round(colMeans(X.3D), 3) != 0, na.rm=TRUE) | any(check.sd != FALSE)) stop('X.3D must be standardized to have mean=0 and sd=1')

   if (is_string_vector(ss)==FALSE) stop('ss must be a character vector')
   if (length(ss) < 1 | length(ss) > 100 ) stop('number of ss must be between 1-100')

   if (increment <= 0) stop('increment must be greater than 0')

   if (tolerance <= 0) stop('tolerance must be greater than 0')

   if (is_string_vector(col.plot)==FALSE) stop('col.plot must be a character vector')
   if (length(col.plot) != length(ss)) stop('col.plot must be of same length as length of ss')
}


# Description: Evaluates whether the max|cor(r,X)| < tolerance and returns flag indicator as TRUE or FALSE.
stop.stage.ss <- function(r, X.cand, seq.v, i, k, tolerance, verbose) {
   flag.stack <- NULL
   maxcor <- NULL
   col <- NULL

   for (j in seq.v) {
      flag.stack[j] <- ifelse(sum(is.na(X.cand[,,j]))==dim(X.cand[,,j])[1]* dim(X.cand[,,j])[2], 1, 0)
   }

   flag.stack <- na.omit(flag.stack)

   if (sum(flag.stack)!=0) {
       seq.v <- seq.v[-(which(flag.stack==1))]
   }

   for (m in seq.v){
      maxcor[m] <- max(abs(cor(X.cand[,,m],r)), na.rm=TRUE)
      col[m] <- which.max(abs(cor(X.cand[,,m],r)))
   }

   index <- cbind(maxcor, col)
   maxcor.final <- max(na.omit(index[,1]))

   flag <- ifelse(maxcor.final < tolerance, TRUE, FALSE) 
   if (flag==TRUE & verbose==TRUE) {
      print(paste("Iteration = ", i, sep=""))
      print(paste("Max|cor| = ", format(maxcor.final, digits=3, nsmall=2), sep=""))   
      print(paste("Tolerance = ", tolerance, sep=""))           
   }
   return(list(flag=flag, seq.v=seq.v))      
}


# Description: Finds the predictor that has the greatest absolute correlation with the residuals and then updates the betas and r.
pickvar.stage.ss <- function(r, X.cand, seq.v, ss, beta.old, i, p, k, names.X, stack.ss, increment, tolerance, verbose) {
   # Step 2: Find the predictor x.j that has the greatest absolute correlation with the residuals.
   maxcor <- NULL
   col <- NULL
   
   for(j in seq.v){
      maxcor[j] <- max(abs(cor(r,X.cand[,,j])), na.rm=TRUE)
      col[j] <- which.max(abs(cor(r,X.cand[,,j])))
   }

   index <- cbind(maxcor, col)
   index.final <- c(index[which.max(index[,1]), 2], which.max(index[,1]))
   stack.ss[index.final[1]] <- which.max(index[,1])

   seq.v.remove <- rep(1:k)[-(index.final[2])]
   for (j in seq.v.remove) {
      X.cand[,index.final[1],j] <- NA
   }

   # Step 3: Let beta.j.hat <- beta.j.hat + delta.j, where delta.j=increment*sign[cor(r,x.j)].
   ind <- rep(0,p) 
   delta <- increment*sign(cor(r,X.cand[,index.final[1],index.final[2]]))
   ind[index.final[1]] <- 1
   beta.new <- beta.old[i,] + delta*ind

   if (index.final[[2]] != 1) {
      name.X.update <- paste(ss[index.final[2]], colnames(X.cand)[index.final[1]], sep="_")
   } else if (index.final[[2]] == 1) {
      name.X.update <- colnames(X.cand)[index.final[1]]
   }

   names.X[index.final[1]] <- name.X.update
  
   # Step 4: Let r <- r - delta.j*x.j.
   r <- r - (delta*X.cand[,index.final[1], index.final[2]])
   
   results <- stop.stage.ss(r, X.cand, seq.v, i, k, tolerance, verbose)
   return(list(r=r, beta.new=beta.new, names.X=names.X, X.cand=X.cand, stack.ss=stack.ss, flag=results$flag, seq.v=results$seq.v))
}


# Description: Plots the coefficient profile.
path.plot.stage.ss <- function(beta.nonzero, i, stack.ss, increment, tolerance, path.index, col.plot) {  
   iteration <- rep(0:i)
   colors.v <- col.plot
   plot(c(0,i), c(min(beta.nonzero), max(beta.nonzero)), type="n", xlab="Iteration", ylab="Coefficients") 
   abline(h=0, lty=3)
   title(print(paste("Increment = ", increment, ", Tolerance = ", tolerance, sep="")))
   for (j in 1:dim(beta.nonzero)[2]) {
      lines(iteration, beta.nonzero[,j], type="s", col=colors.v[stack.ss[j]])
      #axis(4, at=beta.nonzero[(i+1),j], labels=colnames(beta.nonzero)[j], las=2, cex.axis=0.7, tck=-.01)
      axis(4, at=beta.nonzero[(i+1),j], labels=path.index[j], las=2, cex.axis=0.7, tck=-.01)
   }
}







### Spatial scale LARS ###

# Helper function to check user inputs for SS LARS
check.ss.LARS <- function(y, X, ss, a.lst, S.v, C.v, col.plot) {
   if (is_numeric_vector(y)==FALSE) stop('y must be a numeric vector')
   if (round(mean(y), 3) != 0 | round(sd(y), 3) != 1) stop('y must be standardized to have mean=0 and sd=1')

   if (is_numeric_dataframe(X)==FALSE) stop('X must be a numeric data frame')
   if (any(round(colMeans(X), 3) != 0) | any(round(apply(X, 2, sd), 3) != 1)) stop('X must be standardized to have mean=0 and sd=1')

   if (is_string_vector(ss)==FALSE) stop('ss must be a character vector')
   if (length(ss) < 1 | length(ss) > 100 ) stop('number of ss must be between 1-100')

   if (is.list(a.lst)==FALSE) stop('a.lst must be a list')

   if (length(S.v) != length(a.lst) | any(S.v < 0)) stop('S.v must have same length as a.lst and have positive values')

   if (length(C.v) != length(a.lst) | any(C.v != 0)) stop('C.v must have same length as a.lst and have all values initialized to 0')
   
   if (is_string_vector(col.plot)==FALSE) stop('col.plot must be a character vector')
   if (length(col.plot) != length(ss)) stop('col.plot must be of same length as length of ss')
}


# Description: Evaluates whether all variables are in active set 
# and then returns flag indicator as TRUE or FALSE.
stop.lars_v2.ss <- function(X.in, S.v, i, p, verbose) {
   flag <- ifelse(dim(X.in)[2]==length(S.v), TRUE, FALSE) 
   return(flag)      
}


# Description: Finds gamma and d.v and then updates the betas and r.
pickvar.lars_v2.ss <- function(r, y, X, X.cand, X.in, X.out, ss, a.lst, A, S.v, C.v, beta.old, ind.v, col.ind, gamma.hat.index, i, p, stack.ss, verbose) {
   a.lst.ind <- NULL
   for (j in 1:length(a.lst)) {
      if (length(which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[gamma.hat.index]))==0) {
         a.lst.ind[[j]] <- -1
      } else {
         a.lst.ind[[j]] <- which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[gamma.hat.index])
        }
   }

   index.lst <- which(a.lst.ind != -1)
   index.col <- unlist(a.lst.ind[a.lst.ind != -1])

   col.ind.v <- NULL
   sub <- NULL

   for (j in 1:S.v[index.lst]) {
      col.ind.v[j] <- which(colnames(X)==dimnames(a.lst[[index.lst]])[2][[1]][j] )
      sub[j] <- which(col.ind==col.ind.v[j])
   }

   col.ind <- col.ind[-sub]
   C.v[[index.lst]] <- index.col
   a.lst[[index.lst]] <- a.lst[[index.lst]][,index.col, drop=FALSE]
   ind.v <- c(ind.v, which(colnames(X)==dimnames(a.lst[[index.lst]])[2]))

   if (verbose) {
      cat("LARS Iteration ", i, sep="", "\n")
      cat("\t", "Variable ", tail(ind.v, n=1), " added", sep="", "\n")
   }

   if (S.v[index.lst]==1) {
      stack.ss <- c(stack.ss, 1)
   } else {
      stack.ss <- c(stack.ss, which(pmatch(ss, dimnames(a.lst[[index.lst]])[2])!="NA"))      
   }
 
   A <- do.call(adiag, a.lst)
   X.cand <- X%*%A

   if (verbose) {
      cat("\t", " Step 3: Compute c.hat.v and C.hat.", "\n")
   }
   c.hat.v <- t(X.cand) %*%r
   print(c.hat.v)

   C.hat <- max(abs(c.hat.v))

   if (verbose) {
      cat("\t", " Step 4: Set s.j=sign{c.hat.v}, and calculate X.in, A.in, u.in, and a.", "\n")
   }

   X.in <- X[,ind.v]

   gamma.ind.v <- NULL
   for (i in 1:dim(X.in)[2]) {
      gamma.ind.v[i] <- which(colnames(X.cand)==colnames(X.in)[i])
   } 

   X.in <- t(t(X.in)*sign(c.hat.v[gamma.ind.v]))
   X.out <- X[,col.ind, drop=FALSE]

   g.in <- t(X.in) %*%X.in
   inv.g.in <- solve(t(X.in) %*%X.in)
   A.in <- as.vector((t(rep(1,dim(X.in)[2])) %*% inv.g.in %*% rep(1,dim(X.in)[2]))^(-1/2))
   w.in <- A.in * inv.g.in %*% rep(1,dim(X.in)[2]) 
   u.in <- X.in %*% w.in
   a <- t(X.cand) %*% u.in
 
   if (verbose) {
      cat("\t", " Step 5: Calculate gamma.hat and d.v.", "\n")
   }

   if (dim(X.in)[2]!=length(S.v)) {
      gamma.hat.v <- rep(0,p)
      for (j in col.ind) {
         comp1 <- (C.hat-(A%*%c.hat.v)[j])/(A.in-(A%*%a)[j])
         comp2 <- (C.hat+(A%*%c.hat.v)[j])/(A.in+(A%*%a)[j])
         comp <- c(comp1,comp2)
         gamma.hat.v[j] <- min(comp[comp>0])
      }
      gamma.hat <- min(gamma.hat.v[gamma.hat.v>0])
      gamma.hat.index <- which.min(gamma.hat.v[gamma.hat.v>0])
   } else if (dim(X.in)[2]==length(S.v)) {
      gamma.hat <- C.hat/A.in
   }

   if (verbose) {
      if (dim(X.in)[2]!=length(S.v)) {
         cat("\t", "\t", " Gamma values:", gamma.hat.v, "\n")
         cat("\t", "\t", " Gamma hat:", gamma.hat, "\n")
      } else if (dim(X.in)[2]==length(S.v)) {
         cat("\t", "\t", " Gamma hat:", gamma.hat, "\n")
      }
   }

   sign.v <- sign(c.hat.v)
   d.v <- rep(0,p)
   d.v[ind.v] <- sign.v[gamma.ind.v] * w.in
   beta.new <- beta.old[i,] + gamma.hat*d.v

   if (verbose) {
      cat("\t", "  Step 6: Update betas and r.", sep="", "\n")
      cat("\t", "\t", " Beta.new:", beta.new, "\n")  
   }

   r <- y - X.cand%*%t(beta.new%*%A)

   flag.stop <- stop.lars_v2.ss(X.in, S.v, i, p, verbose)
   return(list(r=r, X.in=X.in, X.out=X.out, X.cand=X.cand, a.lst=a.lst, A=A, C.v=C.v, beta.new=beta.new, ind.v=ind.v, col.ind=col.ind, gamma.hat.index=gamma.hat.index, stack.ss=stack.ss, flag.stop=flag.stop))  
}


# Description: Plots the coefficient profile.
path.plot.lars_v2.ss <- function(beta.nonzero, ind.v, i, stack.ss2, aic.sol, col.plot) {
   iteration <- rep(0:i)
   colors.v <- col.plot
   plot(c(0,i), c(min(beta.nonzero), max(beta.nonzero)), type="n", xlab="Iteration", ylab="Coefficients") 
   abline(h=0, lty=3)
   abline(v=aic.sol, lty=3)
   for (j in 1:dim(beta.nonzero)[2]) {
      lines(iteration, beta.nonzero[,j], col=colors.v[stack.ss2[j]])
      axis(4, at=beta.nonzero[(i+1),j], labels=sort(ind.v)[j], las=2, cex.axis=0.7, tck=-.01)
   }
}







### Spatial scale lasso ###

# Helper function to check user inputs for SS lasso
check.ss.lasso <- function(y, X, ss, a.lst, S.v, C.v, col.plot) {
   if (is_numeric_vector(y)==FALSE) stop('y must be a numeric vector')
   if (round(mean(y), 3) != 0 | round(sd(y), 3) != 1) stop('y must be standardized to have mean=0 and sd=1')

   if (is_numeric_dataframe(X)==FALSE) stop('X must be a numeric data frame')
   if (any(round(colMeans(X), 3) != 0) | any(round(apply(X, 2, sd), 3) != 1)) stop('X must be standardized to have mean=0 and sd=1')

   if (is_string_vector(ss)==FALSE) stop('ss must be a character vector')
   if (length(ss) < 1 | length(ss) > 100 ) stop('number of ss must be between 1-100')

   if (is.list(a.lst)==FALSE) stop('a.lst must be a list')

   if (length(S.v) != length(a.lst) | any(S.v < 0)) stop('S.v must have same length as a.lst and have positive values')

   if (length(C.v) != length(a.lst) | any(C.v != 0)) stop('C.v must have same length as a.lst and have all values initialized to 0')

   if (is_string_vector(col.plot)==FALSE) stop('col.plot must be a character vector')
   if (length(col.plot) != length(ss)) stop('col.plot must be of same length as length of ss')
}


# Description: Evaluates whether all variables are in active set 
# and then returns flag indicator as TRUE or FALSE.
stop.lasso.ss <- function(X.in, S.v, gamma.hat, gamma.tilda, i, p, verbose) {
   flag <- ifelse(dim(X.in)[2]==length(S.v) & gamma.hat<gamma.tilda, TRUE, FALSE) 
   return(flag)      
}


# Description: Finds gamma and d.v and then updates the betas and r.
pickvar.lasso.ss <- function(r, y, X, X.cand, X.in, X.out, ss, a.lst, A, S.v, S.v2, C.v, beta.old, ind.v, col.ind, gamma.hat.index, gamma.tilda.index, flag.lasso, i, p, stack.ss, stack.ss2, verbose) {
   if (flag.lasso==TRUE) {
      if (verbose) {
         cat("LASSO Iteration ", i, sep="", "\n")
         cat("\t", "Variable ", ind.v[gamma.tilda.index], " dropped", sep="", "\n")
      }

      a.lst.ind <- NULL
      for (j in 1:length(a.lst)) {
         if (length(which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.in)[gamma.tilda.index]))==0) {
            a.lst.ind[[j]] <- -1
         } else {
            a.lst.ind[[j]] <- which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.in)[gamma.tilda.index])
            }
      }

      index.lst <- which(a.lst.ind != -1)
      index.col <- unlist(a.lst.ind[a.lst.ind != -1])

      col.ind <- sort(c(col.ind, ind.v[gamma.tilda.index]))
      ind.v <- ind.v[-gamma.tilda.index]
      cat("ind.v", ind.v, "\n") 
      cat("col.ind", col.ind, "\n")

      S.v[index.lst] <- 1 
      C.v[index.lst] <- 0

      stack.ss <- stack.ss[-gamma.tilda.index]
      stack.ss2 <- stack.ss2[-gamma.tilda.index]

      A <- do.call(adiag, a.lst)
      X.cand <- X%*%A
    
      if (verbose) {
         cat("\t", " Step 3: Compute c.hat.v and C.hat.", "\n")
      }
      c.hat.v <- t(X.cand) %*%r
      print(c.hat.v)
      C.hat <- max(abs(c.hat.v))

      if (verbose) {
         cat("\t", " Step 4: Set s.j=sign{c.hat.v}, and calculate X.in, A.in, u.in, and a.", "\n")
      }

      X.in <- X[,ind.v, drop=FALSE]

      gamma.ind.v <- NULL
      for (j in 1:dim(X.in)[2]) {
         gamma.ind.v[j] <- which(colnames(X.cand)==colnames(X.in)[j])
      } 

      X.in <- t(t(X.in)*sign(c.hat.v[gamma.ind.v]))
      X.out <- X[,col.ind, drop=FALSE]

   } else if (flag.lasso==FALSE) {
      a.lst.ind <- NULL
      for (j in 1:length(a.lst)) {
         if (length(which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[gamma.hat.index]))==0) {
            a.lst.ind[[j]] <- -1
         } else {
            a.lst.ind[[j]] <- which(unlist(dimnames(a.lst[[j]])[2])==colnames(X.out)[gamma.hat.index])
           }
      }

      index.lst <- which(a.lst.ind != -1)
      index.col <- unlist(a.lst.ind[a.lst.ind != -1])

      col.ind.v <- NULL
      sub <- NULL

      for (j in 1:S.v[index.lst]) {
         col.ind.v[j] <- which(colnames(X)==dimnames(a.lst[[index.lst]])[2][[1]][j] )
         sub[j] <- which(col.ind==col.ind.v[j])
      }

      col.ind <- col.ind[-sub]
      C.v[index.lst] <- index.col
      a.lst[[index.lst]] <- a.lst[[index.lst]][,index.col, drop=FALSE]
      ind.v <- c(ind.v, which(colnames(X)==dimnames(a.lst[[index.lst]])[2]))
      cat("ind.v", ind.v, "\n") 
      cat("col.ind", col.ind, "\n")

      if (verbose) {
         cat("LARS Iteration ", i, sep="", "\n")
         cat("\t", "Variable ", tail(ind.v, n=1), " added", sep="", "\n")
      }

      if (S.v[index.lst]==1) {
         stack.ss <- c(stack.ss, 1)
      } else {
         stack.ss <- c(stack.ss, which(pmatch(ss, dimnames(a.lst[[index.lst]])[2])!="NA"))    
      }

      if (S.v2[index.lst]==1) {
         stack.ss2 <- c(stack.ss2, 1)
      } else {
         stack.ss2 <- c(stack.ss2, which(pmatch(ss, dimnames(a.lst[[index.lst]])[2])!="NA"))    
      }
 
      A <- do.call(adiag, a.lst)
      X.cand <- X%*%A

      if (verbose) {
         cat("\t", " Step 3: Compute c.hat.v and C.hat.", "\n")
      }
      c.hat.v <- t(X.cand) %*%r
      print(c.hat.v)
      C.hat <- max(abs(c.hat.v))

      if (verbose) {
         cat("\t", " Step 4: Letting s.j=sign{c.hat.v}, update X.in, A.in, and u.in and calculate a=t(X)*u.in.", "\n")
      }

      X.in <- X[,ind.v]

      gamma.ind.v <- NULL
      for (j in 1:dim(X.in)[2]) {
         gamma.ind.v[j] <- which(colnames(X.cand)==colnames(X.in)[j])
      } 

      X.in <- t(t(X.in)*sign(c.hat.v[gamma.ind.v]))
      X.out <- X[,col.ind, drop=FALSE]
   }

   g.in <- t(X.in) %*%X.in
   inv.g.in <- solve(t(X.in) %*%X.in)
   A.in <- as.vector((t(rep(1,dim(X.in)[2])) %*% inv.g.in %*% rep(1,dim(X.in)[2]))^(-1/2))
   w.in <- A.in * inv.g.in %*% rep(1,dim(X.in)[2]) 
   u.in <- X.in %*% w.in
   a <- t(X.cand) %*% u.in
 
   if (verbose) {
      cat("\t", " Step 5: Calculate gamma.hat and d.v.", "\n")
   }

   if (dim(X.in)[2]!=length(S.v)) {
      gamma.hat.v <- rep(0,p)
      for (j in col.ind) {
         comp1 <- (C.hat-(A%*%c.hat.v)[j])/(A.in-(A%*%a)[j])
         comp2 <- (C.hat+(A%*%c.hat.v)[j])/(A.in+(A%*%a)[j])
         comp <- c(comp1, comp2)
         gamma.hat.v[j] <- min(comp[comp>0])
      }

      gamma.hat <- min(gamma.hat.v[gamma.hat.v>0])
      cat("gamma.hat.v", gamma.hat.v, "\n")
      gamma.hat.index <- which.min(gamma.hat.v[gamma.hat.v>0])
      cat("gamma.hat", gamma.hat, "\n")

   } else if (dim(X.in)[2]==length(S.v)) {
      gamma.hat <- C.hat/A.in
      print(gamma.hat)
   }

   cat("maximal abs current cor", C.hat-gamma.hat*A.in, "\n")

   if (verbose) {
      if (dim(X.in)[2]!=length(S.v)) {
         cat("\t", "\t", " Gamma values:", gamma.hat.v, "\n")
         cat("\t", "\t", " Gamma hat:", gamma.hat, "\n")
      } else if (dim(X.in)[2]==length(S.v)) {
         cat("\t", "\t", " Gamma hat:", gamma.hat, "\n")
      }
   }

   sign.v <- sign(c.hat.v)
   d.v <- rep(0,p)
   d.v[ind.v] <- sign.v[gamma.ind.v] * w.in
   gamma.j <- -beta.old[i,][ind.v]/d.v[ind.v]
   gamma.tilda <- min(gamma.j[gamma.j>0])
   gamma.tilda.index <- which(gamma.j==min(gamma.j[gamma.j>0]))

   flag.lasso <- ifelse(gamma.tilda<gamma.hat, TRUE, FALSE)
   cat("beta.old[i,][ind.v]", beta.old[i,][ind.v], "\n")
   cat("gamma.j", gamma.j, "\n")
   cat("gamma.tilda", gamma.tilda, "\n")
   cat("gamma.tilda.index", gamma.tilda.index, "\n")

   if (flag.lasso==TRUE) {
      beta.new <- beta.old[i,] + gamma.tilda*d.v
      cat("beta.new[ind.v]", beta.new[ind.v], "\n")
      beta.new <- round(beta.new, 15)
      cat("beta.new[ind.v]", beta.new[ind.v], "\n")
   } else if (flag.lasso==FALSE) {
      beta.new <- beta.old[i,] + gamma.hat*d.v
      cat("beta.new[ind.v]", beta.new[ind.v], "\n")
      beta.new <- round(beta.new, 15)
      cat("beta.new[ind.v]", beta.new[ind.v], "\n")
   }

   if (verbose) {
      cat("\t", "  Step 6: Update betas and r.", sep="", "\n")
      cat("\t", "\t", " Beta.new:", beta.new, "\n")  
   }

   r <- y - X.cand%*%t(beta.new%*%A)

   flag.stop <- stop.lasso.ss(X.in, S.v, gamma.hat, gamma.tilda, i, p, verbose)
   return(list(r=r, X.in=X.in, X.out=X.out, X.cand=X.cand, a.lst=a.lst, A=A, S.v=S.v, S.v2=S.v2, C.v=C.v, beta.new=beta.new, ind.v=ind.v, col.ind=col.ind, gamma.hat.index=gamma.hat.index, gamma.tilda.index=gamma.tilda.index, stack.ss=stack.ss, stack.ss2=stack.ss2, flag.lasso=flag.lasso, flag.stop=flag.stop))  
}


# Description: Plots the coefficient profile.
path.plot.lasso.ss <- function(beta.nonzero, ind.v, i, stack.ss3, aic.sol, col.plot) {
   iteration <- rep(0:i)
   colors.v <- col.plot
   plot(c(0,i), c(min(beta.nonzero), max(beta.nonzero)), type="n", xlab="Iteration", ylab="Coefficients") 
   abline(h=0, lty=3)
   abline(v=aic.sol, lty=3) 
   for (j in 1:dim(beta.nonzero)[2]) {
      lines(iteration, beta.nonzero[,j], col=colors.v[stack.ss3[j]])
      axis(4, at=beta.nonzero[(i+1),j], labels=sort(ind.v)[j], las=2, cex.axis=0.7, tck=-.01)
   }
}





