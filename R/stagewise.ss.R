stagewise.ss <- 
function(y, X, X.3D, ss, increment, tolerance, col.plot, verbose=TRUE, plot=TRUE) {
   check.ss.stage(y, X, X.3D, ss, increment, tolerance, col.plot)   

   i <- 0
   p <- dim(X.3D)[2]
   k <- length(ss)
   seq.v <- rep(1:length(ss))
   names.X <- dimnames(X.3D)[[2]]
   flag.stop <- FALSE
   stack.ss <- rep(NA,p)

   X.cand <- X.3D
 
   # Step 1: Initialize all regression coefficient estimates equal to 0, and let r=y-ybar.
   r <- y # Note: ybar=0 since all variables were standardized to have mean=0 and SD=1.
   beta.old <- array(0, dim=c(1,p))

   while(flag.stop!=TRUE) {
      i <- i + 1

      # Steps 2-4
      v1 <- pickvar.stage.ss(r, X.cand, seq.v, ss, beta.old, i, p, k, names.X, stack.ss, increment, tolerance, verbose)
      r <- v1$r
      beta.new <- v1$beta.new
      names.X <- v1$names.X
      beta.old <- rbind(beta.old, beta.new)
      colnames(beta.old) <- names.X
      X.cand <- v1$X.cand
      stack.ss <- v1$stack.ss
      flag.stop <- v1$flag
      seq.v <- v1$seq.v
     
      # Step 5: Repeat steps 2-4 until none of the predictors are correlated with the residuals.
      if (flag.stop==TRUE) {
         beta.nonzero <- beta.old[,which(!beta.old[(i+1),] == 0)]
         stack.ss <- na.omit(stack.ss)
         print(beta.nonzero[(i+1), ])
         print(paste("No. of vars in model = ", dim(beta.nonzero)[2], sep=""))
         print(paste("Total no. of vars possible = ", p, sep=""))
         if (plot==TRUE) {
            varnames.X <- names(X)
            names.beta.nonzero <- colnames(beta.nonzero)
            path.index <- which(!is.na(match(varnames.X, names.beta.nonzero)))
            #print(path.index)
            path.plot.stage.ss(beta.nonzero, i, stack.ss, increment, tolerance, path.index, col.plot)
         }
         break
      }
   }
return(list(beta.final=beta.nonzero[(i+1), ], stack.ss=stack.ss))
}
