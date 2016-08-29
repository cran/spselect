stepwise.ss <- 
function(y, X.3D, y.name, ss, epsilon, verbose=TRUE) {
   check.ss.step(y, X.3D, y.name, ss, epsilon)

   out <- as.list(NULL)
   aic.v <- as.list(NULL)
   beta <- as.list(NULL)
   r <- as.list(NULL)
   stack.ss <- NULL
   i <- 0
   n <- dim(X.3D)[1]
   p <- dim(X.3D)[2]
   k <- length(ss)
   seq.v <- rep(1:length(ss))
   diff.aic <- 2*epsilon

   null <- lm(y ~ 0)
   if (verbose) {
      print(paste("Results are with no intercept."))
      print(paste("Start: AIC = ", round(extractAIC(null)[2], 4), sep=""))
      print(paste(y.name, " ~ 0", sep=""))
   }

   X.out <- X.3D
   X.in <- matrix(0, nrow(X.out),1)

   while(diff.aic>=epsilon) {
      i <- i + 1
      if (i==1) {
         # Step 1: Find the predictor that has the greatest absolute correlation with the response.
         v1 <- pickvar.step.ss(y, y, X.in, X.out, y.name, seq.v, ss, stack.ss, verbose, i, k)

         # Step 2: Compute the resultant regression model and the residuals.
         out[[i]] <- v1$out
         aic.v[[i]] <- v1$aic.v
         r[[i]] <- v1$r
         beta[[i]] <- v1$beta
         X.in <- v1$X.in
         X.out <- v1$X.out
         seq.v <- v1$seq.v
         stack.ss <- v1$stack.ss


      } else if (i>1) {
         # Step 3: Of the candidate variables, find the predictor that has the greatest absolute correlation with the residuals.
         v2 <- pickvar.step.ss(r[[i-1]], y, X.in, X.out, y.name, seq.v, ss, stack.ss, verbose, i, k)

         # Step 4: Add that predictor to the working design matrix, and compute the resultant regression model and the residuals.
         out[[i]] <- v2$out
         aic.v[[i]] <- v2$aic.v
         r[[i]] <- v2$r
         beta[[i]] <- v2$beta
         X.in <- v2$X.in
         X.out <- v2$X.out
         seq.v <- v2$seq.v
         stack.ss <- v2$stack.ss

         diff.aic <- aic.v[[i-1]][2]-aic.v[[i]][2]
         if (verbose) {
            print(paste("Diff in AIC = ", round(diff.aic, 4), sep=""))
         }
   
         # When diff.aic < epsilon, flag1 is set to TRUE, which tells R to break out of the while loop.
         flag1 <- ifelse(diff.aic<epsilon, TRUE, FALSE) 
   
         # When all predictors have been added, flag2 is set to TRUE, which tells R to break out of the while loop.
         flag2 <- ifelse(i==p, TRUE, FALSE)
         
         # Step 5: Check to see if there is an inadequate improvement in the performance of the model or if all predictors have been added to the model.
         if (flag1==TRUE) {
            final <- out[[(i-1)]]
            stack.ss <- stack.ss[-i]
            summary.final <- summary(out[[(i-1)]])
            final.AIC <- round(AIC(final),2)
            print(final)
   	      print(paste("Diff in AIC = ", round(diff.aic, 4), sep="")) 
            print(paste("Epsilon = ", epsilon, sep=""))           
            print(paste("No. of vars in model = ", dim(X.in)[2]-1, sep=""))
            print(paste("Total no. of vars possible = ", p, sep=""))
            print(paste("Last var we considered when we stopped = ", names(X.in[i]), sep=""))
            break
         } else if (flag2==TRUE) {
            final <- out[[i]]
            stack.ss <- stack.ss
            summary.final <- summary(out[[i]])
            final.AIC <- round(AIC(final),2)
            print(final)
   	      print(paste("Diff in AIC = ", round(diff.aic, 4), sep=""))
            print(paste("Epsilon = ", epsilon, sep=""))           
            print(paste("No. of vars in model = ", dim(X.in)[2], sep=""))
            print(paste("Total no. of vars possible = ", p, sep=""))
            break
         }
      } 
   }
   return(list(beta.final=final, aic.final=final.AIC, summary.final=summary.final, stack.ss=stack.ss))
}
