## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------
BiocStyle::latex(width=78, use.unsrturl=FALSE)

## ----echo=FALSE-------------------------------------------------------------
library(ccrepe)

## ----eval=FALSE-------------------------------------------------------------
#  sim.score.args = list(method="spearman", use="complete.obs")

## ----eval=FALSE-------------------------------------------------------------
#  ccrepe(
#   x = NA,
#   y = NA,
#   sim.score = cor,
#   sim.score.args = list(),
#   min.subj = 20,
#   iterations = 1000,
#   subset.cols.x = NULL,
#   subset.cols.y = NULL,
#   errthresh  = 1e-04,
#   verbose = FALSE,
#   iterations.gap = 100,
#   distributions = NA,
#   compare.within.x = TRUE,
#   concurrent.output = NA,
#   make.output.table = FALSE)

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(c(
 "Sample 1", "Sample 2","Sample 3","Sample 4","Sample 5",
 "Sample 6","Sample 7","Sample 8","Sample 9","Sample 10"),
 c("Feature 1", "Feature 2", "Feature 3","Feature 4"))

test.output <- ccrepe(x=test.input, iterations=20, min.subj=10)

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2.  In this case we would expect feature 1 and feature 2 to be associated. In the output we see this by the positive sim.score value in the [1,2] element of test.output\\$sim.score and the small q-value in the [1,2] element of test.output\\$q.values.", fig.width=7, fig.height=3.5,fig.pos="H"----
par(mfrow=c(1,2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm


data2 <- matrix(rlnorm(105,meanlog=0,sdlog=1),nrow=15,ncol=7)
aligned.rows <- c(seq(1,4),seq(6,9),11,12)  # The datasets dont need 
                                            # to have subjects line up exactly
data2[aligned.rows,1] <-  2*data[,3] + rnorm(10,0,0.01)
data2.rowsum <- apply(data2,1,sum)
data2.norm <- data2/data2.rowsum
apply(data2.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input.2 <- data2.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))
dimnames(test.input.2) <- list(paste("Sample",c(seq(1,4),11,seq(5,8),12,9,10,13,14,15)),paste("Feature",seq(1,7)))

test.output.two.datasets <- ccrepe(x=test.input, y=test.input.2, iterations=20, min.subj=10)

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2.  In this case we would expect feature 1 and feature 2 to be associated. In the output we see this by the positive sim.score value in the [1,2] element of test.output\\$sim.score and the small q-value in the [1,2] element of test.output\\$q.values.", fig.width=7, fig.height=3.5, fig.pos="H"----
par(mfrow=c(1,2))
plot(data2[aligned.rows,1],data[,3],xlab="dataset 2: Feature 1",ylab="dataset 1: Feature 3",main="Non-normalized")
plot(data2.norm[aligned.rows,1],data.norm[,3],xlab="dataset 2: Feature 1",ylab="dataset 1: Feature 3",
     main="Normalized")
test.output.two.datasets

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

test.output.nc.score     <- ccrepe(x=test.input, sim.score=nc.score, iterations=20, min.subj=10)

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2.  In this case we would expect feature 1 and feature 2 to be associated. In the output we see this by the positive sim.score value in the [1,2] element of test.output\\$sim.score and the small q-value in the [1,2] element of test.output\\$q.values. In this case, however, the sim.score represents the NC-Score between two features rather than the Spearman correlation.", fig.width=7, fig.height=3.5, fig.pos="H"----
par(mfrow=c(1,2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output.nc.score

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

my.test.sim.score <- function(x,y=NA,constant=0.5){
       if(is.vector(x) && is.vector(y)) return(constant)
          if(is.matrix(x) && is.na(y)) return(matrix(rep(constant,ncol(x)^2),ncol=ncol(x)))
          if(is.data.frame(x) && is.na(y)) return(matrix(rep(constant,ncol(x)^2),ncol=ncol(x)))
          else stop('ERROR')
   }

test.output.sim.score    <- ccrepe(x=test.input, sim.score=my.test.sim.score, iterations=20, min.subj=10, sim.score.args = list(constant = 0.6))

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2.  In this case we would expect feature 1 and feature 2 to be associated. Note that the values of sim.score are all 0.6 and none of the p-values are very small because of the arbitrary definition of the similarity score.", fig.width=7, fig.height=3.5, fig.pos="H"----
par(mfrow=c(1,2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output.sim.score

## ----<----------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data.rowsum <- apply(data,1,sum)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

test.output.1.3     <- ccrepe(x=test.input, iterations=20, min.subj=10, subset.cols.x=c(1,3))
test.output.1       <- ccrepe(x=test.input, iterations=20, min.subj=10, subset.cols.x=c(1), compare.within.x=FALSE)
test.output.12.3    <- ccrepe(x=test.input, iterations=20, min.subj=10, subset.cols.x=c(1,2),subset.cols.y=c(3), compare.within.x=FALSE)
test.output.1.3$sim.score
test.output.1$sim.score
test.output.12.3$sim.score

## ----eval=FALSE-------------------------------------------------------------
#  nc.score(
#   x,
#   y = NULL,
#   use = "everything",
#   nbins = NULL,
#   bin.cutoffs=NULL)

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data.rowsum <- apply(data,1,sum)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

test.output.matrix <- nc.score(x=test.input)
test.output.num    <- nc.score(x=test.input[,1],y=test.input[,2])

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2 of the second example.  Again, we expect to observe a positive association between feature 1 and feature 2.  In terms of generalized checkerboard scores, we would expect to see more co-variation patterns than co-exclusion patterns.  This is shown by the positive and relatively high value of the [1,2] element of test.output.matrix (which is identical to test.output.num)", fig.height=3, fig.pos="H"----
par(mfrow=c(1, 2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output.matrix
test.output.num

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data.rowsum <- apply(data,1,sum)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

test.output <- nc.score(x=test.input,nbins=4)

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2 of the second example.  Again, we expect to observe a positive association between feature 1 and feature 2.  In terms of generalized checkerboard scores, we would expect to see more co-variation patterns than co-exclusion patterns.  This is shown by the positive and relatively high value in the [1,2] element of test.output.  In this case, the smaller bin number yields a smaller NC-score because of the coarser partitioning of the data.", fig.height=3, fig.pos="H"----
par(mfrow=c(1, 2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output

## ---------------------------------------------------------------------------
data <- matrix(rlnorm(40,meanlog=0,sdlog=1),nrow=10,ncol=4)
data.rowsum <- apply(data,1,sum)
data[,1] = 2*data[,2] + rnorm(10,0,0.01)
data.norm <- data/data.rowsum
apply(data.norm,1,sum)  # The rows sum to 1, so the data are normalized
test.input <- data.norm

dimnames(test.input) <- list(paste("Sample",seq(1,10)),paste("Feature",seq(1,4)))

test.output <- nc.score(x=test.input,bin.cutoffs=c(0.1,0.2,0.3))

## ----fig.cap="Non-normalized and normalized associations between feature 1 and feature 2 of the second example.  Again, we expect to observe a positive association between feature 1 and feature 2.  In terms of generalized checkerboard scores, we would expect to see more co-variation patterns than co-exclusion patterns.  This is shown by the positive and relatively high value in the [1,2] element of test.output.  The bin edges specified here represent almost absent ([  0,0.001)), low abundance ([0.001,0.1)), medium abundance ([0.1,0.25)), and high abundance ([0.6,1)).", fig.height=3, fig.pos="H"----
par(mfrow=c(1, 2))
plot(data[,1],data[,2],xlab="Feature 1",ylab="Feature 2",main="Non-normalized")
plot(data.norm[,1],data.norm[,2],xlab="Feature 1",ylab="Feature 2",
     main="Normalized")
test.output

