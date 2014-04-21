
source("ffbsolver.R")

posConstraints <- c(1, 1, 1, 1, 1)
mu <- c(rep(10, 9), 9, 8)
sigma <- diag(11)
sigma[11, 11] <- 24
sigma <- as.matrix(nearPD(sigma)$mat)

pos <- cbind(diag(4), c(0, 0, 0, 1), diag(4), c(0, 0, 0, 1), c(0, 0, 0, 1))

oppRoster <- 1:5
oppActive <- 1:5

initSoln <- 6:10
methods <- list(list(fn=ffbSolve_constrOptim)
                ,list(fn=ffbSolve_constrOptim_approx)
                ,list(fn=ffbSolve_DEOptim, randomInitPop=FALSE)
                )
solns <- ffbSolve(mu, sigma, pos, oppRoster, oppActive, posConstraints, initSoln, methods)

# see how well we did
getOutcome <- function(){
  rmvnorm(1, mean=mu, sigma=sigma)
}

outcome <- getOutcome()

intSolns <- sapply(solns$integerSolns, cbind)
intSolns <- cbind(sort(initSoln), intSolns)
myPoints <- apply(intSolns, 2, function(x) sum(outcome[x]))
oppPoints <- sum(outcome[oppActive])
