
source("ffbsolver.R")

posConstraints <- c(1, 1, 1, 1, 1)
mu <- c(rep(10, 9), 9, 8)
sigma <- matrix(0, nrow=11, ncol=11)
sigma[11,11] <- 50
sigma <- as.matrix(nearPD(sigma)$mat)

pos <- cbind(diag(4), c(0, 0, 0, 1), diag(4), c(0, 0, 0, 1), c(0, 0, 0, 1))

oppRoster <- 1:5
oppActive <- 1:5

initSoln <- 6:10
solns <- ffbSolve(mu, sigma, pos, oppRoster, oppActive, posConstraints, initSoln, globalOptim=TRUE)

# see how well we did
getOutcome <- function(){
  rmvnorm(1, mean=mu, sigma=sigma)
}

outcome <- getOutcome()

solns <- cbind(sort(initSoln), solns)
myPoints <- apply(solns, 2, function(x) sum(outcome[x]))
oppPoints <- sum(outcome[oppActive])
