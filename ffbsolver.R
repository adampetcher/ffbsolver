
# A solver for selecting an optimal lineup in fantasy football.  The main function is ffbsolve()
# parameters:
## mu -- a vector of expected point totals for each player
## sigma -- a covariance matrix describing the points of all players
## pos -- a position matrix describing the position of each player
## oppRoster -- the roster of the opponent
## oppActive -- the active players of the opponent
## posConstraints -- a vector descrbing the maximum number of players for each position (QB, TE, WR, RB, WR or RB)
## initSoln -- an initial "good" integer solution -- the solver will try to find a better one

library("Matrix")
library("DEoptim")

epsilon <- 10^-8

objective_approx <- function(m, sigma, obj_mu, o_sigma, o_mu){
  m_mat <- as(m, "matrix")
  i <- 1
  j <- 1
  minm <- matrix(0, nrow=length(m), ncol=length(m))
  while(i <= length(m)){
    while(j <= length(m)){
      minm[i,j] = min(m[i], m[j])
      j <- j + 1
    }
    i <- i + 1
  }
  
  covSum <- o_sigma+sum(sigma*minm)
  (o_mu - t(m_mat)%*%obj_mu)/sqrt(covSum)
}

objective <- function(m, sigma, obj_mu, o_sigma, o_mu){
  m_mat <- as(m, "matrix")
  covSum <- o_sigma+t(m_mat)%*%sigma%*%m_mat
  (o_mu - t(m_mat)%*%obj_mu)/sqrt(covSum)
}

objective_with_constraints <- function(m, sigma, obj_mu, o_sigma, o_mu, A, b, penalty=10){
  m_mat <- as(m, "matrix")
  covSum <- o_sigma + t(m_mat)%*%sigma%*%m_mat
  res <- (o_mu - t(m_mat)%*%obj_mu)/sqrt(covSum)
  feasibleTest <- A %*% m - b 
  if(prod(feasibleTest >= 0)){
    res
    }
  else 
  {
    feasibleTest <- pmin(feasibleTest, 0)
    res + 2 ^ penalty*abs(sum(feasibleTest))
  }
}

# gradient of objective function
gradient <- function(m, sigma, obj_mu, o_sigma, o_mu){
  m_mat <- as(m, "matrix")
  covSum <- as.vector(o_sigma+t(m_mat)%*%sigma%*%m_mat)
  term1 <- -obj_mu/sqrt(covSum)
  term2 <- -(as.vector(-t(m_mat)%*%obj_mu+o_mu)*(sigma%*%m_mat))/covSum^(3/2)
  term1 + term2
}

chooseActivePlayersRandom <- function(dontSelect, pos, posConstraints){
  res <- numeric(0)  
  
  curPos <- 1
  while(curPos <= 4){
    i <- 1
    while(i <= posConstraints[curPos]){  
      options <- setdiff(which(pos[curPos,] == 1), dontSelect)
      curPlayer <- sample(c(options, options), 1)
      dontSelect <- c(dontSelect, curPlayer)
      res <- c(res, curPlayer)
      i <- i + 1
    }
    curPos <- curPos + 1
  }
  
  # choose opponent RB_WR
  i <- 1
  while(i <= posConstraints[5]){  
    options <- setdiff(c(which(pos[3,] == 1), which(pos[4,] == 1)), dontSelect)
    curPlayer <- sample(c(options, options), 1)
    dontSelect <- c(dontSelect, curPlayer)
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  res
}

integerizeSolution <- function(soln, pos, posConstraints) {
  res <- numeric(0)
  curPos <- 1
  while(curPos <= 4){
    i <- 1
    while(i <= posConstraints[curPos]){
      newPlayer <- sample(1:length(soln), 1, prob=soln * pos[curPos,])
      soln[newPlayer] <- 0
      res <- c(res, newPlayer)
      i <- i + 1
    } 
    curPos <- curPos + 1
  }
  
  i <- 1
  while(i <= posConstraints[5]){
    newPlayer <- sample(1:length(soln), 1, prob=soln * (pos[3,] + pos[4,]))
    soln[newPlayer] <- 0
    res <- c(res, newPlayer)
    i <- i + 1
  }
  
  sort(res)
}

ffbSolve_constrOptim <- function(options, m_init, A, b, sigma, mu, o_sigma, o_mu){
  numPlayers <- length(m_init)
  m_init_interior <- m_init + rep(epsilon/numPlayers, numPlayers) - m_init * rep(epsilon)
  res <- constrOptim(m_init_interior, objective, gradient, A, b, outer.iterations=1000, 
                     sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu)
  res$par
}

ffbSolve_constrOptim_approx <- function(options, m_init, A, b, sigma, mu, o_sigma, o_mu){
  numPlayers <- length(m_init)
  m_init_interior <- m_init + rep(epsilon/numPlayers, numPlayers) - m_init * rep(epsilon)
  res <- constrOptim(m_init_interior, objective_approx, NULL, A, b, outer.iterations=1000, 
                     sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu)
  res$par
}

ffbSolve_DEOptim <- function(options, m_init, A, b, sigma, mu, o_sigma, o_mu){
  randomInitPop <- options$randomInitPop
  penalty <- options$penalty
  numPlayers <- length(m_init)
  initPop <- NULL
  popSize <- 10 * numPlayers
  if(randomInitPop){
    oppRoster <- options$oppRoster
    pos <- options$pos
    posConstraints <- options$posConstraints
    
    initPop <- t(m_init)
    dontSelect <- oppRoster
    while(dim(initPop)[1] < popSize){
      newActive <- chooseActivePlayersRandom(dontSelect, pos, posConstraints)
      newMem <- as.vector(sparseVector(i=sort(newActive), x=1, length=numPlayers)) 
      initPop <- rBind(initPop, t(newMem))
    }
  }
  
  res <- DEoptim(objective_with_constraints, rep(0, numPlayers), rep(1, numPlayers), 
                 sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu, A=A, b=b, penalty=penalty,
                 DEoptim.control(iter=2000, trace=FALSE,
                                 NP=popSize,initialpop=initPop, 
                                 parallelType=1, parVar=c()))
  soln <- res$optim$bestmem
  # there is a chance that the solution will be invalid -- use the default solution in this case
  #if(prod(A %*% soln - b >= 0) == 0){
  #  soln <- m_init
  #}
  soln
}

ffbOptMethod <- function(options, m_init, A, b, sigma, mu, o_sigma, o_mu){
  options$fn(options, m_init, A, b, sigma, mu, o_sigma, o_mu)
}

ffbSolve <- function(mu, sigma, pos, oppRoster, oppActive, posConstraints, initSoln, methods, fixedPlayers=c()){
  numPlayers <- length(mu)
  rosterSize <- length(oppRoster)
  
  # opponents roster as an indicator vector
  o <- as.vector(sparseVector(i=sort(as.vector(oppActive)), x=1, length=numPlayers)) 
  # we pre-compute some parts of the objective function up front
  o_mat <- as(o, "matrix")
  o_sigma <- t(o_mat)%*%sigma%*%o_mat
  o_mu <- t(o_mat)%*%mu
  
  # constraints
  A_id <- diag(1, numPlayers) # at least 0 and at most 1 of any player
  
  # position constraints
  A_pos <- -1*pos
  A_pos <- rBind(A_pos, -(pos[3,] + pos[4,]))
  A_pos <- as.matrix(A_pos)  
  
  # don't choose any player on the opponent's roster
  A_unavailable <- -1* as.matrix(sparseMatrix(seq(1, rosterSize), oppRoster, dims=c(rosterSize, numPlayers)))
  
  A <- rbind(A_id, -1 * A_id, A_pos)
  A <- rbind(A, A_unavailable)
  
  # fixed players must be active
  i <- 1
  while(i <= length(fixedPlayers)){
    A <- rbind(A, as.vector(sparseVector(length=numPlayers, i=fixedPlayers[i], x=1)))
    i <- i + 1
  }
  
  b <- c(rep(0, numPlayers), rep(-1, numPlayers),-posConstraints[1],-posConstraints[2],
         -(posConstraints[3] + posConstraints[5]),-(posConstraints[4] + posConstraints[5]),
         -(posConstraints[3] + posConstraints[4] + posConstraints[5]))
  
  # don't select players from the opponent's roster
  b <- c(b, rep(-epsilon, rosterSize))
  
  # always select fixed players
  i <- 1
  while(i <= length(fixedPlayers)){
    b <- c(b, 1 - epsilon)
    i <- i + 1
  }
  
  m_init <- as.vector(sparseVector(i=sort(as.vector(initSoln)), x=1, length=numPlayers)) 

  # methods is a list of optimization methods and parameters.  
  # For each one, compute the fractional solution
  
  solns <- lapply(methods, ffbOptMethod, m_init, A, b, sigma, mu, o_sigma, o_mu) 
  
  solnsVals <- sapply(solns, objective, sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu)
  initVal <- objective(m_init, sigma, mu, o_sigma, o_mu)
  
  integerSolns <- lapply(solns, integerizeSolution, pos=pos, posConstraints=posConstraints)
  integerSolnsVecs <- lapply(integerSolns, sparseVector, x=1, length=length(mu))
  integerSolnsVecs <- lapply(integerSolnsVecs, as.vector)
  intSolnsVals <- sapply(integerSolnsVecs, objective, sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu)
  
  solnsValid <- sapply(solns, function(x) prod(A %*% x - b >= 0))
  intSolnsValid <- sapply(integerSolnsVecs, function(x) prod(A %*% x - b >= 0))
  
  list(solns=solns, integerSolns=integerSolns, initVal=initVal, solnsVals=solnsVals, intSolnsVals=intSolnsVals,
       solnsValid=solnsValid, intSolnsValid=intSolnsValid)
  
}





  