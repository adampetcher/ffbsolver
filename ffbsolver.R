
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

objective <- function(m, sigma, obj_mu, o_sigma, o_mu){
  m_mat <- as(m, "matrix")
  covSum <- o_sigma+t(m_mat)%*%sigma%*%m_mat
  (-t(m_mat)%*%obj_mu + o_mu)/sqrt(covSum)
}

objective_with_constraints <- function(m, sigma, obj_mu, o_sigma, o_mu, A, b){
  m_mat <- as(m, "matrix")
  covSum <- o_sigma + t(m_mat)%*%sigma%*%m_mat
  res <- (-t(m_mat)%*%obj_mu + o_mu)/sqrt(covSum)
  if(prod(A %*% m - b >= 0)) res else 
    res + abs(min(A %*% m - b))
}

# gradient of objective function
gradient <- function(m, sigma, obj_mu, o_sigma, o_mu){
  m_mat <- as(m, "matrix")
  covSum <- as.vector(o_sigma+t(m_mat)%*%sigma%*%m_mat)
  term1 <- -mu/sqrt(covSum)
  term2 <- -(as.vector(-t(m_mat)%*%obj_mu+o_mu)*(sigma%*%m_mat))/covSum^(3/2)
  term1 + term2
}

chooseActivePlayersRandom <- function(dontSelect, pos, posConstraints){
  res <- numeric(0)  
  
  curPos <- 1
  while(curPos <= 4){
    i <- 1
    while(i <= posConstraints[curPos]){  
      curPlayer <- sample(setdiff(which(pos[curPos,] == 1), dontSelect), 1)
      dontSelect <- c(dontSelect, curPlayer)
      res <- c(res, curPlayer)
      i <- i + 1
    }
    curPos <- curPos + 1
  }
  
  # choose opponent RB_WR
  i <- 1
  while(i <= posConstraints[5]){  
    curPlayer <- sample(setdiff(c(which(pos[3,] == 1), which(pos[4,] == 1)), dontSelect), 1)
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

ffbSolve <- function(mu, sigma, pos, oppRoster, oppActive, posConstraints){
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
  A_pos <- -1*allPlayerPos
  A_pos <- rBind(A_pos, -(pos[3,] + pos[4,]))
  A_pos <- as.matrix(A_pos)  
  
  # don't choose any player on the opponent's roster
  A_unavailable <- -1* as.matrix(sparseMatrix(seq(1, rosterSize), oppRoster, dims=c(rosterSize, numPlayers)))
  
  A <- rbind(A_id, -1 * A_id, A_pos)
  A <- rbind(A, A_unavailable)
  
  b <- c(rep(0, numPlayers), rep(-1, numPlayers),-posConstraints[1],-posConstraints[2],
         -(posConstraints[3] + posConstraints[5]),-(posConstraints[4] + posConstraints[5]),
         -(posConstraints[3] + posConstraints[4] + posConstraints[5]))
  
  # we allow the selection of only a small amount of some player that is on the opponent's roster
  epsilon <- 10^-8
  b <- c(b, rep(-epsilon, rosterSize))
  
  m_init <- as.vector(sparseVector(i=sort(as.vector(myActive)), x=1, length=numPlayers)) 
  
  # local optimization
  # perturb the initial point so it is in the interior of the region
  m_init_interior <- m_init - numPlayers*epsilon*m_init + rep(epsilon/2, numPlayers)
  res <- constrOptim(m_init_interior, objective, gradient, A, b, outer.iterations=1000, 
                     sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu)
  soln_local <- res$par
  
  # global optimization using differential evolution
  initPop <- t(m_init)
  popSize <- 10 * numPlayers
  dontSelect <- which(mySelectionPlayerMeans == 0)
  while(dim(initPop)[1] < popSize){
    newActive <- chooseActivePlayersRandom(dontSelect, pos, posConstraints)
    newMem <- as.vector(sparseVector(i=sort(newActive), x=1, length=numPlayers)) 
    initPop <- rBind(initPop, t(newMem))
  }
  
  
  startTime <- proc.time()
  res <- DEoptim(objective_with_constraints, rep(0, numPlayers), rep(1, numPlayers), 
                 sigma=sigma, obj_mu=mu, o_sigma=o_sigma, o_mu=o_mu, A=A, b=b,
                 DEoptim.control(NP=popSize,initialpop=initPop, parallelType=1, parVar=c()))
  runningTime <- proc.time() - startTime
  soln_global <- res$optim$bestmem
  
  solns <- list(soln_local, soln_global)
  
  integerSolns <- sapply(solns, integerizeSolution, pos=pos, posConstraints=posConstraints)
  
  integerSolns
  
}





  