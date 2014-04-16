

library("mvtnorm")
source("ffbsolver.R")


numTeams <- 10 # the total number of teams
rosterSize <- 10
minMean <- 4

numQB <- 1
numTE <- 1
numWR <- 2
numRB <- 2
numRB_WR <- 1

rosterCompleteThreshold <- c(numQB, numTE, numWR + numRB_WR, numRB + numRB_WR)

posConstraints <- c(numQB, numTE, numWR, numRB, numRB_WR)

replaceNA <- function(col)
{
  mu <- mean(col, na.rm=TRUE)
  col[is.na(col)] <- mu
  col
}

data = read.csv("data06pivot.csv", row.names=1)
positions = read.csv("positions.csv",row.names=1)
positions <- positions[1:4,]
# the last row is the mean -- remove it
data= data[1:(nrow(data)-1),] 
data <- apply(data, 2, replaceNA)

# we will assume we are in week 17 -- the first 16 weeks are used to make decisions, 
# and we can use the data in week 17 to test our performance

history <- data[1:16,]
historyMeans <- colMeans(history)

# Trim the data to ignore players that are not valuable
goodPlayers <- which(historyMeans >= minMean)
goodPlayerData <- history[,goodPlayers]
goodPlayerPos <- positions[,goodPlayers]
numGoodPlayers <- dim(goodPlayerData)[2]
goodPlayerMeans <- colMeans(goodPlayerData)

chooseActivePlayers <- function(means, pos){
  res <- numeric(0)  
  
  # choose QBs
  i <- 1
  while(i <= numQB){  
    curPlayer <- which(means == max(means[which(pos[1,] == 1)]))
    means[curPlayer] <- 0
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  # choose TEs
  i <- 1
  while(i <= numTE){  
    curPlayer <- which(means == max(means[which(pos[2,] == 1)]))
    means[curPlayer] <- 0
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  # choose WRs
  i <- 1
  while(i <= numWR){  
    curPlayer <- which(means == max(means[which(pos[3,] == 1)]))
    means[curPlayer] <- 0
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  # choose RBs
  i <- 1
  while(i <= numRB){  
    curPlayer <- which(means == max(means[which(pos[4,] == 1)]))
    means[curPlayer] <- 0
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  # choose opponent RB_WR
  i <- 1
  while(i <= numRB_WR){  
    curPlayer <- which(means == max(means[c(which(pos[3,] == 1), which(pos[4,] == 1))]))
    means[curPlayer] <- 0
    res <- c(res, curPlayer)
    i <- i + 1
  }
  
  res
}

posCount <- function(roster, pos){
  if(length(roster) == 0){
    c(0, 0, 0, 0)
  }
  else{
    rowSums(sapply(roster, function(i) pos[,i]))
  }
}

isComplete <- function(roster, pos){
  if (length(roster) == 0){ 
    FALSE
  }
  else{
    x <- posCount(roster, pos)
    prod(x >= c(numQB, numTE, numWR, numRB)) && (x[3] + x[4] >=numWR + numRB + numRB_WR)
  }
}

draftPlayer <- function(roster, means, pos){ 
  
  posMeans <- as.matrix(pos) %*% as.matrix(means) /rowSums(pos)
  playerValues <- means - (t(posMeans) %*% as.matrix(pos))
  playerValues[which(playerValues< 0)] <- 0
  
  player <- sample(1 : length(playerValues), 1, FALSE, playerValues)
  posIndex <- which(pos[,player]==1)
  posFull <- posCount(roster, pos)[posIndex] >= rosterCompleteThreshold[posIndex]
  rosterIncomplete <- !isComplete(roster, pos)
  while(rosterIncomplete && posFull){
    player <- sample(1 : length(playerValues), 1, FALSE, playerValues)
    posIndex <- which(pos[,player]==1) 
    posFull <- posCount(roster, pos)[posIndex] > rosterCompleteThreshold[posIndex]
  }
  
  player
}

# the repeatable simulation starts here
wins <- matrix(ncol=3, nrow=0)
numSims <- 2
simIter <- 1
while(simIter <= numSims){

# We assume that each manager has a roster -- so allocate some of the players to the rosters
# to make the simulation realistic, we will have the managers draft players

availableMeans <- goodPlayerMeans
rosters <- list()
i <- 1
while(i <= numTeams){
  rosters[[i]] <- numeric(0)
  i <- i + 1
}
startTeam <- sample(0 : numTeams-1, 1)
round <- 1
dir <- 1
while(round <= rosterSize){ 
  i <- 0
  while(i <= numTeams - 1){
    curTeam <- (startTeam + dir * i) %% numTeams + 1
    # the manager drafts a player
    player <- draftPlayer(rosters[[curTeam]], availableMeans, goodPlayerPos)
    # to prevent that player from being drafted again, we'll just make the value 0
    availableMeans[player] <- 0
    rosters[[curTeam]] <- c(rosters[[curTeam]], player)
    i <- i + 1
  }
  dir <- dir * -1
  round <- round + 1
}


#rosters[[1]] is my roster, rosters[[2]] is the opponent's roster, all other players are irrelevant
myPlayerData <- goodPlayerData[, rosters[[1]]]
myPlayerPos <- goodPlayerPos[, rosters[[1]]]
myPlayerMeans <- colMeans(myPlayerData)

oppPlayerData <- goodPlayerData[, rosters[[2]]]
oppPlayerPos <- goodPlayerPos[, rosters[[2]]]
oppPlayerMeans <- colMeans(oppPlayerData)

#build the dataset of players that are available to me and to the opponent
selectedPlayers <- numeric(0)
i <- 1
while(i <= numTeams){
  selectedPlayers <- c(selectedPlayers, rosters[[i]])
  i <- i + 1
}
availablePlayerData <- goodPlayerData[, -selectedPlayers]
availablePlayerPos <- goodPlayerPos[, -selectedPlayers]
availablePlayerMeans <- colMeans(availablePlayerData)

allPlayerPos <- cBind(oppPlayerPos, myPlayerPos, availablePlayerPos)

#choose opponent's active players
oppActive <- chooseActivePlayers(c(oppPlayerMeans, rep(0, length(myPlayerMeans)), availablePlayerMeans), allPlayerPos)
# randomly choose the players that remain on the opponent's roster
oppRoster <- oppActive
keepablePlayers <- setdiff(seq(1, rosterSize), oppActive[which(oppActive <= rosterSize)])
keptPlayers <- sample(keepablePlayers, rosterSize - length(oppRoster))
oppRoster <- c(oppRoster, keptPlayers)

# choose my active players by maximizing expected points
mySelectionPlayerMeans <- c(oppPlayerMeans, myPlayerMeans, availablePlayerMeans)
mySelectionPlayerMeans[oppRoster] <- 0
myActive <- chooseActivePlayers(mySelectionPlayerMeans, allPlayerPos)

# put everything back together and start building the data for the optimization routine
optPlayerData <- cbind(oppPlayerData, myPlayerData, availablePlayerData)
numPlayers <- dim(optPlayerData)[2]
mu <- colMeans(optPlayerData)
sigma <- cov(optPlayerData)
# sigma may only be positive semidefinite, and we want it to be positive definite
sigma <- as.matrix(nearPD(sigma)$mat)


# run the solver and get a set of solutions
solns <- ffbSolve(mu, sigma, allPlayerPos, oppRoster, oppActive, posConstraints)

# see how well we did
getOutcome <- function(){
  rmvnorm(1, mean=mu, sigma=sigma)
}

outcome <- getOutcome()

solns <- cbind(sort(myActive), solns)
myPoints <- apply(solns, 2, function(x) sum(outcome[x]))
oppPoints <- sum(outcome[oppActive])

win <- myPoints > oppPoints

wins <- rbind(wins, win)
simIter <- simIter+1

}

  