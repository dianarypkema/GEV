# script create Markov chain for environmental states

randcomp <- function(r1,x1) {
  #takes a number 0 < r1 < 1 and a cum. discrete prob. distr. x1
  # finds the spot in x1 such that x1(loc-1) <= r1 < x1(loc)

  for(i in 1:length(x1)) {
    if(r1 < x1[i]) {
      loc <- i
      break
    }
  }
  return(loc)}

mchain <- function(p,time1, xx) {
  #REVISED March 20, 2018
  #procedure to simulate from a Markov chain
  #generates sample path, equilibrium dist, sequence of states
  #inputs are (1) transition probability matrix,
  #p(i,j) defined as Prob (i-> j);
  #length of run, T;
  #calls function randcomp;
  #OUTPUT is a vector of integers of length T+1, includes
  #starting state chosen at random according to stationary

  #first make sure that we have Prob(i -->j)
  #which means that rowsums are 1
  rowSumsp <- as.numeric(as.character(signif(sum(rowSums(p))), 6))
  if(rowSumsp != as.numeric(as.character(nrow(p)))){
    print("matrix rows do NOT sum to 1")
    break
  }else{

    #find eqm distribution pvec
    vec <- eigen(p)$vectors
    c <- diag(eigen(p)$values)
    cmax <- max(abs(diag(c)))
    pvec <- vec[,1]
    pvec <- pvec/sum(pvec)

    #start by setting up cumulatve probability distribution for equilibrium vector
    cpvec <- cumsum(pvec)
    cpvec <- Re(cpvec)
    #and for all rows
    cp <- t(apply(p, 1, cumsum)) # this does cumulative sum along rows

    #initialize output vector
    sapath <- rep(0, time1+1);
    #find starting state at time 0
    sapath[1] = randcomp(xx[1],cpvec);
    i0 = sapath[1];
    for(i in 1:time1) {
      i0 = randcomp(xx[i+1],cp[i0,])
      sapath[i+1] = i0
    }
    return(sapath)
  }
}