rm(list = ls())

comp_gmixt <- function(k,x,p,mu,sigma){
  #renvoie la densite du reel x pour la keme classe du melange multiplie par son poids
  return(p[k]*dnorm(x,mu[k],sqrt(sigma[k])))
}

f <- function(x, mu,k){
  #fonction utile dans EM
  return((x-mu[k])**2)
}


est_EM <- function(x,n_em,epsilon,p_0, mu_0, sig_0) {

  # initialisation
  K <- length(mu_0)
  nb_obs <- length(x)
  M<-matrix(0, nrow=nb_obs, ncol=K)
  Y<-matrix(0, nrow=nb_obs, ncol=K)
  p <- matrix(0,n_em+1,K)
  mu <- matrix(0,n_em+1,K)
  sig <- matrix(0,n_em+1,K)
  like <- rep(0,n_em)
  p[1,]<- p_0
  mu[1,] <- mu_0
  sig[1,] <- sig_0
  like[1] <- loglikelyhood_mixtmod(x,p[1,],mu[1,],sqrt(sig[1,]))
  
  # Maximization
  for (t in 1:n_em){
    M <- sapply(1:K,comp_gmixt,x,p[t,],mu[t,],sig[t,])/dmixtmod(x,p[t,],mu[t,],sig[t,])
    N <- colSums(M)
    p[t+1,] <- N/nb_obs
    mu[t+1,] <- (t(M)%*%x)/N
    Y <- sapply(x,f,mu[t+1,],1:K)
    sig[t+1,]<- sqrt(diag(Y%*%M)/N)
    like[t+1] <- loglikelyhood_mixtmod(x,p[t+1,],mu[t+1,],sqrt(sig[t+1,]))
    if (abs(like[t+1]-like[t])<epsilon) {
      #si la vraisemblance ne croit plus "assez" on arrete la fonction
      n_stop <- t 
      break}
  }
  n_stop <- n_em
  return(list(p=p[1:n_stop,], mu=mu[1:n_stop,], sig=sig[1:n_stop,],like=like[1:n_stop]))
}



#em à p et sigma fixé
est_EM_mu <- function(x,n_em,epsilon,p, mu_0, sig) {
  # initialisation
  K <- length(mu_0)
  nb_obs <- length(x)
  M<-matrix(0, nrow=nb_obs, ncol=K)
  Y<-matrix(0, nrow=nb_obs, ncol=K)
  mu <- matrix(0,n_em+1,K)
  mu[1,] <- mu_0
  like <- rep(0,n_em)
  like[1] <- loglikelihood_mm(x,p,mu[1,],sig)
  # Maximization
  for (t in 1:n_em){
    M <- sapply(1:K,comp_gmixt,x,p,mu[t,],sig)/dmixtmod(x,p,mu[t,],sig)
    N <- colSums(M)
    mu[t+1,] <- (t(M)%*%x)/N
    Y <- sapply(x,f,mu[t+1,],1:K)
    like[t+1] <- loglikelihood_mm(x,p,mu[t+1,],sig)
    if (abs(like[t+1]-like[t])<epsilon) {
      n_stop <- t 
      break}
  }
  n_stop <- n_em
  return(list(mu=mu[1:n_stop,],like=like[1:n_stop]))
}
