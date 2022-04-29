rm(list = ls())


gen3_mix <- function(n,p,mu, sigma) {
  z <- sample(c(1:length(p)),n,replace = TRUE ,prob=p)
  return(data.frame(x=rnorm(n,mu[z],sigma[z]),z=z))
}



melange <- gen3_mix(1000,c(0.3,0.7),c(0,2.5),c(1,1))


gibbs <- function(n_gibbs,x,K,gamma,s,l,eps,v,mu0,p0,sig0){
  p <- matrix(0,n_gibbs+1,K)
  mu <- matrix(0,n_gibbs+1,K)
  sig <- matrix(0,n_gibbs+1,K)
  p[1,]<- p0
  mu[1,] <- mu0
  sig[1,] <- sig0
  u <- matrix(0,K,length(x)) #matrice de proba
  like <- rep(0,n_gibbs+1)
  #like[1] <- loglikelihood_mm(x,p[1,],mu[1,],sig[1,])
  z <- rep(0,length(x))
  me <- c(0,0)
  ve <- c(0,0)
  var_b <- c(0,0)
  #moy <- matrix(0,n_gibbs+1,K)
  #vr <- matrix(0,n_gibbs+1,K)
  #nb <- matrix(0,n_gibbs+1,K)
  n <- matrix(0,n_gibbs+1,K)
  n[1,] <- p0*length(x)
  for (t in 1:n_gibbs){
    for (i in 1:length(x)){
      u[,i] <- p[t,]*dnorm(x[i],mu[t,],sqrt(sig[t,]))
      if (sum(u[,i])==0) {
        u[,i]=rep(1/K,K)}
      u[,i] <- u[,i]/sum(u[,i])
      z[i] <- sample(1:K,size = 1,replace=TRUE,prob = u[,i])
    }
    n[t+1,as.numeric(levels(as.factor(z)))] <- table(z)
    p[t+1,] <- rdirichlet(1,gamma+n[t+1,])
    for (j in 1:K){
      me[j] <- ((mean(x[z==j])*n[t+1,j])+(l[j]*eps[j]))/(l[j]+n[t+1,j])
      var_b[j] <- mean(x[z==j]^2)-mean(x[z==j])^2
      ve[j] <- s[j] + n[t+1,j]*var_b[j] + (l[j]*n[t+1,j]*(eps[j]-mean(x[z==j]))^2)/(n[t+1,j]+l[j])
      # if (sum(z==j)==0){
      #   me <- eps[j]
      #   ve <- s[j]
      # }
      #nb[t,j] <- n[j]
      #moy[t,j] <- me
      #vr[t,j] <- ve
    mu[t+1,] <- rnorm(K,me,sqrt(sig[t,]/(n[t+1,]+l)))
    sig[t+1,] <- rinvgamma(K,(v+n[t+1,])/2,ve/2)
    #like[t+1] <- loglikelihood_mm(x,p[t+1,],mu[t+1,],sig[t+1,])
    }
  }
  return(list(p=p,mu=mu,sig=sig,like=like,moy=me,vr=ve,n=n))
  
}



K <- 2
gamma <- rep(1/(K+1),K)
s <- 2*rep(var(melange$x),K) #a definir
l <- rep(1,K) #a definir
eps <- rep(mean(melange$x),K) #a definir
v <- rep(20,K)
p0<- rep(1/K,K)
mu0 <- c(5,1.5)
sig0 <- rep(1,K)
n_gibbs <- 10000
test_ <-gibbs(n_gibbs,melange$x,K,gamma,s,l,eps,v,mu0,p0,sig0)


test_$p[1000:10000,]*rnorm()


dmixtmod(a,p0,mu0,sig0)

dmixtmod <- function(x,p,mu,sig){
  return(colSums(p*sapply(x,dnorm,mu,sqrt(sig))))
}

loglikelihood_mm <-function(x,p,mu,sigma){
  return(sum(log(sapply(x,dmixtmod,p,mu,sigma))))
}



dtest <- function(mu_1, mu_2, x){
  output <- 0
  for(i in 1:length(x)) {
    output <- output + log(.7*dnorm(x[i], mu_1) + 0.3*dnorm(x[i], mu_2)) #proprotions à adapter
  }
  return(output)
}



mu_1 <- seq(-2.5, 4.5, .1)
mu_2 <- seq(-2.5, 4.5, .1)
like <- outer(mu_1, mu_2, dtest, x = melange$x)

par(mfrow = c(1,1))

image(mu_1,mu_2,like,xlab=expression(mu[1]), ylab=expression(mu[2]),col=heat.colors(250))
contour(mu_1,mu_2,like,levels=seq(min(like),max(like),100), add=TRUE,drawlabels=FALSE,legend='gibbs1')
points(test_$mu,pch=".")



a <- seq(-5,7,0.01)
a_ <- seq(0,5,0.01)

nb <- c(0,0)
moy <- c(0,0)
ve <- c(0,0)
var_b <- c(0,0)
for (j in 1:2){
  nb[j] <- sum(melange$z==j)
  moy[j]<- ((mean(melange$x[melange$z==j])*nb[j])+(l[j]*eps[j]))/(l[j]+nb[j])
  var_b[j] <- mean(melange$x[melange$z==j]^2)-mean(melange$x[melange$z==j])^2
  ve[j] <- s[j] + nb[j]*var_b[j] + (l[j]*nb[j]*(eps[j]-mean(melange$x[melange$z==j]))^2)/(nb[j]+l[j])
}



hist(test_$mu[1000:10000,1],breaks=200,freq=FALSE,xlab = "",main="")
lines(a,dnorm(a,test_$mu[1000:10000,1],sqrt(test_$sig[1000:10000,1]/(test_$n[1000:10000,1]+l[1]))),type='l',col='red')

hist(test_$mu[1000:10000,2],breaks=200,freq=FALSE,xlab = "",main="")
lines(a,dnorm(a,moy[2],sqrt(1/(nb[2]+l[2]))),type='l',col='red')

test_$n

his <- c()
for (i in 1:length(test_$p[,1])){
  his <- cbind(his,dmixtmod(a,test_$p[i+998,],test_$mu[i+998,],sqrt(test_$sig[i+998,])))
}

dmixtmod(a,test_$p[998,],test_$mu[998,],sqrt(test_$sig[998,]))

warnings()
hist(his,breaks=200,freq=FALSE,xlab = "",main="")

lines(a,dnorm(a,test_$param[1],sqrt(1/(test_$param[5]+l[1]))),type='l',col='red')

lines(a_,(dinvgamma(a_,(v[1]+nb[1])/2,ve[1]/2)),type='l',col='red')

hist(test_$mu[1000:10000,2],breaks=200,freq=FALSE,xlab = "",main="")






bn <- 1000


hist(gen3_mix(n_gibbs-bn,test_$p[bn:n_gibbs,],test_$mu[bn:n_gibbs,],sqrt(test_$sig[bn:n_gibbs,]))$x,breaks=200,main='',xlab='',freq=FALSE)
lines(a,dmixtmod(a,c(0.3,0.7),c(0,2.5),c(1,1)),col='red',lwd=3)


hist(rnorm(n_gibbs-bn,test_$mu[bn:n_gibbs,1],sqrt(test_$sig[bn:n_gibbs,1]/(test_$n[bn:n_gibbs,1]+l[1]))),breaks=100,freq=FALSE)
lines(a,dnorm(a,moy[2],sqrt(1/(nb[2]+l[2]))),type='l',col='red')


rnorm(2,c(0,2.5),c(1,1))
0.7*2.6836433

test_$mu
hist(test_$sig[1000:10000,1],breaks=200,freq=FALSE,xlab = "",main="")
lines(a_,(dinvgamma(a_,(v[1]+nb[1])/2,ve[1]/2)),type='l',col='red')


hist(test_$sig[1000:10000,2],breaks=200,freq=FALSE,xlab = "",main="")
lines(a_,dinvgamma(a_,(v[2]+nb[2])/2,ve[2]/2),type='l',col='red')




gen3_mix <- function(n,p,mu, sigma) {
  z <- sample(c(1:length(p)),n,replace = TRUE ,prob=p)
  return(data.frame(x=rnorm(n,mu[z],sigma[z]),z=z))
}


melange <- gen3_mix(1000,c(0.3,0.7),c(0,2.5),c(1,1))

rnorm(3,c(0,1,5))

gibbs <- function(n_gibbs,x,K,gamma,s,l,eps,v,mu0,p0,sig0){
  p <- matrix(0,n_gibbs+1,K)
  mu <- matrix(0,n_gibbs+1,K)
  sig <- matrix(0,n_gibbs+1,K)
  p[1,]<- p0
  mu[1,] <- mu0
  sig[1,] <- sig0
  u <- matrix(0,K,length(x)) #matrice de proba
  like <- rep(0,n_gibbs+1)
  #like[1] <- loglikelihood_mm(x,p[1,],mu[1,],sig[1,])
  z <- rep(0,length(x))
  me <- c(0,0)
  ve <- c(0,0)
  var_b <- c(0,0)
  #moy <- matrix(0,n_gibbs+1,K)
  #vr <- matrix(0,n_gibbs+1,K)
  #nb <- matrix(0,n_gibbs+1,K)
  for (t in 1:n_gibbs){
    n <- rep(0,K)
    for (i in 1:length(x)){
      u[,i] <- p[t,]*dnorm(x[i],mu[t,],sqrt(sig[t,]))
      if (sum(u[,i])==0) {
        u[,i]=rep(1/K,K)}
      
      u[,i] <- u[,i]/sum(u[,i])
      z[i] <- sample(1:K,size = 1,replace=TRUE,prob = u[,i])
    }
    n[as.numeric(levels(as.factor(z)))] <- table(z)
    p[t+1,] <- rdirichlet(1,gamma+n)
    for (j in 1:K){
      me[j] <- ((mean(x[z==j])*n[j])+(l[j]*eps[j]))/(l[j]+n[j])
      var_b[j] <- mean(x[z==j]^2)-mean(x[z==j])^2
      ve[j] <- s[j] + n[j]*var_b[j] + (l[j]*n[j]*(eps[j]-mean(x[z==j]))^2)/(n[j]+l[j])
      # if (sum(z==j)==0){
      #   me <- eps[j]
      #   ve <- s[j]
      # }
      #nb[t,j] <- n[j]
      #moy[t,j] <- me
      #vr[t,j] <- ve
      mu[t+1,] <- rnorm(K,me,sqrt(sig[t,]/(n+l)))
      sig[t+1,] <- rinvgamma(K,(v+n)/2,ve/2)
      #like[t+1] <- loglikelihood_mm(x,p[t+1,],mu[t+1,],sig[t+1,])
    }
  }
  return(list(p=p,mu=mu,sig=sig,like=like,z=z))
  
}



K <- 2
gamma <- rep(1/(K+1),K)
s <- 2*rep(var(melange$x),K) #a definir
l <- rep(1,K) #a definir
eps <- rep(mean(melange$x),K) #a definir
v <- rep(20,K)
p0<- rep(1/K,K)
mu0 <- c(1,1.5)
sig0 <- rep(1,K)
n_gibbs <- 10000
test_ <-gibbs(n_gibbs,melange$x,K,gamma,s,l,eps,v,mu0,p0,sig0)

