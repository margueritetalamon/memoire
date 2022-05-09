
#PLOT POUR EM 

plot_contour <- function(melange,mu_reel){
    mu_1 <- seq(-1.5, 4, .1)
    mu_2 <- seq(-1.5, 4, .1)
    like <- outer(mu_1, mu_2, loglikelyhood_psig_fix, x = melange$x)
    image(mu_1,mu_2,like,xlab=expression(mu[1]), ylab=expression(mu[2]),col=heat.colors(150))
    contour(mu_1,mu_2,like,levels=seq(min(like),max(like),100), add=TRUE,drawlabels=FALSE)
    points(mu_reel[1],mu_reel[2],pch=5)
}















#PLOT POUR GIBBS


plot_hist_densite <- function(melange, ech_gibbs, n_gibbs, sequence, burn_in){

#plot l'histogramme de l'echantillon, trace en bleu clair les densité estimées grace à chaque itérations de gibbs à partir de la burn in period, en bleue foncé desnité grace à moyenne a posteriori
  
  hist(melange$x,breaks=20,freq=FALSE,main='',xlab='')
  for (i in bn:n_gibbs){
    lines(sequence,dmixtmod(sequence,ech_gibbs$p[i,],ech_gibbs$mu[i,],ech_gibbs$sig[i,]),col='lightcyan3')
  }
  lines(a,dmixtmod(a,c(mean(ech_gibbs$p[bn:n_gibbs,1]),mean(ech_gibbs$p[bn:n_gibbs,2])),c(mean(ech_gibbs$mu[bn:n_gibbs,1]),mean(ech_gibbs$mu[bn:n_gibbs,2])),c(mean(ech_gibbs$sig[bn:n_gibbs,1]),mean(ech_gibbs$sig[bn:n_gibbs,2]))),lwd=3,col='steelblue')
}



plot_hist_mu_marg <- function(sim_gibbs){
  hist(sim_gibbs$mu[,1],breaks=20,freq=FALSE,main='',xlab='',col='lightcyan3')
  hist(sim_gibbs$mu[,2],breaks=200,freq=FALSE,main='',xlab='',col='darkseagreen3')  
}


plot_traj_mu <- function(estimation){
  plot(estimation$mu[,1],type='l',col='royalblue4',main='',xlab='nombre d\'itérations',ylab=expression(mu[j]),lwd=2)
  lines(estimation$mu[,2],type='l',col='springgreen3')
  legend("bottomleft",legend=c(expression(mu[1]),expression(mu[2])),lty =c(1,1),col=c("royalblue4","springgreen3"),lwd=c(2,2))
}

plot_traj_p <- function(estimation){
  plot(estimation$p[,1],type='l',col='royalblue4',main='',xlab='nombre d\'itérations',ylab=expression(p[j]),lwd=2)
  lines(estimation$p[,2],type='l',col='springgreen3')
  legend("bottomleft",legend=c(expression(p[1]),expression(p[2])),lty =c(1,1),col=c("royalblue4","springgreen3"),lwd=c(2,2))
}

plot_traj_sig <- function(estimation){
  plot(estimation$sig[,1],type='l',col='royalblue4',main='',xlab='nombre d\'itérations',ylab=expression(sig[j]),lwd=2)
  lines(estimation$sig[,2],type='l',col='springgreen3')
  legend("bottomleft",legend=c(expression(sig[1]),expression(sig[2])),lty =c(1,1),col=c("royalblue4","springgreen3"),lwd=c(2,2))
}

