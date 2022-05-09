




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


plot_traj <- function(sim_gibbs,param){
  plot(sim_gibbs$param[,1],type='l',col='royalblue4',main='',xlab='nombre d\'itérations',ylab=expression(param[j]),lwd=2)
  lines(sim_gibbs$param[,2],type='l',col='springgreen3')
  legend("bottomleft",legend=c(expression(param[1]),expression(param[2])),lty =c(1,1),col=c("royalblue4","springgreen3"),lwd=c(2,2))
}






