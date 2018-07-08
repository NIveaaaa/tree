library(rstan)
library(ggplot2)
library(gridExtra)
library(dplyr)
#load data
load("C:/Users/z5187692/Downloads/project/hmmtest.Rdata")

fit_ss=extract(hmm.test,permuted=F,inc_warmup=F)

ss_num=apply(fit_ss,3L,c)
ss_name=expand.grid(dimnames(fit_ss)[2])
ss_df=data.frame(ss_name,ss_num)

ss_df=cbind(ss_df,indexx=rep(1:5000,each=3))
ss_df$chains=as.factor(ss_df$chains)
p1=ggplot(ss_df,aes(indexx,beta0.3.,colour=chains))+geom_line(alpha=0.4)
p2=ggplot(ss_df,aes(indexx,beta1.3.,colour=chains))+geom_line(alpha=0.4)
p3=ggplot(ss_df,aes(indexx,eta.248.,colour=chains))+geom_line(alpha=0.4)
grid.arrange(p1,p2,p3,ncol=3)


chain1=filter(ss_df,chains=="chain:1")
chain2=filter(ss_df,chains=="chain:2")
chain3=filter(ss_df,chains=="chain:3")


plot(hmm.test,plotfun="trace",pars=c("mux[1]","mux[2]","mux[3]"),ncol=1)

plot(hmm.test,plotfun="trace",pars=c("beta2","beta1[3]","beta0[3]","eta[200]"),ncol=2)

plot(hmm.test,plotfun="trace",pars=c("rho","sigbeta1","sigbeta0","sigmaeta"))
