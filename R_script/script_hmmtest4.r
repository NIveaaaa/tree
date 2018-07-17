
library("rstan")
#setwd("~/project/data")
logy=read.table("logy.txt",head=FALSE)
year=read.table("year.txt",head=FALSE)
age=read.table("age.txt",head=FALSE)
temp=read.table("temp.txt",head=FALSE)
sd.temp=sd(as.numeric(temp),na.rm=TRUE)
mean.temp=mean(as.numeric(temp),na.rm=TRUE)
stemp=scale(as.numeric(temp))/2 #standard deviation is 0.5

mlogy=logy
mlogy[is.na(logy)]=0

sage=apply(age,2,function(x) (x-100)/100)

msage=sage
msage[is.na(msage)]=0

#start and end index of tree ring observation
s=apply(sage,1,function(x) min(which(is.na(x)=="FALSE")))
e=apply(sage,1,function(x) max(which(is.na(x)=="FALSE")))

#sort out all trees with support from recorded temperature.
tree_withsup=c(1:247)[s<418 & e>=418]
tree_withsup=c(tree_withsup,c(1:247)[s>418])

mlogy_withsup=mlogy[tree_withsup,]
msage_withsup=msage[tree_withsup,]

s_withsup=apply(sage[tree_withsup,],1,function(x) min(which(is.na(x)=="FALSE")))
e_withsup=apply(sage[tree_withsup,],1,function(x) max(which(is.na(x)=="FALSE")))
# getrowIndex=function(x){
#   index=which(x!=0)
#   index
# }

#nt_1=which(e>=418)

# indexpool=apply(age[,1:500],2,function(x) {which(x!=0)})
# 
# len=as.numeric(do.call(rbind, lapply(indexpool, length)))
# 
# b=matrix(0,nrow=247,500)
# for(i in 1:500){
#   b[1:len[i],i]=unlist(indexpool[i])
# }



standata=list(K=3,span=500,nt=121,start=s_withsup,end=e_withsup,
              age=as.matrix(msage_withsup),N_mis=417,y=as.matrix(mlogy_withsup),
              x_obs=as.vector(stemp)[418:500],alpha=c(1,1,1))



inits_stan=function(){
  mu_x=rnorm(3,0,0.5)
  sigmay=runif(1,0.3,0.4)
  sigmax=runif(1,0.4,0.6)
  sigbeta0=runif(1,0.4,0.6)
  sigbeta1=runif(1,0.5,0.7)
  siggamma=runif(1,0.6,0.8)
  sigmaeta=runif(1,0.3,0.4)
  betamu0=runif(1,6,6.5)
  beta0=rnorm(121,6.5,0.5)
  betamu1=runif(1,-0.5,-0.3)
  beta1=rnorm(121,-0.4,0.1)
  gammamu=runif(1,0.4,0.5)
  PI=c(1/3,1/3,1/3)
  list(mu_x=mu_x,ky=sigmay,kx=sigmax,sigbeta0=sigbeta0,sigbeta1=sigbeta1,
       siggamma1=siggamma,keta=sigmaeta,beta0mu=betamu0,beta1mu=betamu1,
       beta1=beta1,beta0=beta0,gammamu=gammamu,PI=PI
       )
}
inits=list()
for(i in 1:3)inits[[i]] = inits_stan()
#stanfit0<-stan(file="hmm.stan",init=inits,data=standata,chains=1,iter=200)
#control=list(max_treedepth=15)
#setwd("~/project/Basic model")
hmm.test4<-stan(file="hmm_test4.stan",data=standata,chains=3,init=inits,iter=10000)

# 
# save(stanfit5,file="stanfit5.Rdata")
# 
# stan2coda <- function(fit) {
#   mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
# }
# 
# mcmcfit5=stan2coda(stanfit5)


# plot(stanfit3,par=c("lp__","mux[1]"))
# 
# params1 <- as.data.frame(extract(stanfit3, permuted=FALSE)[,1,])
# params2<-as.data.frame(extract(stanfit3, permuted=FALSE)[,2,])
# params3<-as.data.frame(extract(stanfit3, permuted=FALSE)[,3,])
# par(mar = c(4, 4, 0.5, 0.5))
# plot(params1$"mux[2]", params1$"mux[3]", col=2, pch=16, cex=0.8,
#      xlab="mu1", xlim=c(-3, 3), ylab="mu2", ylim=c(-3, 3))
# lines(0.08*(1:100) - 4, 0.08*(1:100) - 4, col="grey", lw=2)
# params1
#   