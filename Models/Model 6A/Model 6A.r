# This is a Bayesian hierarchical model with multivariate normal likelihood, run in JAGS
# With hierarchical covariation in the random effects at the species and popultaion level. 
library(R2jags)
library(R2OpenBUGS)
library(loo)

dat1=read.csv("DataInputJul2020.csv")
plot(dat1$SCL,dat1$K)
plot(log(dat1$SCL),log(dat1$K))
summary(dat1)
table(dat1$Species,dat1$Population)
species=c("Kemp's ridley", "Green", "Hawksbill","Loggerhead","Olive ridley")
dat2=dat1
dat2$n[is.na(dat2$n)]=10
dat2$n[dat2$n>dat2$Turtles]=dat2$Turtles[dat2$n>dat2$Turtles]
dat2$Species=species[dat1$Species]
dat2$Species=factor(dat2$Species,levels=c("Green", "Loggerhead","Hawksbill","Olive ridley","Kemp's ridley"))
library(ggplot2)
jpeg("SampleFigure.jpg",units="in",width=6.5,height=6,res=1024)
ggplot(dat2,aes(x=log10(SCL),y=log10(K),col=Species,size=log10(n)))+geom_point()+
  labs(color = "Species",x=expression(log[10](italic(L)[infinity])),
       y=expression(log[10](italic(K))),size=expression(log[10](n)))+
  scale_color_manual(values=c("green","brown","red","blue","magenta"))+theme_classic()
dev.off()

write("model {
  for (i in 1:N) {
    Y[i, 1:2] ~ dmnorm(mu[i,], Sigma.inv[,])
    LL[i]<-logdensity.mnorm(Y[i,], mu[i,],Sigma.inv[,])
    for (j in 1:2)  {
      mu[i,j]<-alpha[j]+beta[x1[i],j]+population[pop[i],j]  
      resid[i,j]<-Y[i,j]-mu[i,j]
    }
  }
  #Prior for intercept for Linf and K
  for(j in 1:2) {
    alpha[j]~dnorm(0,1.0E-6)  #loop over Linf and K
  }
  # Fixed effect priors(e.g type of sample)
  for(j in 1:2) { #loop over Linf and K
    beta[1,j]<-0  #First level of factor is reference level, set equal to zero for both Linf (alpha1) and K (alpha2)
    for(i in 2:n.beta) { #Loop over number of levels factor has 
      beta[i,j]~dnorm(0,1.0E-6)  #uninformative prior
    }}
  # Random effect priors for species and population
  sp.cor~dunif(-1,1)
  pop.cor~dunif(-1,1)
  for(j in 1:2) {
    sp.tau[j]~dgamma(0.01,0.01)
    sp.sd[j]<-1/sqrt(sp.tau[j])
    pop.tau[j]~dgamma(0.01,0.01)  # assumes variance among pops same for all species
    pop.sd[j]<-1/sqrt(pop.tau[j])
    sp.cov[j,j]<-1/sp.tau[j]
    pop.cov[j,j]<-1/pop.tau[j]
    sp.mean[j]<-0
 }
  sp.cov[1,2]<-sp.cor*sp.sd[1]*sp.sd[2]
  sp.cov[2,1]<-sp.cor*sp.sd[1]*sp.sd[2]
  pop.cov[1,2]<-pop.cor*pop.sd[1]*pop.sd[2]
  pop.cov[2,1]<-pop.cor*pop.sd[1]*pop.sd[2]
  sp.prec[1 : 2 , 1 : 2]  <- inverse(sp.cov[ , ])
  pop.prec[1 : 2 , 1 : 2]  <- inverse(pop.cov[ , ])
  for(i in 1:n.sp) {
    species[i,1:2]~dmnorm(sp.mean[],sp.prec[,])
  }
  for(i in 1:n.pop) {
    population[i,1:2]~dmnorm(species[sp.pop[i],],pop.prec[,])
  } 
  Sigma.inv[1:2, 1:2] ~ dwish(R[,], 2)
  Sigma[1:2, 1:2]<- inverse(Sigma.inv[,])
  # predicted values for each species, for reference type of sample
  for(i in 1:n.sp) {
    for(j in 1:2) {
      mean.mu[i,j]<-alpha[j]+species[i,j]
      exp.mean[i,j]<-exp(mean.mu[i,j])
    }}   
  # predicted values for each species, for reference type of sample
  for(i in 1:n.pop) {
    for(j in 1:2) {
      mean.mupop[i,j]<-alpha[j]+population[i,j]
      exp.meanpop[i,j]<-exp(mean.mupop[i,j])
    }}  
  # population intercept without log
  for(j in 1:2) {intercept[j]<-exp(alpha[j])}
}

",file="MVNpop.txt")

#Set up data to pass to JAGS
x=dat1[!duplicated(dat1$Population),]
x=x[order(x$Population),]
sp.pop=x$Species  #This gives the species associated with each population

jagsdat6=list(Y=cbind(log(dat1$SCL),log(dat1$K)),x1=as.numeric(dat1$Method)
              ,N=dim(dat1)[1],R=matrix(c(2,0,0,2),2,2,byrow=TRUE),
              pop=as.numeric(dat1$Population),
              n.sp=length(unique(dat1$Species)),
              n.pop=length(unique(dat1$Population)),
              sp.pop=sp.pop,
              n.beta=length(unique(dat1$Method)))

# Set up initial values to pass to JAGS (one list of values for each of 2 chains)
init1=list(list(alpha=c(4,-1)),
  list(alpha=c(5,-0.5)))

# List parameters to save from MCMC
params6=c("alpha","beta","species","population","exp.mean",
          "sp.sd","pop.sd","sp.cor","pop.cor","exp.meanpop","intercept"
          ,"mu","resid","Sigma","LL")

res6=jags(jagsdat6,init1,params6,model.file="MVNpop.txt",
  n.chains=2,n.iter=210000,n.burnin=10000,n.thin=4)
write.csv(res6$BUGSoutput$summary,file="Model6.csv")
max(res6$BUGSoutput$summary[,"Rhat"])
min(res6$BUGSoutput$summary[,"n.eff"][res6$BUGSoutput$summary[,"n.eff"]>1])

#Make row in summary table for the run
run=6
restab=data.frame(Model=run)
restab$n.eff[1]=min(res6$BUGSoutput$summary[,"n.eff"][res6$BUGSoutput$summary[,"n.eff"]>1])
restab$Rhat[1]=round(max(res6$BUGSoutput$summary[,"Rhat"]),3)
restab$P.value[1]=NA
restab$DIC[1]=res6$BUGSoutput$DIC 
restab$pDIC[1]=res6$BUGSoutput$pD 
LLrows=paste0("LL[",1:jagsdat6$N,"]")
a=res6$BUGSoutput$sims.matrix[,LLrows]
b=waic(a)
restab$WAIC[1]=b$estimates["waic",1]
restab$pWAIC[1]=b$estimates["p_waic",2]
b=loo(a)
restab$LOOIC[1]=b$looic
restab$pLOOIC[1]=b$p_loo
restab$deviance[1]=res6$BUGSoutput$summary["deviance","mean"]
restab
write.csv(restab,paste0("waicres",run,".csv"))
pareto_k_table(b)
write.csv(pareto_k_table(b),file=paste0("paretoK",run,".csv"))


#Plot residuals
residrows=paste0("resid[",rep(1:jagsdat6$N,2),",",c(rep(1,jagsdat6$N),rep(2,jagsdat6$N)),"]")
murows=paste0("mu[",rep(1:jagsdat6$N,2),",",c(rep(1,jagsdat6$N),rep(2,jagsdat6$N)),"]")
df=data.frame(Expected=res6$BUGSoutput$summary[murows,"mean"],
              Residuals=res6$BUGSoutput$summary[residrows,"mean"],
              Variable=c(rep("Linf",jagsdat6$N),rep("K",jagsdat6$N)))
ggplot(df)+geom_point(aes(x=Expected,y=Residuals))+geom_abline(intercept=0,slope=0)+facet_wrap(~Variable,nrow=2, scales = "free")

# Look at variance parameters
sigrows=c("Sigma[1,1]","Sigma[1,2]","Sigma[2,1]","Sigma[2,2]")
a=res6$BUGSoutput$summary[sigrows,"50%"]
covmatError=matrix(a,2,2)
covmatSpecies=covmatPopulation=matrix(0,2,2)
covmatSpecies[1,1]=res6$BUGSoutput$summary["sp.sd[1]","50%"]^2
covmatSpecies[2,2]=res6$BUGSoutput$summary["sp.sd[2]","50%"]^2
covmatSpecies[1,2]=covmatSpecies[2,1]=res6$BUGSoutput$summary["sp.sd[2]","50%"]*res6$BUGSoutput$summary["sp.sd[1]","50%"]*res6$BUGSoutput$summary["sp.cor","50%"]
covmatPopulation[1,1]=res6$BUGSoutput$summary["pop.sd[1]","50%"]^2
covmatPopulation[2,2]=res6$BUGSoutput$summary["pop.sd[2]","50%"]^2
covmatPopulation[1,2]=covmatPopulation[2,1]=res6$BUGSoutput$summary["pop.sd[2]","50%"]*res6$BUGSoutput$summary["pop.sd[1]","50%"]*res6$BUGSoutput$summary["pop.cor","50%"]
round(covmatError,3)
round(covmatSpecies,3)
round(covmatPopulation,3)

