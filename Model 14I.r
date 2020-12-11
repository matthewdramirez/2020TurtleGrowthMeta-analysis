# This is a Bayesian hierarchical model with multivariate normal likelihood, run in JAGS
# With hierarchical covariation in the random effects at the species and popultaion level. 
library(R2jags)
library(R2OpenBUGS)
library(loo)
library(ggplot2)
library(tidyverse)

dat1=read.csv("DataInputJul2020.csv")
plot(dat1$SCL,dat1$K)
plot(log(dat1$SCL),log(dat1$K))
dat1$minRange<-dat1$SCL_min/dat1$SCL
dat1$maxRange<-dat1$SCL_max/dat1$SCL
summary(dat1)
ggplot(dat1,aes(x=minRange,y=maxRange,color=factor(Species),fill=factor(Species)))+geom_point()+stat_smooth(method="lm")
ggplot(dat1,aes(x=maxRange,y=SCL,color=factor(Species),fill=factor(Species)))+stat_smooth(method="lm")+theme_bw()
ggplot(dat1,aes(x=minRange,y=SCL,color=factor(Species),fill=factor(Species)))+stat_smooth(method="lm")+theme_bw()
ggplot(dat1,aes(x=abs(Lat),y=SCL,color=factor(Species),fill=factor(Species)))+stat_smooth(method="lm")+theme_bw()
ggplot(dat1,aes(x=abs(Lat),y=K,color=factor(Species),fill=factor(Species)))+stat_smooth(method="lm")+theme_bw()
table(dat1$minRange<1/4)
table(dat1$maxRange>3/4)
table(dat1$minRange<1/4,dat1$maxRange>3/4)
dat1$SizeRange<-rep(0,dim(dat1)[1])
dat1$SizeRange[dat1$minRange<1/4 & dat1$maxRange>3/4]=1  
dat1$SizeRange[dat1$minRange<1/4 & dat1$maxRange<3/4]=2
dat1$SizeRange[dat1$minRange>1/4 & dat1$maxRange>3/4]=3
dat1$SizeRange[dat1$minRange>1/4 & dat1$maxRange<3/4]=4
table(dat1$SizeRange)
# 1 = large and small, 2=small no large, 3=large no small, 4= no large or small
table(dat1$Species,dat1$Population)
species=c("Kemp's ridley", "Green", "Hawksbill","Loggerhead","Olive ridley")
dat2=dat1
dat2$n[is.na(dat2$n)]=10
dat2$n[dat2$n>dat2$Turtles]=dat2$Turtles[dat2$n>dat2$Turtles]
dat2$Species=species[dat1$Species]
dat2$Species=factor(dat2$Species,levels=c("Green", "Loggerhead","Hawksbill","Olive ridley","Kemp's ridley"))
#jpeg("SampleFigure.jpg",units="in",width=6.5,height=6,res=1024)
ggplot(dat2,aes(x=log(SCL),y=log(K),col=Species,size=log10(n)))+geom_point()+
  labs(color = "Species",x=expression(log[10](italic(L)[infinity])),
       y=expression(log[10](italic(K))),size=expression(log[10](n)))+
  scale_color_manual(values=c("green","brown","red","blue","magenta"))+theme_classic()
#dev.off()

write("model {
  for (i in 1:N) {
    Y[i, 1:2] ~ dmnorm(mu[i,], Sigma.inv[,])
    LL[i]<-logdensity.mnorm(Y[i,], mu[i,],Sigma.inv[,])
    for (j in 1:2)  {
      mu[i,j]<-alpha[j]+beta3[j]*x3[i]+species[sp[i],j]
      resid[i,j]<-Y[i,j]-mu[i,j]
    }
  }
  #Prior for intercept for Linf and K
  for(j in 1:2) {
    alpha[j]~dnorm(0,1.0E-6)  #loop over Linf and K
    beta3[j]~dnorm(0,1.0E-6)
  }
  # Random effect priors for species
  sp.cor~dunif(-1,1)
  for(j in 1:2) {
    sp.tau[j]~dgamma(0.01,0.01)
    sp.sd[j]<-1/sqrt(sp.tau[j])
    sp.cov[j,j]<-1/sp.tau[j]
    sp.mean[j]<-0
 }
  sp.cov[1,2]<-sp.cor*sp.sd[1]*sp.sd[2]
  sp.cov[2,1]<-sp.cor*sp.sd[1]*sp.sd[2]
  sp.prec[1 : 2 , 1 : 2]  <- inverse(sp.cov[ , ])
  for(i in 1:n.sp) {
    species[i,1:2]~dmnorm(sp.mean[],sp.prec[,])
  }
  Sigma.inv[1:2, 1:2] ~ dwish(R[,], 2)
  Sigma[1:2, 1:2]<- inverse(Sigma.inv[,])
  # predicted values for each species, for reference type of sample
  for(i in 1:n.sp) {
    for(j in 1:2) {
      mean.mu[i,j]<-alpha[j]+species[i,j]
      exp.mean[i,j]<-exp(mean.mu[i,j])
      for(lat in 1:42) {
       exp.lat.mean[i,j,lat]<-exp(alpha[j]+species[i,j]+beta3[j]*lat)
    }
    }}   
  # population intercept without log
  for(j in 1:2) {intercept[j]<-exp(alpha[j])}
}",file="MVNspLatNoMethod.txt")

#Set up data to pass to JAGS
jagsdat14=list(Y=cbind(log(dat1$SCL),log(dat1$K)),
              x3=abs(dat1$Lat),
              N=dim(dat1)[1],R=matrix(c(2,0,0,2),2,2,byrow=TRUE),
              sp=as.numeric(dat1$Species),
              n.sp=length(unique(dat1$Species)))

# Set up initial values to pass to JAGS (one list of values for each of 2 chains)
init1=list(list(alpha=c(4,-1)),
           list(alpha=c(5,-0.5)))

# List parameters to save from MCMC
params14=c("alpha","beta3","species","exp.mean","exp.lat.mean",
          "sp.sd","sp.cor","intercept"
          ,"mu","resid","Sigma","LL")

res14=jags(jagsdat14,init1,params14,model.file="MVNspLatNoMethod.txt",
          n.chains=2,n.iter=110000,n.burnin=10000,n.thin=4)
write.csv(res14$BUGSoutput$summary,file="Model14.csv")
max(res14$BUGSoutput$summary[,"Rhat"])
min(res14$BUGSoutput$summary[,"n.eff"][res14$BUGSoutput$summary[,"n.eff"]>1])

#Make row in summary table for the run
run=14
restab=data.frame(Model=run)
restab$n.eff[1]=min(res14$BUGSoutput$summary[,"n.eff"][res14$BUGSoutput$summary[,"n.eff"]>1])
restab$Rhat[1]=round(max(res14$BUGSoutput$summary[,"Rhat"]),3)
restab$P.value[1]=NA
restab$DIC[1]=res14$BUGSoutput$DIC 
restab$pDIC[1]=res14$BUGSoutput$pD 
LLrows=paste0("LL[",1:jagsdat14$N,"]")
a=res14$BUGSoutput$sims.matrix[,LLrows]
b=waic(a)
restab$WAIC[1]=b$estimates["waic",1]
restab$pWAIC[1]=b$estimates["p_waic",2]
b=loo(a)
restab$LOOIC[1]=b$looic
restab$pLOOIC[1]=b$p_loo
restab$deviance[1]=res14$BUGSoutput$summary["deviance","mean"]
restab
write.csv(restab,paste0("waicres",run,".csv"))
pareto_k_table(b)
write.csv(pareto_k_table(b),file=paste0("paretoK",run,".csv"))


#Plot residuals
residrows=paste0("resid[",rep(1:jagsdat7$N,2),",",c(rep(1,jagsdat7$N),rep(2,jagsdat7$N)),"]")
murows=paste0("mu[",rep(1:jagsdat7$N,2),",",c(rep(1,jagsdat7$N),rep(2,jagsdat7$N)),"]")
df=data.frame(Expected=res7$BUGSoutput$summary[murows,"mean"],
              Residuals=res7$BUGSoutput$summary[residrows,"mean"],
              Variable=c(rep("Linf",jagsdat7$N),rep("K",jagsdat7$N)))
ggplot(df)+geom_point(aes(x=Expected,y=Residuals))+geom_abline(intercept=0,slope=0)+facet_wrap(~Variable,nrow=2, scales = "free")

# Look at variance parameters
sigrows=c("Sigma[1,1]","Sigma[1,2]","Sigma[2,1]","Sigma[2,2]")
a=res7$BUGSoutput$summary[sigrows,"50%"]
covmatError=matrix(a,2,2)
covmatSpecies=covmatPopulation=matrix(0,2,2)
covmatSpecies[1,1]=res7$BUGSoutput$summary["sp.sd[1]","50%"]^2
covmatSpecies[2,2]=res7$BUGSoutput$summary["sp.sd[2]","50%"]^2
covmatSpecies[1,2]=covmatSpecies[2,1]=res7$BUGSoutput$summary["sp.sd[2]","50%"]*res7$BUGSoutput$summary["sp.sd[1]","50%"]*res7$BUGSoutput$summary["sp.cor","50%"]
covmatError
covmatSpecies

##Calculate latitude plots
latres<-expand.grid(Species=species,Variable=c("Linf","K"), Latitude=1:42)
latrows<-paste0("exp.lat.mean[",as.numeric(latres$Species),",",as.numeric(latres$Variable),",",latres$Latitude,"]")
latres<-cbind(latres,res14$BUGSoutput$summary[latrows,c("mean","2.5%","97.5%")])
latres<-latres %>% mutate(Variable2 = recode_factor(Variable, "K"="italic(K)",
                                                   "Linf"="italic(L)[infinity]"))
latdat<-rbind(dat1[,c("Species","Lat")],dat1[,c("Species","Lat")])
latdat$Data<-c(dat1$SCL,dat1$K)
latdat$Variable=rep(c("Linf","K"),each=dim(dat1)[1])
latdat$Species<-species[latdat$Species]
latdat<-latdat %>% mutate(Variable2 = recode_factor(Variable, "K"="italic(K)","Linf"="italic(L)[infinity]"),
                          Lat=abs(Lat))

latrange<-data.frame(spname=levels(dat1$Spname),min=0,max=0)
for(i in 1:5) latrange[i,c("min","max")]<-range(abs(dat1$Lat[as.character(dat1$Spname)==as.character(latrange$spname[i])])) 
latrange$Species=sort(levels(latres$Species))
latrange
for(i in 1:5) {
  latres[(latres$Latitude<(latrange$min[i]-5) | (latres$Latitude>latrange$max[i]+5)) &latres$Species==latrange$Species[i],c("2.5%","mean","97.5%")]<-NA
}

latres$Species=factor(latres$Species,levels=c("Green", "Loggerhead","Hawksbill","Olive ridley","Kemp's ridley"))
cbPalette3 <- c("#228833", "#EE6677",  "#4477AA","#CCBB44", "#AA3377")

#Both K and Linf
ggplot(latres,aes(x=Latitude,y=mean,col=Species,fill=Species))+
  geom_line(size=1.5)+geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3)+theme_bw()+
  facet_grid(rows=vars(Variable2),scales = "free",
             label="label_parsed")+ylab("Mean parameter value")+
  geom_point(data=latdat,aes(x=Lat,y=Data,col=Species),size=2)+
  scale_fill_manual(values=cbPalette3)+
  scale_color_manual(values=cbPalette3)

#Just K

g1<-ggplot(filter(latres,Variable=="K"),aes(x=Latitude,y=mean,col=Species,fill=Species))+
  geom_line(size=1.5)+geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3)+theme_classic()+
  ylab(expression(italic(K)))+
  geom_point(data=filter(latdat,Variable=="K"),aes(x=Lat,y=Data,col=Species),size=2)+
  scale_fill_manual(values=cbPalette3)+
  scale_color_manual(values=cbPalette3)
g1
ggsave(filename = "Figure4.eps",
       width=6.5,height=4.5,units="in",
       plot = print(g1),
       device = cairo_ps)

#Just Linf
ggplot(filter(latres,Variable=="Linf"),aes(x=Latitude,y=mean,col=Species,fill=Species))+
  geom_line(size=1.5)+geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3)+theme_classic()+
  ylab(expression(italic(L)[infinity] ))+
  geom_point(data=filter(latdat,Variable=="Linf"),aes(x=Lat,y=Data,col=Species),size=2)+
  scale_fill_manual(values=cbPalette3)+
  scale_color_manual(values=cbPalette3)

save(file="temp",list=c("latres","latdat","cbPalette3"))
