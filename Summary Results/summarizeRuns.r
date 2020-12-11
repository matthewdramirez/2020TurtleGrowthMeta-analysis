library(tidyverse)
library(ggplot2)
library(scales)
library(loo)
residrows=paste0("resid[",rep(1:jagsdat6$N,2),",",c(rep(1,jagsdat6$N),rep(2,jagsdat6$N)),"]")
murows=paste0("mu[",rep(1:jagsdat6$N,2),",",c(rep(1,jagsdat6$N),rep(2,jagsdat6$N)),"]")
LLrows=paste0("LL[",1:jagsdat6$N,"]")

resSummary<-list(res6$BUGSoutput,
                 res7$BUGSoutput,
                 res8$BUGSoutput,
                 res9$BUGSoutput,
                 res10$BUGSoutput,
                 res11$BUGSoutput,
                 res12$BUGSoutput,
                 res13$BUGSoutput,
                 res14$BUGSoutput)
residall<-data.frame(Expected=NULL,Residuals=NULL,Variable=NULL,Source=NULL)
modnames<-c("A","B","C","D","E","F","G","H","I")
for(i in 1:9)  {                
  residall<-bind_rows(residall, data.frame(Expected=resSummary[[i]]$summary[murows,"mean"],
              Residuals=resSummary[[i]]$summary[residrows,"mean"],
              Variable=c(rep("Linf",jagsdat10$N),rep("K",jagsdat10$N)),
#              VarLabel=c(rep("expression(log(italic(L)[infinity]))",jagsdat10$N),rep("expression(log(italic(K)))",jagsdat10$N)),
#            VarLabel=c(rep("log(italic(L)[infinity])",jagsdat10$N),rep("log(italic(K))",jagsdat10$N)),
          Model=as.character(rep(modnames[i],length(murows)))))
}              
head(residall)
write.csv(residall,"ResidAll.csv")
residall<-residall %>% mutate(Variable=factor(Variable))%>% mutate(Variable = recode_factor(Variable, "K"="log(italic(K))",
                                          "Linf"="log(italic(L)[infinity])"))
ggplot(residall)+geom_point(aes(x=Expected,y=Residuals))+geom_abline(intercept=0,slope=0)+theme_bw()+
  facet_grid(vars(Model),vars(Variable),scales = "free",
             label="label_parsed")
ggplot(residall,aes(sample=Residuals))+geom_qq()+geom_qq_line()+theme_bw()+
  facet_grid(vars(Model),vars(Variable),scales = "free",
             label="label_parsed")

             
ParetoTable<-NULL
for(i in 1:9) {
  ParetoTable<-cbind(ParetoTable,pareto_k_table(loo(resSummary[[i]]$sims.matrix[,LLrows]))[,"Proportion"])
}
dimnames(ParetoTable)[[2]]<-modnames
ParetoTable
write.csv(round(ParetoTable,2),file="ParetoTable.csv")

AllRes<-NULL
for(i in 1:9) {
  x<-read.csv(paste0("waicres",i+5,".csv"))
  AllRes<-bind_rows(AllRes,x)
}
AllRes$deltaDIC<-AllRes$DIC-min(AllRes$DIC)
AllRes$deltaWAIC<-AllRes$WAIC-min(AllRes$WAIC)
AllRes$deltaLOOIC<-AllRes$LOOIC-min(AllRes$LOOIC)
for(i in 1:9) AllRes$PercentK[i]<-ParetoTable[1,i]
AllRes$Model<-modnames
AllRes$Iterations<-c(res6$n.iter,res7$n.iter,res8$n.iter,res9$n.iter,res10$n.iter,res11$n.iter,res12$n.iter,res13$n.iter,res14$n.iter)-10000
AllRes$LOO.wts<-round(exp(-AllRes$deltaLOOIC/2)/sum(exp(-AllRes$deltaLOOIC/2)),2)
AllRes
write.csv(AllRes,file="ALlRes.csv")

######
#Plot vonbert curvef

age=seq(0,70)
vbvals<-expand.grid(sim=1:50000,sp=1:5,Age=age)
head(vbvals)
for(i in 1:5) {
  vbvals$K[vbvals$sp==i]=rep(res9$BUGSoutput$sims.matrix[,paste0("exp.mean[",i,",2]")],length(age))
  vbvals$Linf[vbvals$sp==i]=rep(res9$BUGSoutput$sims.matrix[,paste0("exp.mean[",i,",1]")],length(age))
}
vbvals<-vbvals %>% mutate(Length=Linf*(1-exp(-K*Age))) 
head(vbvals)

vbplot<-vbvals %>% group_by(Age,sp) %>% 
  summarize(Mean=mean(Length),
            min=quantile(Length,0.025),
            max=quantile(Length,0.975))
head(vbplot)
dim(vbplot)
x=levels(latres$Species)[c(5,1,3,2,4)]
x
vbplot$Species=x[vbplot$sp]
vbplot$Species=factor(vbplot$Species,levels=c("Green", "Loggerhead","Hawksbill","Olive ridley","Kemp's ridley"))
#for(i in 1:5) {
#  vbplot[vbplot$Mean>max(dat1$SCL[dat1$Species==i])&vbplot$sp==i,c("Mean","min","max")]<-NA
#}
summary(vbplot)
cbPalette3 <- c("#228833", "#EE6677",  "#4477AA","#CCBB44", "#AA3377")


ggplot(vbplot,aes(x=Age,y=Mean,col=Species,fill=Species))+
  geom_line(size=1.5)+geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2)+theme_classic()+
  scale_fill_manual(values=cbPalette3)+
  scale_color_manual(values=cbPalette3)+
  ylab("Length (cm)")
