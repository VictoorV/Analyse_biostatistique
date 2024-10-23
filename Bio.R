data(golub)
row.names(golub)=golub.gnames[,3]
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))
golub[1,gol.fac=="ALL"]
boxplot(golub[720,] ~ gol.fac, method="jitter")
mean(golub[720,1:27])
mean(golub[720,28:38])
cat<-rep("A",27)
cat<-append(cat,rep("B",11))
plot(golub[720,],col=ifelse(cat == "A","blue","red"),pch=15)
legend("topright", legend = c("ALL","AML"), col = c("blue","red"), pch = c(15, 15))
##
plot(golub[,1], golub[,38], xlab = 'Patient 1 (ALL)', ylab = 'Patient 38 (AML)',col=c("blue","red"),pch=c(15,15))
legend("topleft", legend = c("Patient 1","Patient 38"), col = c("blue","red"), pch = c(15, 15))
##
hist(golub[,1:27],col="blue",xlab="ALL",ylab="Frequence",breaks=10)
hist(golub[,28:38],col="blue",xlab="AML",ylab="Frequence",breaks=10)
##
X<-rep(0,3051)
pval<-rep(0,3051)
for(i in 1:3051){
  test<-t.test(golub[i,1:27],golub[i,28:38], alternative = "two.sided",var.equal=FALSE,conf.level=0.95 )
  X[i]<-test$statistic
  pval[i]<-test$p.value
}
hist(X,breaks=100,xlim=c(-10,10),col="blue",xlab="T",ylab="Frequence")
hist(pval,breaks=100,col="blue",ylab="Frequence",xlab="pval")
##
pval_s<-sort(pval)
pval_bonf<-rep(0,3051)
res_bonf<-0
for(i in 1:3051){
  pval_bonf[i]<-pval_s[i]*3051
  if(pval_s[i]<=0.05/3051){
    res_bonf<-res_bonf+1
  }
}
pval_sid<-rep(0,3051)
res_sid<-0
for(i in 1:3051){
  pval_sid[i]<-1-(1-pval_s[i])^3051
  if(pval_s[i]<=1-(1-0.05)^(1/3051)){
    res_sid<-res_sid+1
  }
}
pval_HB<-rep(0,3051)
res_HB<-0
for(i in 1:3051){
  pval_HB[i]<-pval_s[i]*(3051+1-i)
  if(pval_s[i]<=0.05/(3051+1-i)){
    res_HB<-res_HB+1
  }
}
pval_BH<-rep(0,3051)
res_BH<-0
for(i in 1:3051){
  pval_BH[i]<-pval_s[i]*3051/i
  if(pval_s[i]<=0.05*i/3051){
    res_BH<-res_BH+1
  }
}
pval_BY<-rep(0,3051)
res_BY<-0
for(i in 1:3051){
  L<-3051*sum(c(1/(1:3051)))
  pval_BY[i]<-pval_s[i]*L/i
  if(pval_s[i]<=i*0.05/L){
    res_BY<-res_BY+1
  }
}
##
p.adjust(pval_s,method="bonferroni")
sum((p.adjust(pval_s,method="BH")<0.05))
sum((p.adjust(pval_s,method="BY")<=0.05))
sum((p.adjust(pval_s,method="bonferroni")<0.05))
##
plot(pval_s,ifelse(pval_bonf>=1,1,pval_bonf),col="#33FF00",type='l',xlab="p-valeurs triées",ylab="p-valeurs corrigées")
points(pval_s,ifelse(pval_sid>=1,1,pval_sid),col="#CCCC33",type='l')
points(pval_s,ifelse(pval_HB>=1,1,pval_HB),col="#333399",type='l')
points(pval_s,ifelse(pval_BH>=1,1,pval_BH),col="#990000",type='l')
points(pval_s,ifelse(pval_BY>=1,1,pval_BY),col="#FF9933",type='l')
legend("bottomright", legend = c("BONF","SID","HB","BH","BY"), col = c("#33FF00","#CCCC33","#333399","#990000","#FF9933"), lty = c(1,1,1,1,1))

res<-data.frame(matrix(nrow=1,ncol=5))
rownames(res)<-c("Resultats")
colnames(res)<-c("BONF","SID","HB","BH","BY")
res[1,1]<-res_bonf
res[1,2]<-res_sid
res[1,3]<-res_HB
res[1,4]<-res_BH
res[1,5]<-res_BY