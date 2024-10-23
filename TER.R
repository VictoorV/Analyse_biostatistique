# Fonction p-valeurs
liste_pval<-function(n,N,mu=0,mu0=0){
  a<-c()
  for(i in 1:N){
    X<-rnorm(n,mu,1)
    T<-sqrt(n)*(mean(X)-mu0)/sqrt(var(X))
    p<-2*(1-pt(abs(T),n-1))
    a=append(a,2*(1-pt(abs(T),n-1)))
  }
  return(a)
}
# Affichage des resultats
hist(liste_pval(1000,1000),freq=FALSE,main ="Histogramme des p-valeurs",xlab="p-valeurs",ylab="Densité")
curve(dunif(x), add=TRUE, col="red") # La loi des p-val est une loi uniforme sur [0,1]





# Fonction test de Student/Intervalle de confiance
Student<-function(n,mu,mu0,B){
  # Intervalle de confiance a 95% pour mu
  X<-rnorm(n,mu,1)
  mean<-mean(X)
  T<-sqrt(n)*(mean(X)-mu0)/sqrt(var(X))
  pval<-2*(1-pt(abs(T),n-1))
  print("Moyenne, T, pval")
  print(c(mean,T,pval))
  inter1<-c(mean-qt(0.975,n-1)*sqrt(var(X))/sqrt(n),mean+qt(0.975,n-1)*sqrt(var(X))/sqrt(n))
  print("Intervalle de confiance a 95% pour la moyenne")
  print(inter1)
  # Intervalle de confiance a 95% bootstrap pour la p-valeur
  a<-c()
  VecBoot<-replicate(B,sample(X,n,replace=TRUE))
  for(i in 1:B){
    Tb<-sqrt(n)*(mean(VecBoot[,i])-mu0)/sqrt(var(VecBoot[,i]))
    pb<-2*(1-pt(abs(Tb),n-1))
    a=append(a,pb)
  }
  a=sort(a) # Rangement de a par ordre croissant
  inter2<-c(a[ceiling(B*0.05/2)],a[ceiling(B*(1-0.05/2))]) # Intervalle de confiance bootstrap a 95% pour la pval
  print("Intervalle de confiance a 95% pour la pval")
  print(inter2)
  res=c(0,0) # Pour verifier si nos parametres sont dans les intervalles de confiance
  if(inter1[1] < mu & mu < inter1[2]){
    res[1]=1
  }
  if(inter2[1] < pval & pval < inter2[2]){
    res[2]=1
  }
  return(res)
}
#
res_mean=0
res_pval=0
for(i in 1:1000){
  a=Student(100,0,0,3000)
  res_pval=res_pval+a[2]
}
#
Student(10,0,0,200)
#
n<-100
mu=5
mu0=8
Y<-rexp(n,1/mu) # Moyenne de l'echantillon : mu 
mean<-mean(Y)
t.test(Y,mu=mu0,conf.level=0.95)
T<-sqrt(n)*(mean(Y)-mu0)/sqrt(var(Y))
pval<-2*(1-pt(abs(T),n-1))
print(c(mean,T,pval))
sort(Y)
X<-rnorm(n,0,1)
sort(X)
sum(Y[Y<6])/9

### puissance
par(mfrow=c(1,3))
n<-30
f<-function(mu){
  X<-rnorm(n,mu,1)
  T<-sqrt(n)*(mean(X)-mu)/sqrt(var(X))
  #ifelse(mu>=0,1-pt(qt(0.95,n-1)-sqrt(n)*(mu)/sqrt(var(X)),n-1),(1-pt(qt(0.95,n-1)+sqrt(n)*(mu)/sqrt(var(X)),n-1)))
  return(1-pt(qt(0.975,n-1)-sqrt(n)*(mu)/sqrt(var(X)),n-1)+pt(-qt(0.975,n-1)-sqrt(n)*(mu)/sqrt(var(X)),n-1))
}
plot(seq(-2,2,0.01),f(seq(-2,2,0.01)),type='l',xlab='mu',ylab='P_H1(|T|>k_a)',ylim=c(0,1))
points(0,0.05,col='red',pch=19)
legend("bottomright", legend = c("n=30"))

curve(dnorm(x,5,1),from=0,to=10,col='red')
curve(dexp(x,1/5),add=TRUE,col='blue')
legend("topright", legend = c("Densite de N(5,1)","Densite de E(1/5)"),col=c("red","blue"),lty=c(1,1))
