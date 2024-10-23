J<-5000
theta<-rbinom(J,1,0.05)
epsilon<-rnorm(J,0,1)
Y<-theta+epsilon
test<-rep(0,J)


for(i in 1:J){
  if(Y[i]>1.96){test[i]=1}}

#test-theta
S=0
for(j in 1:J){
  
  if(theta[j]<test[j]){S=S+1}
}
S
S/J
sum(theta)
sum(test)
t.test(c(Y[1],mu=0))

###
n<-10
J<-5000
T<-rep(0,J)
resultat=rep(0,J)
theta<-rep(0,J)
s<-0
alpha<-0.005
for(i in 1:J){
  a=rnorm(n,0,1)
  b=rbinom(1,1,0.05)
  theta[i]<-b
  X<-a+b
  T[i]<-sqrt(n)*(mean(X))/sqrt(var(X))
  #pval
  #p<-2*(1-pt(abs(T),n-1))
  
  resultat[i]<-ifelse(T[i]<=qt(1-alpha,n-1),0,1)
  s<-ifelse(resultat[i]>b,s+1,s)
}
s # les faux positif
sum(theta) #le nombre de vrai site
sum(resultat) #le nombre de fois ou on conclut H1 a tord ou a raison
theta