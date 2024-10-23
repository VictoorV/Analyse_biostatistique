
f<-function(x){
  return(1-(2*pnorm(0.05*sqrt(log(1000)*2)))**1000)
}
(2*pnorm(0.05*sqrt(log(5)*2))-1)**5
2*150*(pnorm(0.05*(150))-1)
curve(pnorm(x),from=-10,to=10,add=TRUE,col="red")
curve(pnorm(log(x)),from=-10,to=10,add=TRUE)



a<-0.2
J<-500
P<-integrate(dnorm,lower=-a*sqrt(2*log(J)),upper=a*sqrt(2*log(J)))
1-P[["value"]]**J
g<-function(x){
  return(1-(2*integrate(dnorm,lower=0,upper=a*sqrt(2*log(x)))[["value"]])**x)
}
h<-function(x){
  return(1-exp(0.5*x*log(1-x**(-4*a*a/pi))))
}
i<-function(x){
  return(g(x)-h(x))
}

curve(Vectorize(g)(x),from=1,to=5,col="red")
curve(Vectorize(h)(x),from=1,to=5,add=TRUE)
curve(Vectorize(i)(x),from=1,to=50)

i(11)
