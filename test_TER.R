# Tests multiples
methode1<-function(n,alpha,mu,sd,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pval<-rep(0,J)
  decision<-rep(0,J)
  theta<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sd) #vecteur
    theta[i]<-rbinom(1,1,p)
    Y<-theta[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pval[i]<-1-pt(T[i],n-1)
    if(pval[i] <= alpha){
      decision[i]<-1
    }
    else{
      decision[i]<-0
    }
    if(decision[i]==1 & theta[i]==0){
      faux_pos=faux_pos+1
      #print(c(sum(Y)/n,pval[i],T[i]))
    }
    if(decision[i]==0 & theta[i]==1){
      faux_neg=faux_neg+1
      #print(c(sum(Y)/n,pval[i],T[i]))
    }
  }
  print("Nb vrais sites, Resultats, Nb faux pos, Nb faux neg")
  return(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision)-faux_pos+faux_neg))
}

J<-5000
alpha<-0.05
mu<-0
sd<-1
p<-0.05
methode1(10,alpha,mu,sd,J,p)
J*alpha*(1-p)

methode1(10,1-(1-alpha)^(1/J),mu,sd,J,p)

fp<-c()
fn<-c()
for(i in 1:100){
  res<-methode1(10,0.05,mu,sd,J,p)
  res
  fp<-append(fp,res[3])
  fn<-append(fn,res[4])
}
plot(c(1:100),fn,col="red",type='l')
points(c(1:100),fp,col="blue",type='l')

methode1(10,1-(1-alpha)^(1/10),mu,sd,10,p)

### HB
HB<-function(n,alpha,mu,sd,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pval<-rep(0,J)
  decision<-rep(0,J)
  theta<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sd) #vecteur
    theta[i]<-rbinom(1,1,p)
    Y<-theta[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pval[i]<-1-pt(T[i],n-1)
  }
  pval_s<-sort(pval)
  for(j in (J):1){
    if(pval_s[j]<alpha/(J+1-j)){
      decision[j]=1
      if(theta[order(pval)[j]]==0){
        faux_pos<-faux_pos+1
      }
    }
    else{
      decision[j]=0
      if(theta[order(pval)[j]]==1){
        faux_neg<-faux_neg+1
      }
    }
  }
  print("Nb vrais sites, Resultats, Nb faux pos, Nb faux neg")
  print(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision) - faux_pos + faux_neg))
}

HB(10,0.05,0,1,5000,0.05)
















methode3 <- function(n, alpha, mu, sd, J, p) {
  faux_pos <- 0
  faux_neg <- 0
  T <- rep(0, J)
  pval <- rep(0, J)
  decision <- rep(0, J)
  theta <- rep(0, J)
  
  for (i in 1:J) {
    bruit <- rnorm(n, mu, sd)  # Vecteur
    theta[i] <- rbinom(1, 1, p)
    Y <- theta[i] + bruit
    T[i] <- sqrt(n) * (mean(Y) - 0) / sqrt(var(Y))
    pval[i] <- 1 - pt(T[i], n - 1)
  }
  
  # Correction de Benjamini-Hochberg avec un seuil plus strict
  sorted_indices <- order(pval)
  m <- length(pval)
  thresholds <- 4*(1:m) * alpha / (m)  # Seuil plus strict
  
  for (i in 1:m) {
    index <- sorted_indices[i]
    if (pval[index] <= thresholds[i]) {
      decision[index] <- 1
    } else {
      break  # Arrêter le processus de rejet si le seuil n'est pas atteint
    }
  }
  
  for (i in 1:J) {
    if (decision[i] == 1 && theta[i] == 0) {
      faux_pos <- faux_pos + 1
    }
    
    if (decision[i] == 0 && theta[i] == 1) {
      faux_neg <- faux_neg + 1
    }
  }
  
  cat("Nb vrais sites, Résultats, Nb faux pos, Nb faux neg\n")
  cat(sum(theta), sum(decision), faux_pos, faux_neg, sum(decision) - faux_pos + faux_neg, "\n")
  
  return(pval)
}

pval<-methode2(10,0.05,0,1,5000,0.05)
westfall_young <- function(n, alpha, mu, sd, J, p) {
  faux_pos <- 0
  faux_neg <- 0
  T <- rep(0, J)
  pval <- rep(0, J)
  decision <- rep(0, J)
  theta <- rep(0, J)
  
  for (i in 1:J) {
    bruit <- rnorm(n, mu, sd)  # Vecteur
    theta[i] <- rbinom(1, 1, p)
    Y <- theta[i] + bruit
    T[i] <- sqrt(n) * (mean(Y) - 0) / sqrt(var(Y))
    pval[i] <- 1 - pt(T[i], n - 1)
  }
  
  # Correction de Westfall-Young
  sorted_indices <- order(pval)
  m <- length(pval)
  thresholds <- rep(0, m)
  thresholds[m] <- alpha
  
  for (i in (m - 1):1) {
    thresholds[i] <- min((i/m) * alpha, thresholds[i + 1])
  }
  
  for (i in 1:m) {
    index <- sorted_indices[i]
    if (pval[index] <= thresholds[i]) {
      decision[index] <- 1
    } else {
      break  # Arrêter le processus de rejet si le seuil n'est pas atteint
    }
  }
  
  for (i in 1:J) {
    if (decision[i] == 1 && theta[i] == 0) {
      faux_pos <- faux_pos + 1
    }
    
    if (decision[i] == 0 && theta[i] == 1) {
      faux_neg <- faux_neg + 1
    }
  }
  
  cat("Nb vrais sites, Résultats, Nb faux pos, Nb faux neg\n")
  cat(sum(theta), sum(decision), faux_pos, faux_neg, sum(decision) - faux_pos + faux_neg, "\n")
}
westfall_young(10,0.05,0,1,J,0.05)

methode2 <- function(n, alpha, mu, sd, J, p) {
  faux_pos <- 0
  faux_neg <- 0
  T <- rep(0, J)
  pval <- rep(0, J)
  decision <- rep(0, J)
  theta <- rep(0, J)
  
  for (i in 1:J) {
    bruit <- rnorm(n, mu, sd)  # Vecteur
    theta[i] <- rbinom(1, 1, p)
    Y <- theta[i] + bruit
    T[i] <- sqrt(n) * (mean(Y) - 0) / sqrt(var(Y))
    pval[i] <- 1 - pt(T[i], n - 1)
  }
  
  # Correction de Benjamini-Hochberg pour réduire les faux négatifs
  sorted_indices <- order(pval)
  m <- length(pval)
  thresholds <- (1:m) * alpha / m
  
  for (i in 1:m) {
    index <- sorted_indices[i]
    if (pval[index] <= thresholds[i]) {
      decision[index] <- 1
    } else {
      break  # Arrêter le processus de rejet si le seuil n'est pas atteint
    }
  }
  
  for (i in 1:J) {
    if (decision[i] == 1 && theta[i] == 0) {
      faux_pos <- faux_pos + 1
    }
    
    if (decision[i] == 0 && theta[i] == 1) {
      faux_neg <- faux_neg + 1
    }
  }
  
  cat("Nb vrais sites, Résultats, Nb faux pos, Nb faux neg\n")
  cat(sum(theta), sum(decision), faux_pos, faux_neg, sum(decision) - faux_pos + faux_neg, "\n")
}
methode2(10,0.05,0,1,5000,0.05)




#################
benj<-function(n,alpha,mu,sd,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pval<-rep(0,J)
  decision<-rep(0,J)
  theta<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sd) #vecteur
    theta[i]<-rbinom(1,1,p)
    Y<-theta[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pval[i]<-1-pt(T[i],n-1)
  }
  pval_s<-sort(pval)
  for(j in (J-1):1){
    if(pval_s[j]*J/j<alpha){
      decision[j]=1
      if(theta[order(pval)[j]]==0){
        faux_pos<-faux_pos+1
      }
    }
    else{
      decision[j]=0
      if(theta[order(pval)[j]]==1){
        faux_neg<-faux_neg+1
      }
    }
  }
  print("Nb vrais sites, Resultats, Nb faux pos, Nb faux neg")
  print(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision) - faux_pos + faux_neg))
}
benj(10,0.05,0,1,5000,0.05)
