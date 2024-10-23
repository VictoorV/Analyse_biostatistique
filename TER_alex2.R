# Fonction methode de plusieurs tests multiples
methode1<-function(n,alpha,mu,sigma,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pval<-rep(0,J)
  decision<-rep(0,J)
  theta<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sigma**2) #vecteur
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
  #print("Nb vrai sites, Resultats, Nb faux pos, Nb faux neg")
  return(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision)-faux_pos+faux_neg))
}

J<-5000
alpha<-0.05
mu<-0
sigma<-0.5
p<-0.05
Bonf<-(alpha/J)
methode1(10,alpha,mu,sigma**2,J,p)
J*alpha*(1-p)
methode1(10,Bonf,mu,sigma**2,J,p)
##sidak<-1-(1-alpha/J)**J            verifier cette formule qui semble pas marcher, surement que c sidak2
sidak2<- 1- (1-alpha)**(1/J)
methode1(10,sidak2,mu,sigma**2,J,p)
quant_Bonf<-qt(Bonf,n-1)
quant_Bonf

### Aucun faux neg avec sidak et a peu prés autant de faux pos que de vrai pos

##########################################"

n<-10
J<-5000

Stud_bonf<-function(n,J,alpha){
  T<-rep(0,J)
  resultat=rep(0,J)
  theta<-rep(0,J)
  epsilon<-rep(0,n*J)
  s<-0
  for(i in 1:J){
    a=rnorm(n,0,1)
    b=rbinom(1,1,0.05)
    theta[i]<-b
    epsilon[i]<-a
    T[i]<-sqrt(n)*(mean(a+b))/sqrt(var(a+b))
    resultat[i]<-ifelse(T[i]<=qt(n-1,alpha/J),0,1)
    s<-ifelse(resultat[i]>b,s+1,s)
  }
  return(c("Nbr faux positif:",sum(s),"Nbr vrai site:",sum(theta),"Nbr de ccl H1",sum(resultat)))
}
s # les faux positif
sum(theta) #le nombre de vrai site
sum(resultat) #le nombre de fois ou on conclut H1 a tord ou a raison
theta
Stud_bonf(10,5000,0.05)

########################################################

Stud_BH<-function(n,alpha,mu,sigma,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pvalm<-rep(0,J)
  decision<-rep(0,J)
  thetam<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sigma**2) #vecteur
    thetam[i]<-rbinom(1,1,p)
    Y<-thetam[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pvalm[i]<-1-pt(T[i],n-1)}
  tableau2 <- data.frame(thetam,pvalm)
  
  pval2<-sort(pvalm)
  pval_BH<-rep(0,J)
  pval_HB<-rep(0,J)
  Nbr_accept<-0
  Nbr_accept_sidak<-0
  Nbr_accept_sidak2<-0
  Nbr_accept_HB<-0
  
  
  attach(tableau2)
  tableau <- tableau2[order(pvalm),]
  detach(tableau2)
  
  alpha_sid<-1-(1-alpha)**(1/J)
  alpha_sid2<-alpha/(1-(1-alpha)**(1/J))
  
  for (i in 1:J){
    pval_BH[i]<-pval2[i]*J/i
    pval_HB[i]<-pval2[i]*(J-i+1)}
  tableau$pval_BH<-pval_BH
  tableau$pval_HB<-pval_HB
  
  for ( j in 1:J){
    if(pval_BH[j]>=0.05){
      #      print(j)
      Nbr_accept<-j-1
      break}}
  for ( j in 1:J){
    if(tableau[j,2]<=alpha_sid){
      Nbr_accept_sidak<-Nbr_accept_sidak+1} 
    if(tableau[j,2]<=alpha_sid2){
      Nbr_accept_sidak2<-Nbr_accept_sidak2+1}
    if(tableau[j,4]<=alpha){
      Nbr_accept_HB<-Nbr_accept_HB+1}
  }
  
  fo_pos_sid<-0
  fo_pos_sid2<-0
  fo_pos_BH<-0
  fo_pos_HB<-0
  
  fo_neg_sid<-0
  fo_neg_sid2<-0
  fo_neg_BH<-0
  fo_neg_HB<-0
  
  for (i in 1:J){
    if(tableau[i,1]==0){
      if(tableau[i,2]<alpha_sid){fo_pos_sid<-fo_pos_sid+1}
      if(tableau[i,2]<alpha_sid2){fo_pos_sid2<-fo_pos_sid2+1}
      if(i<=Nbr_accept){fo_pos_BH<-fo_pos_BH+1}
      if(tableau[i,4]<alpha){fo_pos_HB<-fo_pos_HB+1}
    }
    if(tableau[i,1]==1){
      if(tableau[i,2]>=alpha_sid){fo_neg_sid<-fo_neg_sid+1}
      if(tableau[i,2]>=alpha_sid2){fo_neg_sid2<-fo_neg_sid2+1}
      if(i>Nbr_accept){fo_neg_BH<-fo_neg_BH+1}
      if(tableau[i,3]>=alpha){fo_neg_HB<-fo_neg_HB+1}
    }}
  
  #print("Nb vrai sites, Resultats, Nb faux pos, Nb faux neg")
  Sidak<-c(Nbr_accept_sidak,fo_pos_sid,fo_neg_sid)
  Sidak2<-c(Nbr_accept_sidak2,fo_pos_sid2,fo_neg_sid2)
  BH<-c(Nbr_accept,fo_pos_BH,fo_neg_BH)
  HB<-c(Nbr_accept_HB,fo_pos_HB,fo_neg_HB)
  Vrai<-c(sum(thetam),0,0)
  retour<-data.frame(Vrai,Sidak,Sidak2,BH,HB)
  return(retour)
  #  return(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision)-faux_pos+faux_neg))
}

####################    C EST ICI      #############################"
####################    C EST ICI      #############################"
####################    C EST ICI      #############################"
####################    C EST ICI      #############################"
####################    C EST ICI      #############################"

Stud_BH(10,0.05,0,0.75,5000,0.05)

## Lecture résultat: ligne 1: Nbr accepté
##                   ligne 2: Nbr Faux positif
##                   ligne 3 : Nbr faux négatif

##                   colonne 1: Les vrai théta
##                   colonne 2: Sidak-Bonferonni
##                   colonne 4: Benjamini-Hochberg ( le plus efficace)
##                   colonne 5: Holm-Bonferonni

## On noteras que sidak2 (colonne 3) est éclaté au sol ( méthode sidak classique)

#######################################

Stud_HB<-function(n,alpha,mu,sigma,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)
  pvalm<-rep(0,J)
  decision<-rep(0,J)
  thetam<-rep(0,J)
  for(i in 1:J){
    bruit<-rnorm(n,mu,sigma**2) #vecteur
    thetam[i]<-rbinom(1,1,p)
    Y<-thetam[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pvalm[i]<-1-pt(T[i],n-1)}
  tableau2 <- data.frame(thetam,pvalm)
  
  pval2<-sort(pvalm)
  pval_HB<-rep(0,J)
  pval_BH<-rep(0,J)
  Nbr_accept_HB<-0
  Nbr_accept_BH<-0
  Nbr_accept_SB<-0
  Nbr_accept_sid<-0
  
  attach(tableau2)
  tableau <- tableau2[order(pvalm),]
  detach(tableau2)
  
  alpha_SB<-1-(1-alpha)**(1/J)
  alpha_sid<-
    
    for (i in 1:J){
      pval_HB[i]<-alpha/(J-i+1)}
  tableau$pval_HB<-pval_HB
  
  fo_pos_sid<-0
  fo_pos_BH<-0
  fo_pos_HB<-0
  fo_pos_SB<-0
  
  fo_neg_sid<-0
  fo_neg_BH<-0
  fo_neg_HB<-0
  fo_neg_SB<-0
  
  
  for (i in 1:J){
    if(tableau[i,1]==0){
      if(tableau[i,2]<alpha_sid){fo_pos_sid<-fo_pos_sid+1}
      if(tableau[i,3]<alpha){fo_pos_BH<-fo_pos_BH+1}
      
    }
    if(tableau[i,1]==1){
      if(tableau[i,2]>=alpha_sid){fo_neg_sid<-fo_neg_sid+1}
      if(tableau[i,2]>=tableau[i,3]){fo_neg_HB<-fo_neg_HB+1}
      
    }}
  for ( j in 1:J){
    if(pval_BH[j]>=0.05){
      #      print(j)
      Nbr_accept<-j
      break}}
  for ( j in 1:J){
    if(pval2[j]>=1-(1-alpha)**(1/J)){
      #      print(j)
      Nbr_accept_sidak<-j
      break
    }  
  }
  
  #print("Nb vrai sites, Resultats, Nb faux pos, Nb faux neg")
  return(c(Nbr_accept-1,fo_pos_BH,fo_neg_BH,Nbr_accept_sidak-1,fo_pos_sid,fo_neg_sid))
  #  return(c(sum(theta),sum(decision),faux_pos,faux_neg,sum(decision)-faux_pos+faux_neg))
}
