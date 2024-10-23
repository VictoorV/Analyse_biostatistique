Stud_BH<-function(n,alpha,mu,sigma,J,p){
  faux_pos<-0
  faux_neg<-0
  T<-rep(0,J)   #Stat de test 
  alpha_sid<-1-(1-alpha)**(1/J)  # Le nouveau niveau pour la correction de Sidak
  pvalm<-rep(0,J) 
  decision<-rep(0,J)
  thetam<-rep(0,J)
  for(i in 1:J){ # boucle pour créer les échatillons
    bruit<-rnorm(n,mu,sigma**2) #vecteur
    thetam[i]<-rbinom(1,1,p)
    Y<-thetam[i]+bruit
    T[i]<-sqrt(n)*(mean(Y)-0)/sqrt(var(Y))
    pvalm[i]<-1-pt(T[i],n-1)}
  tableau2 <- data.frame(thetam,pvalm) # tableau avec p-valeur et valeur de théta
  
  pval2<-sort(pvalm) # On ré ordonne les p-val pour les correction BH et HB
  pval_BH<-rep(0,J)
  pval_HB<-rep(0,J)
  Nbr_accept<-0
  Nbr_accept_sidak<-0
  Nbr_accept_HB<-0
  
  attach(tableau2)
  tableau <- tableau2[order(pvalm),] # on réordonne le tableau en fonction des p-val
  detach(tableau2)
  
  for (i in 1:J){ # on créer les p-val ajusté pour BH et HB
    pval_BH[i]<-pval2[i]*J/i
    pval_HB[i]<-pval2[i]*(J-i+1)}
  tableau$pval_BH<-pval_BH
  tableau$pval_HB<-pval_HB
  
  for ( j in 1:J){    # Recherche du seuil pour BH
    if(pval_BH[j]>=0.05){
      #      print(j)
      Nbr_accept<-j-1
      break}}
  for ( j in 1:J){   #On effectue les test pour Sidak et HB
    if(tableau[j,2]<=alpha_sid){
      Nbr_accept_sidak<-Nbr_accept_sidak+1} 
    if(tableau[j,4]<=alpha){
      Nbr_accept_HB<-Nbr_accept_HB+1}}
  fo_pos_sid<-0
  fo_pos_BH<-0
  fo_pos_HB<-0
  
  fo_neg_sid<-0
  fo_neg_BH<-0
  fo_neg_HB<-0
  
  for (i in 1:J){ # Calcul des faux positif et négatif pour chaque correction
    if(tableau[i,1]==0){
      if(tableau[i,2]<alpha_sid){fo_pos_sid<-fo_pos_sid+1}
      if(i<=Nbr_accept){fo_pos_BH<-fo_pos_BH+1}
      if(tableau[i,4]<alpha){fo_pos_HB<-fo_pos_HB+1}}
    if(tableau[i,1]==1){
      if(tableau[i,2]>=alpha_sid){fo_neg_sid<-fo_neg_sid+1}
      if(i>Nbr_accept){fo_neg_BH<-fo_neg_BH+1}
      if(tableau[i,4]>=alpha){fo_neg_HB<-fo_neg_HB+1}}}
  
  Sidak<-c(Nbr_accept_sidak,fo_pos_sid,fo_neg_sid)
  BH<-c(Nbr_accept,fo_pos_BH,fo_neg_BH)
  HB<-c(Nbr_accept_HB,fo_pos_HB,fo_neg_HB)
  Vrai<-c(sum(thetam),0,0)
  retour<-data.frame(Vrai,Sidak,BH,HB) # on renvoi un tableau 
  
  return(retour)
  
}
Stud_BH(10,0.05,0,1,5000,0.05)
