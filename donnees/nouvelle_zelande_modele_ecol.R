model{
  
  # priors des covariables
  beta~dnorm(0,0.0001)
  delta~dnorm(0,0.0001)
  
  # priors de l'effet années
  for (k in 1:nyears){
    alpha[k]~dnorm(0,1000)
  }
  
  for(i in 1:npoints){
    
    # vraisemblance de la couche de process
    N[i] ~ dpois(esp[i])
    log(esp[i])<-alpha[year[i]]+beta*PROPNATFOR[i]+delta*altitude[i]
  
    # valeurs répliquées de N (pour l'ajustement)
    N.rep[i] ~ dpois(esp[i])
    
    E[i]<-pow((N[i]-esp[i]),2)/(esp[i]+0.5)
    E.rep[i]<-pow((N.rep[i]-esp[i]),2)/(esp[i]+0.5)
    
    # résidus
    r.naif[i] <- N[i] - esp[i]
  
      
      
    }	# i

}