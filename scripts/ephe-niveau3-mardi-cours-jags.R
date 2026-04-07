model{

  
# vraisemblance
  
  for(i in 1 : n.parcelles){
    
consommation[i] ~ dbern(psi[i])
  
  logit(psi[i]) <- alpha + beta*div_spe_veg[i]
  
  }
  
# priors
  
alpha ~ dnorm(0,0.5)
beta ~ dnorm(0,0.5)
  
  
}