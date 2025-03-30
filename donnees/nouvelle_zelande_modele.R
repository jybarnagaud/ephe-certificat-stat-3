model{

# priors des covariables
det.alpha~dnorm(0,0.0001)
beta~dnorm(0,0.0001)
delta~dnorm(0,0.0001)
muyear~dnorm(0,0.0001) # prior de l'effet aléatoire années

taualpha<-1/(sigalpha*sigalpha)
sigalpha~dunif(0,100)

tauyear<-1/(sigyear*sigyear)
sigyear~dunif(0,100)

for (k in 1:nyears){
	beta.det[k]~dnorm(0,0.0001)
	mualpha[k]~dnorm(muyear,tauyear)
}

for(i in 1:npoints){
	
# vraisemblance de la couche de process
	N[i] ~ dpois(esp[i])
	log(esp[i])<-alpha[i]+beta*PROPNATFOR[i]+delta*altitude[i]

# prior hiérarchique sur alpha: tous les sites faits une même année ont un prior commun (<=> effet aléatoire année)
alpha[i]~dnorm(mualpha[year[i]],taualpha)
	
# hypothèse: l'abondance est fonction du site, de la proportion de forêt native sur le site et de l'altitude sur le site

# vraisemblance de la couche de proba de détection
for(j in 1:nrepl){
obs[i,j] ~ dbin(p[i],N[i])
}
# obs correspond aux données brutes

# modèle pour la couche de détection
logit(p[i])<-det.alpha+beta.det[year[i]]*julian_date[i]

# hypothèse: la proba de détection a une valeur de base modulée par la date, et cette modulation change entre les deux années de suivi
		
#----------------#
### ajustement ###
#----------------#

# valeurs attendues de obs
for(j in 1:nrepl){
eval[i,j]<-N[i]*p[i]
E[i,j]<-pow((obs[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
	
# données répliquées
obs.new[i,j]~dbin(p[i],N[i])
E.new[i,j]<-pow((obs.new[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
	
}	# i
}	# j

fit<-sum(E[,])
fit.new<-sum(E.new[,])

#--------------------------------------------------------#
### quantité dérivée: probabilité de détection moyenne ###
#--------------------------------------------------------#

pmean<-mean(p[])

}