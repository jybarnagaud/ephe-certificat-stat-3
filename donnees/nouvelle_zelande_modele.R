model{

# priors des covariables
det.alpha~dnorm(0,0.0001)
beta~dnorm(0,0.0001)
delta~dnorm(0,0.0001)
muyear~dnorm(0,0.0001) # prior de l'effet al�atoire ann�es

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

# prior hi�rarchique sur alpha: tous les sites faits une m�me ann�e ont un prior commun (<=> effet al�atoire ann�e)
alpha[i]~dnorm(mualpha[year[i]],taualpha)
	
# hypoth�se: l'abondance est fonction du site, de la proportion de for�t native sur le site et de l'altitude sur le site

# vraisemblance de la couche de proba de d�tection
for(j in 1:nrepl){
obs[i,j] ~ dbin(p[i],N[i])
}
# obs correspond aux donn�es brutes

# mod�le pour la couche de d�tection
logit(p[i])<-det.alpha+beta.det[year[i]]*julian_date[i]

# hypoth�se: la proba de d�tection a une valeur de base modul�e par la date, et cette modulation change entre les deux ann�es de suivi
		
#----------------#
### ajustement ###
#----------------#

# valeurs attendues de obs
for(j in 1:nrepl){
eval[i,j]<-N[i]*p[i]
E[i,j]<-pow((obs[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
	
# donn�es r�pliqu�es
obs.new[i,j]~dbin(p[i],N[i])
E.new[i,j]<-pow((obs.new[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
	
}	# i
}	# j

fit<-sum(E[,])
fit.new<-sum(E.new[,])

#--------------------------------------------------------#
### quantit� d�riv�e: probabilit� de d�tection moyenne ###
#--------------------------------------------------------#

pmean<-mean(p[])

}