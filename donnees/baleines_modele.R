model{

#------------------#	
# a prioris (priors)
#------------------#

# sur la variance de la variable de réponse
tau<-1/(sig*sig)
sig~dunif(0,100)
	
# sur les effets des variables continues
alpha~dnorm(0,0.0001) # intercept
beta~dnorm(0,0.0001) # pente de l'effet nb_indiv

# sur les effets des variables facteur
gamma[1]<-0 # on fixe le premier niveau de facteur à 0, il servira de référence !! attention à bien repérer la correspondance entre niveaux de facteurs et indices

# note sur les niveaux fixés: ici, on fixe un niveau parce qu'il y a déjà un intercept (alpha) dans le modèle. S'il n'y en avait pas, il serait inutile de fixer ce niveau

for(k in 2:nchemical) {
gamma[k]~dnorm(0,0.0001) # on estime les autres niveaux	
}

delta[1]<-0

for(n in 2:nplace){
delta[n]~dnorm(0,0.0001)
}

#-------------#
# vraisemblance
#-------------#

for(i in 1:ndata){
	# on fait tourner le modèle sur toutes les données
time_z12[i]~dnorm(mu[i],tau)

#------#
# modèle
#------#

mu[i]<-alpha+beta*Lnb_indiv[i]+gamma[chemical[i]]+delta[place[i]]

#-------#
# résidus
#-------#

residual[i]<-time_z12[i]-mu[i] # résidus
predicted[i]<-mu[i]  # valeurs prédites
sq[i]<-pow(residual[i],2) # résidus au carré pour les données observées

#------------------#
# données répliquées
#------------------#

time_z12.new[i]~dnorm(mu[i],sig) # on crée un nouveau jeu de données
sq.new[i]<-pow(time_z12.new[i]-predicted[i],2) # résidus du nouveau jeu de données

} # fermer la boucle i

}
