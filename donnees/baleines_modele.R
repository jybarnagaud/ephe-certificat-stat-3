model{

#------------------#	
# a prioris (priors)
#------------------#

# sur la variance de la variable de r�ponse
tau<-1/(sig*sig)
sig~dunif(0,100)
	
# sur les effets des variables continues
alpha~dnorm(0,0.0001) # intercept
beta~dnorm(0,0.0001) # pente de l'effet nb_indiv

# sur les effets des variables facteur
gamma[1]<-0 # on fixe le premier niveau de facteur � 0, il servira de r�f�rence !! attention � bien rep�rer la correspondance entre niveaux de facteurs et indices

# note sur les niveaux fix�s: ici, on fixe un niveau parce qu'il y a d�j� un intercept (alpha) dans le mod�le. S'il n'y en avait pas, il serait inutile de fixer ce niveau

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
	# on fait tourner le mod�le sur toutes les donn�es
time_z12[i]~dnorm(mu[i],tau)

#------#
# mod�le
#------#

mu[i]<-alpha+beta*Lnb_indiv[i]+gamma[chemical[i]]+delta[place[i]]

#-------#
# r�sidus
#-------#

residual[i]<-time_z12[i]-mu[i] # r�sidus
predicted[i]<-mu[i]  # valeurs pr�dites
sq[i]<-pow(residual[i],2) # r�sidus au carr� pour les donn�es observ�es

#------------------#
# donn�es r�pliqu�es
#------------------#

time_z12.new[i]~dnorm(mu[i],sig) # on cr�e un nouveau jeu de donn�es
sq.new[i]<-pow(time_z12.new[i]-predicted[i],2) # r�sidus du nouveau jeu de donn�es

} # fermer la boucle i

}
