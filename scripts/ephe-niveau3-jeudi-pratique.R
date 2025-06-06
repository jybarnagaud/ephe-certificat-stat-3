#---------------------------------------------------------------------#
### Formation "Analyse de donn�es pour les �cologues avanc�" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# TD du jeudi : pratiquer la RLQ, interpr�ter, faire des graphiques 
# lisibles.
# NB : seul le premier exemple est comment� - les codes sont tous les 
# m�mes. A noter que les �tapes interm�diaires (analyses � 1 tableau, 
#  coinertie) ne sont cod�es que dans le premier cas, mais elles doivent
# �tre faites syst�matiquement afin de progresser graduellement dans 
# l'interpr�tation.
#---------------------------------------------------------------------#

rm(list=ls())
library(ade4)
setwd("C:/Users/jeany/OneDrive/Documents/EPHE_enseignement/Statistiques/2024/niveau3/enseignant/donnees")

#-----------------------------------------------------------------#
#### Cas pratique 1 : communaut�s de papillons dans les Landes ####
#-----------------------------------------------------------------#

# index des noms d'esp�ces : index_lepidos.csv

# donn�es

traits=read.table("table_traits_lepidos.txt",header=T,sep="\t")
rhopa=read.table("table_especes_lepidos.txt",header=T,sep="\t")
habitat=read.table("table_habitat_lepidos.txt",header=T,sep="\t",row.names=1)

# homog�n�iser les esp�ces 

list.sp = intersect(rownames(traits),colnames(rhopa))
rhopa.sub = rhopa[,list.sp]
traits.sub = traits[list.sp,]
traits.sub$polyphagie = factor(traits.sub$polyphagie)
traits.sub$type_host_plant = factor(traits.sub$type_host_plant)

# exploration matrice par matrice : faire les projections graphiques de chaque table

afc.L=dudi.coa(rhopa.sub,scannf=F,nf=5)
pca.R=dudi.pca(habitat,scannf=F,nf=3) 
hs.Q=dudi.hillsmith(traits.sub,scannf=F,nf=4)

# exploration de la coinertie

pca.Rw=dudi.pca(habitat,scannf=F,nf=3,row.w = afc.L$lw) 
coi.rhopa = coinertia(afc.L,pca.Rw,scannf=F,nf=2)

# analyse RLQ --> on garde 2 axes, l'un des deux domine nettement

hs.Qw=dudi.hillsmith(traits.sub,scannf=F,nf=4,row.w=afc.L$cw)

RLQ=rlq(pca.Rw,afc.L,hs.Q,scannf=F,nf=2)


# graphique global (pas tr�s lisible)

plot(RLQ)

# figure par figure

s.arrow(RLQ$l1) # variables de la matrice R dans l'espace RLQ
s.arrow(RLQ$c1) # traits de la matrice Q dans l'espace RLQ
s.label(RLQ$lQ, boxes = FALSE) # esp�ces dans la matrice RLQ
s.label(RLQ$lR, boxes = FALSE) # sites dans la matrice RLQ

# on peut tester la significativit� globale des relations traits - environnement, � partir de l'inertie de la RLQ

testrlq=randtest(RLQ, modeltype = 6, nrepet = 999)
plot(testrlq)

# on essaie un autre test

testrlq2=randtest(RLQ, modeltype = 5, nrepet = 999) # on permute les esp�ces avant de permuter les sites
plot(testrlq2)

## analyse fourth-corner
# on fait 99 permutations, on peut aussi faire 999 

fourth.rhopalo=fourthcorner(data.frame(habitat),data.frame(rhopa.sub),data.frame(traits.sub),modeltype = 6, p.adjust.method.G = "none",p.adjust.method.D = "none", nrepet = 99)
plot(fourth.rhopalo, alpha = 0.05, stat = "D2") # si un trait est cat�gorique, on mesure s�par�ment les associations entre chaque cat�gorie et les variables environnementales
plot(fourth.rhopalo, alpha = 0.05, stat = "G") # "", on mesure l'association entre la variable environnementale et toutes les cat�gories d'un coup (F-test, comme dans une ANOVA)
plot(fourth.rhopalo, alpha = 0.05, stat = "D") # "", on mesure l'association entre la variable environnementale et chaque cat�gorie s�par�ment � partir d'un test d'homog�n�it� intra-groupe

# Et si on applique une correction ("fdr") pour tests multiples

fourth.rhopalo.correct=fourthcorner(data.frame(habitat),data.frame(rhopa.sub),data.frame(traits.sub),modeltype = 6, p.adjust.method.G = "fdr",p.adjust.method.D = "fdr", nrepet = 99)

plot(fourth.rhopalo.correct, alpha = 0.05, stat = "D2")
plot(fourth.rhopalo.correct, alpha = 0.05, stat = "G") 
plot(fourth.rhopalo.correct, alpha = 0.05, stat = "D") 

#----------------------------------------------#
#### Cas pratique 2 : reptiles en Occitanie ####
#----------------------------------------------#

# donn�es

L=read.table("Reptiles_matrice_sites_especes.txt",header=T,sep="\t",row.names=1)
R=read.table("Reptiles_matrice_habitats.txt",header=T,sep="\t",row.names=1)
Q=read.table("Reptiles_matrice_especes_traits.txt",header=T,sep="\t",row.names=1)
Q$reprod_mode = factor(Q$reprod_mode)

# analyse sur la matrice L

Lcoa=dudi.coa(L,scannf=F,nf=2)
Qhs=dudi.hillsmith(Q,row.w =Lcoa$cw,scannf=F,nf=2)
Rpca=dudi.pca(R,row.w = Lcoa$lw,scannf=F,nf=2)

# RLQ
rlq=rlq(Rpca,Lcoa, Qhs)

# test par permutations

testrlq=randtest(RLQ, modeltype = 6, nrepet = 999)
plot(testrlq)

# graphiques

plot(rlq)

# fen�tre par fen�tre

s.arrow(rlq$l1)
s.arrow(rlq$c1)
s.label(rlq$lQ, boxes = FALSE)

# voir les scores

summary(rlq)

# fourth corner

nrepet=999
fc.rept=fourthcorner(R,L,Q, modeltype = 6, p.adjust.method.G = "none",p.adjust.method.D = "none", nrepet = nrepet)
plot(fc.rept, alpha = 0.05, stat = "D2")

fc.adj=p.adjust.4thcorner(fc.rept,p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
plot(fc.adj, alpha = 0.05, stat = "D2")

test.qaxes=fourthcorner.rlq(rlq,modeltype=6,typetest="Q.axes")
print(test.qaxes,stat="D")

#----------------------------------------------------#
#### Cas pratique 3 : oiseaux de Nouvelle Z�lande ####
#----------------------------------------------------#

# index des noms d'esp�ces : NZ_index_especes

# donn�es (NB : les sites sont d�j� ordonn�s dans les matrices habitats et sites x esp�ces)

NZ_habitats = read.table("NZ_habitats.txt",header=T,sep="\t")
NZ_sitsp = read.table("NZ_oiseaux.txt",header=T,sep="\t")
NZ_traits = read.table("NZ_traits.txt",header=T,sep="\t",stringsAsFactors = T)
NZ.index = read.table("NZ_index_especes.txt",header=T,sep="\t")

# analyse sur la matrice L

Lcoa=dudi.coa(NZ_sitsp,scannf=F,nf=2)
Qhs=dudi.hillsmith(NZ_traits,row.w =Lcoa$cw,scannf=F,nf=2)
Rpca=dudi.pca(NZ_habitats,row.w = Lcoa$lw,scannf=F,nf=2)

# RLQ

rlq=rlq(Rpca,Lcoa, Qhs,scannf=F,nf=2)

# test par permutations

testrlq=randtest(rlq, modeltype = 6, nrepet = 999)
plot(testrlq)

# graphiques

plot(rlq)

# fen�tre par fen�tre

s.arrow(rlq$l1)
s.arrow(rlq$c1)
s.label(rlq$lQ, boxes = FALSE)

# voir les scores

summary(rlq)

# fourth corner

nrepet=999
fc.nz=fourthcorner(NZ_habitats,NZ_sitsp,NZ_traits, modeltype = 6, p.adjust.method.G = "none",p.adjust.method.D = "none", nrepet = nrepet)
plot(fc.nz, alpha = 0.05, stat = "D2")

fc.nz.adj=p.adjust.4thcorner(fc.nz,p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
plot(fc.nz.adj, alpha = 0.05, stat = "D2")

test.qaxes=fourthcorner.rlq(rlq,modeltype=6,typetest="Q.axes")
print(test.qaxes,stat="D")

# diff�rencier exotiques et natives dans l'espace de RLQ

rownames(NZ.index) = NZ.index$species.acronym
NZ.index = NZ.index[rownames(rlq$lQ),]
s.class(rlq$lQ, fac=factor(NZ.index$Status.in.NZ),col=c("darkred","steelblue","orange"))

#----------------------------------------------#
#### Cas pratique 4 : oiseaux de La R�union ####
#----------------------------------------------#

# donn�es

acronymes.sp = read.table("oiseaux_reunion_acronymes.txt",header=T)
traits.reu = read.table("oiseaux_lareunion_traits.txt",header=T,row.names=1,sep="\t")
traits.reu$Source = factor(traits.reu$Source)
traits.reu$Diet = factor(traits.reu$Diet)
traits.reu$Vegetation_strata = factor(traits.reu$Vegetation_strata)
habitats.reu = read.table("oiseaux_la_reunion_habitats.txt",header=T,row.names=1)
sites.sp.reu = read.table("oiseaux_la_reunion.txt",header=T,sep="\t")

# analyse sur la matrice L

Lcoa=dudi.coa(sites.sp.reu,scannf=F,nf=2)
Qhs=dudi.hillsmith(traits.reu,row.w =Lcoa$cw,scannf=F,nf=2)
Rpca=dudi.pca(habitats.reu,row.w = Lcoa$lw,scannf=F,nf=2)

# RLQ

rlq=rlq(Rpca,Lcoa, Qhs,scannf=F,nf=2)

# test par permutations

testrlq=randtest(rlq, modeltype = 6, nrepet = 999)
plot(testrlq)

# graphiques

plot(rlq)

# fen�tre par fen�tre

s.arrow(rlq$l1)
s.arrow(rlq$c1)
s.label(rlq$lQ, boxes = FALSE)

# voir les scores

summary(rlq)

# fourth corner

nrepet=999
fc.reunion=fourthcorner(habitats.reu,sites.sp.reu,traits.reu, modeltype = 6, p.adjust.method.G = "none",p.adjust.method.D = "none", nrepet = nrepet)
plot(fc.reunion, alpha = 0.05, stat = "D2")

fc.reunion.adj=p.adjust.4thcorner(fc.reunion,p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")
plot(fc.reunion.adj, alpha = 0.05, stat = "D2")

test.qaxes=fourthcorner.rlq(rlq,modeltype=6,typetest="Q.axes")
print(test.qaxes,stat="D")

# diff�rencier exotiques et natives dans l'espace de RLQ

rownames(acronymes.sp) = acronymes.sp$acronym
acronymes.sp = acronymes.sp[rownames(rlq$lQ),]
s.class(rlq$lQ, fac=factor(acronymes.sp$Status),col=c("darkred","steelblue"))

#-----------------------------------------------------------#
#### Ce qu'il faut comprendre de ces diff�rents exemples ####
#-----------------------------------------------------------#

# - La RLQ est descriptive et ne fait pas d'estimation de param�tres : les relations traits-environnements existent, m�me si le test par permutations les trouve non significatives
# - Le test de permutation dit � quel point les relations traits - environnement observ�es pourraient �merger sous un mod�le nul. L'int�r�t du mod�le en deux �tapes est qu'il permet 
#   d'identifier si le patron trait-environnement observ� est plut�t d� aux r�partitions d'esp�ces sur les gradients environnementaux, ou s'il y a plus de structures qu'attendu avec ces associations
# - attention aux tests ajust�s, qui peuvent �tre trop conservatifs
