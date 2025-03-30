#----------------------------------------------------------------------------#
### Formation "Certificat en analyse de données pour les écologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# ce script permet la réplication des exemples du cours sur les modèles mixtes
#----------------------------------------------------------------------------#

rm(list=ls())
setwd("F:/certificat_2023/niveau3/supports/donnees_final")

# packages utiles
library(lme4)
library(ggeffects)
library(patchwork)
library(DHARMa)
library(maptools)
library(questionr)
library(ggplot2)
library(viridis)
library(sf)
library(sjPlot)
library(cowplot)
library(lme4)
library(nlme)
library(mgcv)
library(ggeffects)
library(MuMIn)
library(performance)
library(lmerTest)
library(glmmTMB)
library(see)
library(visreg)
library(reshape2)

#--------------------------------------------------------#
#### simulation : pourquoi la non indépendance résiduelle
# doit être prise au sérieux ####
#--------------------------------------------------------#

# 20 blocs (= strates)
blocs = paste("B",1:20,sep="")

# moyenne de l'intercept de chaque bloc
a = rnorm(20,5,10)

# pente unique pour tous les blocs
b = -5

# une variable explicative mesurée de la même manière dans tous les blocs, avec
# 30 points par bloc
x = rnorm(30,0,1)

# un résidu non structuré avec une certaine variance 
eps = rnorm(30,0,1)

# on génère la variable de réponse- dans chaque bloc
var.rep = NULL
for(i in 1:length(a)){
y = a[i]+b*x+eps
var.rep = c(var.rep,y)
}

data.sim = data.frame(reponse = var.rep,
           blocs = rep(blocs,each=30),
           predicteur = rep(x,times=length(a)))

# données "observées"
library(ggplot2)
ggplot(data.sim)+
  aes(x = blocs, y = reponse)+
      geom_boxplot()

ggplot(data.sim)+
  aes(x = predicteur, y = reponse,col=blocs)+
  geom_point()+
    labs(x = "Profondeur (centrée réduite)",y = "Concentration planctonique",title = "simulation")

# on ajuste un modèle linéaire sans tenir compte des blocs
m.naif = lm(reponse~predicteur,data=data.sim)

# diagnostic : résidus corrects
par(mfrow=c(2,2)) ; plot(m.naif)

# la prédiction parait correcte
library(ggeffects)
plot(ggpredict(m.naif),add.data=T)

# variabilité des résidus par blocs : sous l'hypothèse du modèle linéaire, cette
# variation ne devrait pas exister
boxplot(residuals(m.naif)~data.sim$blocs)

# impact sur les coefficients : 
coef(m.naif)
confint(m.naif)

# vraie valeur populationnelle de l'intercept : l'estimation est peu décalée, 
# mais l'intervalle de confiance de l'estimation n'est pas assez large
moy.sim.a = mean(a)
ic.sim.a = paste(round(mean(a) - 1.96 * (sd(a)/sqrt(length(a))),2), ";",
round(mean(a) + 1.96 * (sd(a)/sqrt(length(a))),2))


# en revanche, la pente est bonne (autour de -5) : de toute façon elle n'était pas
#  supposée varier d'un bloc à l'autre

# on tient compte de l'effet strates
m.strate = lm(reponse~predicteur+blocs,data=data.sim)
confint(m.strate)

# est-ce que l'intervalle de confiance de la strate B1 est bon? 
data.b1 = subset(data.sim,blocs == "B1")

moy.sim = mean(data.b1$reponse)
sd.sim = sd(data.b1$reponse)

ic.sim.b1 = paste(round(moy.sim - 1.96 * (sd.sim/sqrt(nrow(data.b1))),2), ";",
               round(moy.sim + 1.96 * (sd.sim/sqrt(nrow(data.b1))),2))

# l'estimation est à nouveau un peu étroite
# Il faut réfléchir un instant à la manière dont le jeu de données est construit : 
# on a simulé des valeurs populationnelles pour a et b, donc il faut que notre 
# inférence se situe aussi à échelle populationnelle

library(lme4) 
m.rand = lmer(reponse~predicteur+(1|blocs),data=data.sim)
confint(m.rand)

# cette fois l'intervalle de confiance de a est quasiment le bon

# maintenant, la pente b va varier d'une strate à l'autre
b2 =  rnorm(20,-5,5)

var.rep2 = NULL
for(i in 1:length(a)){
  y = a[i]+b2[i]*x+eps
  var.rep2 = c(var.rep2,y)
}

data.sim$reponse2 = var.rep

# on relance le modèle mixte avec un effet aléatoire sur l'intercept(= on ne tient 
# pas compte du fait que l'effet de la variable explicative change d'un bloc à l'autre)
m.rand2 = lmer(reponse2~predicteur+(1|blocs),data=data.sim)
confint(m.rand2)

# le vrai intervalle de confiance de b

moy.sim.b2 = mean(b2)
ic.sim.b2 = paste(round(mean(b2) - 1.96 * (sd(b2)/sqrt(length(b2))),2), ";",
                 round(mean(b2) + 1.96 * (sd(b2)/sqrt(length(b2))),2))

moy.sim.b2
ic.sim.b2

m.rand3 = lmer(reponse2~predicteur+(predicteur|blocs),data=data.sim)
confint(m.rand3)

#-----------------------------------------------------------#
#### Le modèle linéaire mixte : chenille processionnaire ####
#-----------------------------------------------------------#

# données
chenilles <- read.csv2("chenille_processionnaire.csv",row.names=1, encoding="UTF-8")
chenilles$prop_attaq <- chenilles$nbpinsattac / chenilles$nbpins

## exploration des données

# carte des régions
regions <- st_read("regions.geojson")

ggplot(regions)+
  geom_sf(aes(fill=jointure_g))+
  geom_point(data=chenilles,mapping=aes(x=longitude,y=latitude,colour=region),size=2)+
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) + 
  labs(color="Régions")+ 
  theme_classic()

##### proportion d'arbres infestés par an

ggplot(chenilles) +
  aes(x = factor(annee), 
      y = prop_attaq, 
      group = placette, 
      color = region) +
  geom_line(show.legend = FALSE) +
  labs(x = "Années", y = "Proportion de pins infestés") +
  facet_wrap(~region, ncol = 1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_viridis_d()


# on rescale les années par commodité
chenilles$annee_resc <- chenilles$annee-min(chenilles$annee)+1

#### avec un GLM ####

# on commence juste sur une région
chenilles.nord <- subset(chenilles,region=="Centre Nord semi-oceanique")
mod.nord <- glm(cbind(nbpinsattac,nbpins-nbpinsattac)~annee_resc,family=binomial,data=chenilles.nord)

par(mfrow = c(2,2))
plot(mod.nord)

summary(mod.nord)
p.nord <- ggpredict(mod.nord,terms = "annee_resc")%>%
  plot(residuals=T,log.y=T)
p.nord

# GLM sur toutes les régions
mod0 <- glm(cbind(nbpinsattac,nbpins-nbpinsattac)~annee_resc,family=binomial,data=chenilles)

# résidus
par(mfrow=c(2,2))
plot(mod0)

# interprétation
summary(mod0)
library(questionr)
odds.ratio(mod0)

# cartographie des résidus
chenilles$res.mod0 <- residuals(mod0)
p0 <- plot(ggpredict(mod0, terms = c("annee_resc")), residuals = T)

p1a <- ggplot(regions) +
  geom_sf() +
  geom_point(data = chenilles, 
             mapping = aes(x = longitude, y = latitude, colour = res.mod0),
             size = 2) +
  scale_color_viridis(discrete = F) + 
  labs(color = "Résidus du \n GLM binomial") + 
  theme_classic()

p1b <- ggplot(chenilles) +
  aes(x = placette, y = res.mod0) +
  geom_boxplot() +
  labs(x = "Placettes", y = "Résidus du GLM binomial") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))

p1a + p1b

## modèle mixte

# GLMM avec effet placettes dans GRECO 
mod2 <- glmer(cbind(nbpinsattac,nbpins-nbpinsattac)~annee_resc+(1|region/placette),family=binomial,data=chenilles)
summary(mod2)

# résidus
plot(mod2, type=c("p","smooth"), col.line=1)
plot(mod2,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1,xlab="Valeurs prédites",ylab="sqrt(abs(résidus))")

plot(mod2,xlab="Valeurs prédites",ylab="Résidus de Pearson du GLMM")
qqnorm(scale(residuals(mod2)))
abline(0,1)

# interprétation
summary(mod2)

library(broom.mixed)
tidy(mod2,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="Wald")
tidy(mod2,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

# surdispersion
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod2)
testDispersion(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = mod2, re.form = NULL)
testDispersion(simulationOutput,alternative="greater")
testDispersion(simulationOutput, alternative = "less", plot = FALSE) # seulement sous-dispersion
testDispersion(simulationOutput, alternative = "greater", plot = FALSE) # seulement surdispersion

# graphiques de résultats : effets fixes
p3 <- plot(ggpredict(mod2,terms="annee_resc"),residuals=T)+
  labs(x="Années",y="Taux d'infestation",title="")
p3

p4 <- plot(ggpredict(mod2,type="re"),residuals=T)
p4

p5 <- plot(ggpredict(mod2,terms=c("annee_resc","region"),type="re"),residuals=T)
p5

p6 <- plot(ggpredict(mod2,terms=c("annee_resc","placette"),type="re"),residuals=T)
p6

# effets aléatoires
library(sjPlot)
plot_model(mod2, type = "re", show.values = TRUE)

# R²
library(MuMIn)
r.squaredGLMM(mod2)

# figure sur effets aléatoires
rpoints <- st_sample(regions, 700) %>% # random points, as a list ...
  st_sf() %>%  # ... to data frame ...
  st_transform(4326)  # ... and a metric CRS

ggplot(regions) +
  geom_sf() + 
  geom_sf(data = subset(regions, jointure_g %in% c("A", "B", "F", "G", "J")), 
          aes(fill = jointure_g)) +
  geom_point(data = chenilles, 
             mapping = aes(x = longitude, y = latitude), 
             colour = "red", 
             pch = 17, 
             size = 2) +
  geom_sf(data = rpoints, pch = 3, col = 'white', alpha = 0.67) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) + 
  labs(color = "Régions") + 
  theme_classic() +
  theme(legend.position = "none")


#----------------------------------------#
#### Communautés d'oiseaux forestiers ####
# effets aléatoires sur l'intercept #
#----------------------------------------#

# charger les données
dperche1=read.table("richesse_avifaune_perche.txt",header=T,sep="\t")
dperche1$massif=factor(dperche1$massif)

## Premier modèle : GLM

# Richesse spécifique en fonction des variables visées, tous massifs confondus
mod.perche1=glm(rs~Hdom+ST+essence,data=dperche1,family=poisson)

# contrôle résidus
par(mfrow=c(2,2))
plot(mod.perche1) # quelques points extrêmes mais qui ne sont pas dérangeants (on reste bien en dessous des niveaux de leverage qui influent sur les coefficients)

# l'effet massif est bien visible
boxplot(residuals(mod.perche1,type="pearson")~dperche1$massif, main="résidus ~ massif")
summary(aov(residuals(mod.perche1,type="pearson")~dperche1$massif))

# Richesse spécifique en fonction des variables visées, avec effet fixe massif
mod.perche2=glm(rs~Hdom+ST+essence+massif,data=dperche1,family=poisson)
summary(mod.perche2) # suffisamment de points par massif pour un effet fixe, tous les effets sont bien estimés malgré les hétérogénéités d'échantillonnage

## on passe en modèle mixte pour représenter la hiérarchie de l'échantillonnage et concentrer cette hiérarchie en un terme unique de variance

# ce modèle génère un warning
mod.perche3=glmer(rs~Hdom+ST+essence+(1|massif),data=dperche1,family=poisson)

# étude des variables
summary(dperche1$Hdom)
summary(dperche1$ST)

# le warning est généré par le fait qu'on a des variables sur de grandes échelles de valeurs. Une solution est de centrer - réduire (recommandé)
dperche1$sHdom=as.vector(scale(dperche1$Hdom))
dperche1$sST=as.vector(scale(dperche1$ST))

# variables centrées-réduites
summary(dperche1$sHdom)
summary(dperche1$sST)

# on refait le modèle avec les variables centrées réduites
mod.perche3=glmer(rs~sHdom+sST+essence+(1|massif),data=dperche1,family=poisson)


# diagnostic sur résidus
p1=plot(mod.perche3) 

res=data.frame(res=residuals(mod.perche3))
p2=ggplot(res, aes(sample = res)) +
	stat_qq() +
	stat_qq_line() 


p3=ggplot(data.frame(lev=hatvalues(mod.perche3),pearson=residuals(mod.perche3,type="pearson")),
					aes(x=lev,y=pearson)) +
	geom_point() +
	theme_bw() 

plot_grid(p1,p2,p3)

# effet blocs bien traité?
boxplot(residuals(mod.perche3)~dperche1$massif)

# R²
r.squaredGLMM(mod.perche3)

# sortie numérique
summary(mod.perche3)

# sorties graphiques des effets avec sjPlot : avec type="eff", prédicteurs constants à 0
p.mix1=plot_model(mod.perche3,type="eff",terms="sST")
p.mix2=plot_model(mod.perche3,type="eff",terms="sHdom")
p.mix3=plot_model(mod.perche3,type="eff",terms="essence")
p.mix4=plot_model(mod.perche3,type="re")
cowplot::plot_grid(p.mix1,p.mix2,p.mix3,p.mix4)

# avec ggeffects, fonction ggemmeans : prédicteurs maintenus à leur valeur moyenne
p.eff1=plot(ggemmeans(mod.perche3,terms="sST"),residuals=T)
p.eff2=plot(ggemmeans(mod.perche3,terms="sHdom"),residuals=T)
p.eff3=plot(ggemmeans(mod.perche3,terms="essence"),residuals=T)
cowplot::plot_grid(p.eff1,p.eff2,p.eff3)

# modèle initial (GLM) avec variables centrées-réduites
mod.perche1b=glm(rs~sHdom+sST+essence,data=dperche1,family=poisson)

# comparaison des intervalles de confiance (sauf intercept pour faciliter la représentation graphique)
conf.massif=confint(mod.perche3,method="Wald")[-c(1,2),] # on utilise les IC de Wald parce qu'il y a un gros df mais il vaudrait mieux calculer les IC profile (! long)
conf.nomassif=confint(mod.perche1b,method="Wald")[-1,]

# graphiquement : les différences d'IC ne sont pas très fortes parce que l'effet massif est faible, mais on voit quand même de petits décalages
# ces décalages seraient bien plus forts si le ddl était plus faible et avec des gros déséquilibres d'échantillonnage entre massifs
plot(x=0,y=0,type="n",xlim=c(1,13),xaxt="n",xlab="pentes",ylab="intervalle de confiance (Wald)",ylim=c(-0.4,0.4))
axis(side=1,at=seq(1.5,12.5,2),labels=c("Hdom","ST","DOU","EPC","PS","SP"))
abline(h=0,col="steelblue",lty="dotted")
segments(x0=seq(1,12,2),x1=seq(1,12,2),y0=conf.massif[,1],y1=conf.massif[,2],lwd=2)
segments(x0=seq(2,13,2),x1=seq(2,13,2),y0=conf.nomassif[,1],y1=conf.nomassif[,2],lwd=2,col="darkred")
legend("topright",bty="n",lty="solid",col=c("black","darkred"),legend=c("effet aléatoire massif","pas d'effet massif"))

#--------------------------------------#
### Modèle mixte à pentes aléatoires ###
#--------------------------------------#

## simulation
x=rnorm(30,0,1) # variable explicative
a=c(0,1,2,5,10) # intercept par blocs
b=c(-0.3,3,0,1,-0.2) # pente par blocs
eps=rnorm(30,0,0.5) # erreur intra-bloc commune

y=matrix(NA,nrow=30,ncol=5)
colnames(y)=paste("A",1:5,sep="")
for(i in 1:5){
	y[,i]=a[i]+b[i]*x+eps
	}
y2=melt(y)
y2$x=rep(x,times=5)

m1=lm(value~x,data=y2) # modèle linéaire sans effet bloc
m2=lmer(value~x+(1|Var2),data=y2)
m3=lmer(value~x+(x|Var2),data=y2) # forte sensibilité dès que les pentes sont trop hétérogènes

confint(m1)
confint(m2)
confint(m3) # comparez les intervalles de confiance

p0=ggplot(y2)+aes(x=x,y=value,col=Var2)+	geom_point()+labs(title="données brutes")
p1=plot(ggpredict(m1,terms="x"),add.data=T)+labs(title="sans effet bloc")
p2=plot(ggpredict(m2,terms="x"),add.data=T)+labs(title="effet aléatoire sur l'intercept")
p3=plot(ggpredict(m3,terms="x"),add.data=T)+labs(title="effet aléatoire intercept et pente")
# la pente de x et son intercept ne changent pas, mais les intervalles de confiance sont fortement affectés par la structure en blocs

cowplot::plot_grid(p0,p1,p2,p3)

## GLM classique : terme d'interaction
mod.perche4=glm(rs~massif*(sHdom+sST+essence),data=dperche1,family=poisson)

par(mfrow=c(2,2))
plot(mod.perche4)
summary(mod.perche4)

# interactions non estimables
table(dperche1[,c("massif","essence")])

# GLMM à pentes aléatoires
mod.perche5=glmer(rs~sHdom+sST+essence+(sHdom+sST+essence|massif),data=dperche1,family=poisson,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # identique à (1+sHdom+sST+essence|massif). Le changement d'algorithme est rendu nécessaire par la complexité du modèle

# contrôlez vous-mêmes les résidus (voir fonctions plus haut) - NB : quand ils sont corrects pour un modèle simple, il le sont aussi pour un 
# modèle plus complexe, sauf artefact qui suggèrerait une erreur dans la conception du modèle

# intervalles de confiance
confint(mod.perche5,method="Wald") # à nouveau : penser à passer en "profile" pour avoir des IC plus corrects, mais !! durée de calcul

# sorties graphiques du modèle à pentes aléatoires
pmixslope1=plot_model(mod.perche5,type="eff",terms="sST",title="Surface terrière, intercept et pentes aléatoires")+ylim(10,20)
pmixslope2=plot_model(mod.perche5,type="eff",terms="sHdom",title="Hauteur dominante, intercept et pentes aléatoires")+ylim(0,20)
pmixslope3=plot_model(mod.perche5,type="eff",terms="essence",title="Essence, intercept et pentes aléatoires")

# sorties graphiques du modèle à intercept aléatoire
pmix1=plot_model(mod.perche3,type="eff",terms="sST",title="Surface terrière, intercept aléatoire")+ylim(10,20)
pmix2=plot_model(mod.perche3,type="eff",terms="sHdom",title="Hauteur dominante, intercept aléatoire")+ylim(0,20)
pmix3=plot_model(mod.perche3,type="eff",terms="essence",title="Essence, intercept aléatoire") # VID n'est présent que sur Perche-Trappe

cowplot::plot_grid(pmix1,pmix2,pmix3,pmixslope1,pmixslope2,pmixslope3)

#--------------------------------#
### Phénologie du débourrement ###
#--------------------------------#

# phénologie du débourrement des chênes et des pins en forêt d'Orléans
d_pheno=read.table("phenologie_debourrement.txt",header=T,sep="\t")

# on explore les données
head(d_pheno)
dim(d_pheno)
str(d_pheno)
d_pheno$annee=factor(d_pheno$annee)

# on représente les données brutes
p1=ggplot(d_pheno)+
	aes(x=essence,y=date_debourrement)+
		geom_boxplot()

p2=ggplot(d_pheno)+
	aes(x=annee,y=date_debourrement)+
	geom_boxplot()

p3=ggplot(d_pheno)+
	aes(x=diametre,y=date_debourrement)+
	geom_point()

cowplot::plot_grid(p1,p2,p3)

# interaction diamètre*essence
ggplot(d_pheno)+
	aes(x=diametre,y=date_debourrement)+
		geom_point()+
			facet_wrap(~essence)


hist(d_pheno$date_debourrement) # on ne pourra pas appliquer de transformation dans ce cas (aucune interprétabilité bio + pas de transfo permettant
# de normaliser une bimodalité), mais ce n'est pas grave car l'origine de la bimodalité est prise en compte dans le modèle (essence)

## modèle sans effet hiérarchique
mod.db0=lm(date_debourrement~essence+diametre,data=d_pheno)

#  plan d'échantillonnage
table(d_pheno[,c("code_placette","code_parcelle")])

## modèle mixte à effet placette
mod.db1=lmer(date_debourrement~essence+diametre+(1|code_placette),data=d_pheno)

# on tient compte du niveau d'emboîtement "parcelles"
mod.db2=lmer(date_debourrement~essence+diametre+(1|code_parcelle/code_placette),data=d_pheno) 
summary(mod.db2)

# comparaison graphique avec / sans emboîtement
mod.emb0=plot_model(mod.db0,type="eff",terms="essence",title="sans structure hiérarchique") 
mod.emb1=plot_model(mod.db1,type="eff",terms="essence",title="effet aléatoire placette") 
mod.emb2=plot_model(mod.db2,type="eff",terms="essence",title="effet aléatoire placette/ parcelle") 
cowplot::plot_grid(mod.emb0,mod.emb1,mod.emb2)

# effet de l'année
dev.off()
boxplot(residuals(mod.db2)~ d_pheno$annee,main="résidus ~ année") 

## modèle à effets aléatoires croisés
mod.db3=lmer(date_debourrement~essence+diametre+(1|code_parcelle/code_placette)+(1|annee),data=d_pheno) 
summary(mod.db3)

# attention : dans ce modèle, il n'y a que 2 niveaux pour l'effet année --> insuffisant pour bien estimer le terme de variance
mod.db4=lmer(date_debourrement~essence+diametre+annee+(1|code_parcelle/code_placette),data=d_pheno) 
summary(mod.db4)

mod.db5=lmer(date_debourrement~essence+diametre+(1|annee)+(1|code_parcelle/code_placette),data=d_pheno) 

# effet essence
m.spat=plot(ggemmeans(mod.db3,terms="essence"),residuals=T)+labs(title="sans effet année")
m.spat.an=plot(ggemmeans(mod.db4,terms="essence"),residuals=T)+labs(title="avec effet année fixe")
m.spat.an2=plot(ggemmeans(mod.db5,terms="essence"),residuals=T)+labs(title="avec effet aléatoire année")
cowplot::plot_grid(m.spat,m.spat.an,m.spat.an2)

# effet diamètre
m.spat=plot(ggemmeans(mod.db3,terms="diametre"),residuals=T)+labs(title="sans effet année")
m.spat.an=plot(ggemmeans(mod.db4,terms="diametre"),residuals=T)+labs(title="avec effet année fixe")
m.spat.an2=plot(ggemmeans(mod.db5,terms="diametre"),residuals=T)+labs(title="avec effet aléatoire année")
cowplot::plot_grid(m.spat,m.spat.an,m.spat.an2)

# effets aléatoires
dotplot(ranef(mod.emb3, condVar = TRUE,whichel="code_parcelle")) # effet parcelle
dotplot(ranef(mod.emb3, condVar = TRUE,whichel="code_placette:code_parcelle")) # effet placette
dotplot(ranef(mod.emb3, condVar = TRUE,whichel="annee")) # effet année

## effet fixe année: solution plus réaliste du point de vue de l'échantillonnage
mod.emb4=lmer(date_debourrement~essence+diametre+annee+(1|code_parcelle/code_placette),data=d_pheno) 

# sortie graphique
plot_model(mod.emb4,type="eff",terms=c("essence","annee"))

# pour plus d'info, voir le livre de Douglas Bates, dispo en ligne
# http://lme4.r-forge.r-project.org/book/Ch1.pdf
# pour l'évaluation des effets aléatoires, voir C1.6




