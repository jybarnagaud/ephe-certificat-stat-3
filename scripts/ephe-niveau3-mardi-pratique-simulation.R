# ======================================
# Simulation + JAGS + diagnostics
# ======================================

library(rjags)
library(coda)

set.seed(123)

# -----------------------------
# 1. SIMULATION
# -----------------------------
npoints <- 643
nrepl   <- 3
nyears  <- 2

PROPNATFOR <- rnorm(npoints)
altitude   <- rnorm(npoints)
# 350 sites en année 1, 293 en année 2
year <- c(rep(1, 350), rep(2, 293))
# optionnel: mélanger l'ordre des sites
year <- sample(year)
julian_date <- scale(rnorm(npoints, 150, 20))[,1]

# vrais paramètres
true <- list(
  det.alpha = qlogis(0.65),  # garantit p moyen ~0.65
  beta = 0.8,
  delta = -0.6,
  muyear = 0.2,
  sigalpha = 0.7,
  sigyear = 0.5,
  pmean = 0.65
)

# on fixe beta.det = 0 pour contrôler facilement pmean
beta.det <- rep(0, nyears)
mualpha  <- rnorm(nyears, true$muyear, true$sigyear)

alpha <- rnorm(npoints, mualpha[year], true$sigalpha)

lambda <- exp(alpha + true$beta*PROPNATFOR + true$delta*altitude)
N <- rpois(npoints, lambda)

p <- plogis(true$det.alpha + beta.det[year]*julian_date)

obs <- matrix(NA, npoints, nrepl)
for(i in 1:npoints){
  for(j in 1:nrepl){
    obs[i,j] <- rbinom(1, N[i], p[i])
  }
}

# -----------------------------
# 2. MODELE JAGS
# -----------------------------
model_string <- "
model{

det.alpha~dnorm(0,0.0001)
beta~dnorm(0,0.0001)
delta~dnorm(0,0.0001)
muyear~dnorm(0,0.0001)

taualpha<-1/(sigalpha*sigalpha)
sigalpha~dunif(0,100)

tauyear<-1/(sigyear*sigyear)
sigyear~dunif(0,100)

for (k in 1:nyears){
  beta.det[k]~dnorm(0,0.0001)
  alpha[k]~dnorm(0,tauyear)
}

for(i in 1:npoints){
  N[i] ~ dpois(esp[i])
  log(esp[i])<-alpha[year[i]]+beta*PROPNATFOR[i]+delta*altitude[i]

  
  for(j in 1:nrepl){
    obs[i,j] ~ dbin(p[i],N[i])
  }

  logit(p[i])<-det.alpha+beta.det[year[i]]*julian_date[i]

  for(j in 1:nrepl){
    eval[i,j]<-N[i]*p[i]
    E[i,j]<-pow((obs[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)

    obs.new[i,j]~dbin(p[i],N[i])
    E.new[i,j]<-pow((obs.new[i,j]-eval[i,j]),2)/(eval[i,j]+0.5)
  }
}

fit<-sum(E[,])
fit.new<-sum(E.new[,])

pmean<-mean(p[])
}
"

# -----------------------------
# 3. DONNEES JAGS
# -----------------------------
data_jags <- list(
  npoints=npoints,
  nrepl=nrepl,
  nyears=nyears,
  obs=obs,
  PROPNATFOR=PROPNATFOR,
  altitude=altitude,
  year=year,
  julian_date=julian_date
)

# initial values
inits <- function(){
  list(N = apply(obs,1,max)+1)
}

params <- c("det.alpha","beta","delta","muyear",
            "sigalpha","sigyear","beta.det",
            "fit","fit.new","pmean")

# -----------------------------
# 4. FIT JAGS
# -----------------------------
mod <- jags.model(textConnection(model_string),
                  data=data_jags,
                  inits=inits,
                  n.chains=3,
                  n.adapt=1000)

update(mod, 2000)

samp <- coda.samples(mod,
                     variable.names=params,
                     n.iter=5000)

# -----------------------------
# 5. RECUPERATION PARAMETRES
# -----------------------------
summary_stats <- summary(samp)$statistics
print(summary_stats)

# comparaison vrais vs estimés
cat("\nComparaison:\n")
cat("beta vrai:", true$beta, "\n")
cat("delta vrai:", true$delta, "\n")
cat("det.alpha vrai:", true$det.alpha, "\n")

# -----------------------------
# 6. PPC
# -----------------------------
fit <- as.matrix(samp)[,"fit"]
fit.new <- as.matrix(samp)[,"fit.new"]

p_value <- mean(fit.new > fit)
cat("Bayesian p-value:", p_value, "\n")

# -----------------------------
# 7. GRAPHIQUES
# -----------------------------
par(mfrow=c(2,2))

# traceplots
traceplot(samp[,c("beta","delta","det.alpha")])

# densités
plot(density(as.matrix(samp)[,"beta"]), main="beta")
abline(v=true$beta, col=2, lwd=2)

plot(density(as.matrix(samp)[,"delta"]), main="delta")
abline(v=true$delta, col=2, lwd=2)

# densité de pmean
plot(density(as.matrix(samp)[,"pmean"]), main="pmean")
abline(v=true$pmean, col=2, lwd=2)

# PPC
plot(fit, fit.new,
     xlab="fit", ylab="fit.new",
     main="Posterior predictive check")
abline(0,1,col=2,lwd=2)
