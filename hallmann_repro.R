## Hallmann 2017

# required libraries
library(R2jags)
library(RCurl) # used for loading data from github

# ---- Load & Standardize ----

# Function to standardize variables
standardize <- function(x, mu = mean(x), sig = sd(x)) (x-mu)/sig

# read the data
url_1 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s1_data.csv")
url_2 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s2_data.csv")

data <- read.csv(text = url_1, header = TRUE, sep = ",")
model.frame <- read.csv(text = url_2, header = TRUE, sep = ",")

# summary(data)
# summary(model.frame)

# ---- Null Model ----

## Indices contain the corresponding rows in model.frame of the  starting (1)
## and end (2) position of the potID groups
## plot is used to get plot info for the location-specific random effect (eps)
## loctype is used for the habitat cluster covariate

jagsdataNull <- list(
  m = data$biomass,
  index1 = with(model.frame, tapply(1:nrow(model.frame), potID, min))[data$potID], 
  index2 = with(model.frame, tapply(1:nrow(model.frame), potID, max))[data$potID],
  plot = as.numeric(model.frame$plot),
  loctype = as.numeric(data$location.type[match(model.frame$potID, data$potID)]),
  daynr = as.numeric((model.frame$daynr - mean(data$mean.daynr)) / sd(data$mean.daynr)),
  daynr2 = as.numeric((model.frame$daynr - mean(data$mean.daynr)) / sd(data$mean.daynr))^2,
  ndaily = nrow(model.frame),
  n = nrow(data),
  nrandom = max(as.numeric(model.frame$plot))
)

sink("NullModel.jag")
cat("
model{
## likelihood function for the latent expected daily biomass
for (i in 1:n){
m[i] ~ dnorm(sum(y[index1[i]:index2[i]]), tau[i])
tau[i] <- 1 / Var[i]
Var[i] <- sum(vr[index1[i]:index2[i]])
}

## likelihood function for muHat (y), it's dependent function (z) and variance (vr)
for (i in 1:ndaily){
y[i]  <- exp(z[i]) ## z is actually called y in the paper
z[i]  <- int + c[1] * daynr[i] + c[2] * daynr2[i] + b[loctype[i]] + eps[plot[i]]
## this is called y in the paper
vr[i] <- exp(2 * z[i] + lvar) * (exp(lvar) - 1) ## variance
}

## priors
int ~ dnorm(0,.01) ## global intercept (see eq. 4 in paper)
b[1] <- 0  
for (i in 2:3) {b[i] ~ dnorm(0, .01)}
for (i in 1:2) {c[i] ~ dnorm(0, .01)}
sdhat ~ dunif(0, 5)
lvar <- pow(sdhat, 2) ## variance of daily log-biomass
for(i in 1:nrandom) {
eps[i] ~ dnorm(0, tau.re) ## location-specific random effect (see eq. 4 in paper)
}
tau.re <- pow(sd.re, -2) ## convert sd to percision parameter of normal dist
sd.re ~ dunif(0, 1)
}
")
sink(NULL)

# Run the model
parametersNull <- c("int", "b", "c", "eps", "sdhat", "sd.re")
jagsmod0 <- jags(jagsdataNull, inits = NULL, parametersNull,
                "NullModel.jag", n.chains = 3, n.iter = 120, n.burnin = 20,
                n.thin = 10)
jagsmod0




