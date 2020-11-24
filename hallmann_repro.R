## Hallmann 2017

library(R2jags)
library(RCurl) # used for loading data from github

# ---- Load & Standardize ----
# Read the data
url_1 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s1_data.csv")
url_2 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s2_data.csv")

data <- read.csv(text = url_1, header = TRUE, sep = ",")
model.frame <- read.csv(text = url_2, header = TRUE, sep = ",")

# summary(data)
# summary(model.frame)

# ---- Null Model ----

## Indices contain the corresponding rows in model.frame of the  starting (1)
## and end (2) position of the potID groups; tau_1/2 in the paper
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

## NOT RUN
# sink("NullModel.jag")
# cat("
# model{
# ## likelihood function for the latent expected daily biomass
# for (i in 1:n){
# m[i] ~ dnorm(sum(y[index1[i]:index2[i]]), tau[i])
# tau[i] <- 1 / Var[i] ## 1/variance: using precision parameter not Variance  
# Var[i] <- sum(vr[index1[i]:index2[i]])
# }
# 
# ## likelihood function for muHat (y), it's dependent function (z) and variance (vr)
# for (i in 1:ndaily){
# y[i]  <- exp(z[i]) ## z is actually called y in the paper
# z[i]  <- int + c[1] * daynr[i] + c[2] * daynr2[i] + b[loctype[i]] + eps[plot[i]]
# ## this is called y in the paper
# vr[i] <- exp(2 * z[i] + lvar) * (exp(lvar) - 1) ## variance
# }
# 
# ## priors
# int ~ dnorm(0,.01) ## global intercept (see eq. 4 in paper)
# b[1] <- 0  
# for (i in 2:3) {b[i] ~ dnorm(0, .01)}
# for (i in 1:2) {c[i] ~ dnorm(0, .01)}
# sdhat ~ dunif(0, 5)
# lvar <- pow(sdhat, 2) ## variance of daily log-biomass
# for(i in 1:nrandom) {
# eps[i] ~ dnorm(0, tau.re) ## location-specific random effect (see eq. 4 in paper)
# }
# tau.re <- pow(sd.re, -2) ## convert sd to precision parameter of normal dist
# sd.re ~ dunif(0, 1)
# }
# ")
# sink(NULL)

# Run the model
parametersNull <- c("int", "b", "c", "eps", "sdhat", "sd.re")
jagsmod0 <- jags(jagsdataNull, inits = NULL, parametersNull,
                "NullModel.jag", n.chains = 3, n.iter = 120, n.burnin = 20,
                n.thin = 10)
jagsmod0

# ---- Re-Analysis ----
library(dplyr)

# Filter Biomass data (only first Dataset of every plot)
py <- select(data, c(plot, year)) %>%
  group_by(plot) %>%
  distinct(plot, year) %>%
  arrange(plot, year) # overview plot by year

new_plots <- group_by(py, plot) %>%
 slice(1) # get first sampling year of each plot

new_data <- semi_join(data, new_plots, by = c("plot", "year"))
new_model.frame <- semi_join(model.frame, new_plots, by = c("plot", "year"))

# New Null model (with meaningful/unambiguous variable names)
newjagsdataNull <- list(
  m_bio = new_data$biomass,
  tau1 = with(new_model.frame, tapply(1:nrow(new_model.frame), potID, min)), 
  tau2 = with(new_model.frame, tapply(1:nrow(new_model.frame), potID, max)),
  plot = as.numeric(new_model.frame$plot),
  loctype = as.numeric(new_data$location.type[match(new_model.frame$potID, new_data$potID)]),
  daynr = as.numeric((new_model.frame$daynr - mean(new_data$mean.daynr)) / sd(new_data$mean.daynr)),
  daynr2 = as.numeric((new_model.frame$daynr - mean(new_data$mean.daynr)) / sd(new_data$mean.daynr))^2,
  ndaily = nrow(new_model.frame),
  n = nrow(new_data),
  nrandom = max(as.numeric(new_model.frame$plot))
)

# Jags model
## NEEDS TO RUN ONCE TO CREATE .JAGS FILE
## SET WD TO APROPRIATE DIRECTORY
{
# sink("newNullModel.jag")
# cat("model{
# ## Likelihood function for the latent expected daily biomass
# for (i in 1:n){
# m_bio[i] ~ dnorm(sum(z[tau1[i]:tau2[i]]), sig_sq[i]) # measured daily biomass
# sig_sq[i] <- 1 / Var[i] # Precision
# Var[i] <- sum(vr[tau1[i]:tau2[i]]) # Variance of measured biomass
# }
# 
# ## Likelihood function for muHat, it's dependent function and variance
# for (i in 1:ndaily){
# z[i]  <- exp(y[i]) # Latent expected daily biomass
# y[i]  <- g_intcp + c[1] * daynr[i] + c[2] * daynr2[i] + b[loctype[i]] + eps[plot[i]]
# vr[i] <- exp(2 * y[i] + lvar) * (exp(lvar) - 1) # Variance (see paper eq. 6)
# }
# 
# ## Priors
# g_intcp ~ dnorm(0,.01) # global intercept (see eq. 4 in paper [c])
# b[1] <- 0 
# for (i in 2:3) {b[i] ~ dnorm(0, .01)} # Habitat cluster 2 & 3 covariates
# for (i in 1:2) {c[i] ~ dnorm(0, .01)} # coef for daynr & naynr^2 covariates
# sdhat ~ dunif(0, 5) # Estimated sd
# lvar <- pow(sdhat, 2) # Residual variance of daily log-biomass
# for(i in 1:nrandom) {
# eps[i] ~ dnorm(0, sig_sq.re) # Location-specific random effect (see eq. 4 in paper [u_s])
# }
# sig_sq.re <- pow(sd.re, -2) # Residual [sigma^2_j]; pow():Convert sd to percision parameter
# sd.re ~ dunif(0, 1) # Residual sd
# }
# ")
# sink(NULL)
}

# Run the model
parametersNull <- c("int", "b", "c", "eps", "sdhat", "sd.re")
newjagsmod0 <- jags(newjagsdataNull, inits = NULL, parametersNull,
                 "newNullModel.jag", n.chains = 3, n.iter = 120, n.burnin = 20,
                 n.thin = 10)
newjagsmod0

# ---- Visualization ----
library(ggplot2)

# Costum palette to match colors from paper
my_pal <- c("#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9",
            "#74add1", "#4575b4")

# Recreate Fig. 2a
tot_bm_data <- mutate(new_data, bm_p_day = biomass / (to.daynr - from.daynr)) %>%
  select(year, bm_p_day)

ggplot(tot_bm_data, aes(x = factor(year, levels = seq(1989, 2014)),
                        y = bm_p_day,
                        fill = year)) +
  geom_boxplot() +
  #geom_jitter(width = .1, alpha = .2) + # optional: datapoints
  scale_fill_gradient2(low = "#4575b4",
                        mid = "#fee090",
                        high = "#f46d43",
                        na.value = "grey50",
                        midpoint = 2005,
                        guide = "none") +
  scale_x_discrete(breaks = seq(1990, 2015, by = 5), drop = F) +
  theme_classic() +
  labs(y = "Biomass [g/d]", x = "Year")

# Recreate Fig. 2b
# Season from 1st of April to 30th of October; i.e. Day 91 to 303
s_bm_data <- filter(data, from.daynr >= 91 & to.daynr <= 303) %>%
  mutate(bm_p_day = biomass / (to.daynr - from.daynr)) %>% 
  select(year, bm_p_day, mean.daynr)

ggplot(s_bm_data) +
  geom_point(aes(x = mean.daynr, y = bm_p_day, color = year)) +
  scale_color_gradient2(low = "#4575b4",
                        mid = "#fee090",
                        high = "#f46d43",
                        na.value = "grey50",
                        midpoint = 2005,
                        guide = "none") +
  theme_classic() +
  labs(y = "Biomass [g/d]", x = "Day of year")
  
