## Hallmann 2017

setwd("~/Documents/Uni/Umwi/M.sc./Repro")
load("full_basic")


library(R2jags)
library(RCurl) # used for loading data from github
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# ---- Load ----
# Read the data
url_1 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s1_data.csv")
url_2 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s2_data.csv")

data <- read.csv(text = url_1, header = TRUE, sep = ",", stringsAsFactors = T)
model.frame <- read.csv(text = url_2, header = TRUE, sep = ",", stringsAsFactors = T)

# summary(data)
# summary(model.frame)

# ---- Original basic model ----
jagsdataBasic <- list(
  m_bio = data$biomass,
  tau1 = with(model.frame, tapply(1:nrow(model.frame), potID, min)),
  tau2 = with(model.frame, tapply(1:nrow(model.frame), potID, max)),
  plot = as.numeric(model.frame$plot),
  loctype = as.numeric(data$location.type[match(data$potID,data$potID)]),
  daynr = as.numeric((model.frame$daynr-mean(data$mean.daynr)) / sd(data$mean.daynr)),
  daynr2 = as.numeric((model.frame$daynr-mean(data$mean.daynr)) / sd(data$mean.daynr))^2,
  year = model.frame$year - 1988,
  ndaily = nrow(model.frame),
  n = nrow(data),
  nrandom = max(as.numeric(model.frame$plot))
)

parametersBasic <- c("g_intcp", "log.lambda", "b", "c", "eps", "sdhat", "sd.re")
oldjagsmodBasic <- jags(jagsdataBasic, inits = NULL, parametersBasic, 
                        "BasicModel.jag", n.iter = 12000, n.burnin = 2000, 
                        n.chains = 3, n.thin = 10)

# Get mean parameter values from posterior
parms_mean <- oldjagsmodBasic$BUGSoutput$mean 

# 2.5 % credible interval
parms_mean_lo <- list("b" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[1:3, 3]),
                          "c" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[4:7, 3]),
                          "deviance" = oldjagsmodBasic$BUGSoutput$summary[8, 3],
                          "eps" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[9:71, 3]),
                          "g_intcp" = oldjagsmodBasic$BUGSoutput$summary[72, 3],
                          "log.lambda" = oldjagsmodBasic$BUGSoutput$summary[73, 3],
                          "sd.re" = oldjagsmodBasic$BUGSoutput$summary[74, 3],
                          "sdhat" = oldjagsmodBasic$BUGSoutput$summary[75, 3])

# 97.5 % credible interval
parms_mean_hi <- list("b" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[1:3, 7]),
                          "c" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[4:7, 7]),
                          "deviance" = oldjagsmodBasic$BUGSoutput$summary[8, 7],
                          "eps" = as.numeric(oldjagsmodBasic$BUGSoutput$summary[9:71, 7]),
                          "g_intcp" = oldjagsmodBasic$BUGSoutput$summary[72, 7],
                          "log.lambda" = oldjagsmodBasic$BUGSoutput$summary[73, 7],
                          "sd.re" = oldjagsmodBasic$BUGSoutput$summary[74, 7],
                          "sdhat" = oldjagsmodBasic$BUGSoutput$summary[75, 7])

# ---- Re-Analysis ----
 # ---- Data wrangling ----
# Filter Biomass data (only first data set of every plot)
py <- select(data, c(plot, year)) %>%
  group_by(plot) %>%
  distinct(plot, year) %>%
  arrange(plot, year) # overview plot by year

new_plots <- group_by(py, plot) %>%
 slice(1) # get first sampling year of each plot

# Create new dfs
new_data <- semi_join(data, new_plots, by = c("plot", "year"))
new_model.frame <- semi_join(model.frame, new_plots, by = c("plot", "year"))

 # ---- Basic Model ----
newjagsdataBasic <- list(
  m_bio = new_data$biomass,
  tau1 = with(new_model.frame, tapply(1:nrow(new_model.frame), potID, min)),
  tau2 = with(new_model.frame, tapply(1:nrow(new_model.frame), potID, max)),
  plot = as.numeric(new_model.frame$plot),
  loctype = as.numeric(new_data$location.type[match(new_model.frame$potID,new_data$potID)]),
  daynr = as.numeric((new_model.frame$daynr-mean(new_data$mean.daynr)) / sd(new_data$mean.daynr)),
  daynr2 = as.numeric((new_model.frame$daynr-mean(new_data$mean.daynr)) / sd(new_data$mean.daynr))^2,
  year = new_model.frame$year - 1988,
  ndaily = nrow(new_model.frame),
  n = nrow(new_data),
  nrandom = max(as.numeric(new_model.frame$plot))
)

# Jags model
## NEEDS TO RUN ONCE TO CREATE .JAGS FILE
## SET WD TO APROPRIATE DIRECTORY
{
# sink("BasicModel.jag")
# cat("model{
# ## Likelihood function for the latent expected daily biomass
# for (i in 1:n) {
# m_bio[i] ~ dnorm(sum(z[tau1[i]:tau2[i]]), sig_sq[i])
# sig_sq[i] <- 1/Var[i]
# Var[i] <- sum(vr[tau1[i]:tau2[i]])
# }
# 
# ## Likelihood function for muHat, it's dependent function and variance
# for (i in 1:ndaily) {
# z[i] <- exp(y[i])
# y[i] <- g_intcp + log.lambda * year[i] + c[1] * daynr[i] + c[2] * daynr2[i] +
#   c[3] * daynr[i] * year[i] + c[4] * daynr2[i] * year[i] + b[loctype[i]] +
#   eps[plot[i]]
# vr[i] <- exp(2 * y[i] + lvar) * (exp(lvar) - 1)
# }
# 
# ## Priors
# g_intcp ~ dnorm(0, .01)
# log.lambda ~ dnorm(0, .01)
# b[1] <- 0
# for( i in 2:3) {b[i] ~ dnorm(0, .01)}
# for( i in 1:4) {c[i] ~ dnorm(0, .01)}
# sdhat ~ dunif(0, 5)
# lvar <- pow(sdhat, 2)
# for (i in 1:nrandom) {
# eps[i] ~ dnorm(0, tau.re)
# }
# tau.re <- pow(sd.re, -2)
# sd.re ~ dunif(0, 1)
# }
# ")
# sink(NULL)
}

# Run the model
jagsmodBasic <- jags(newjagsdataBasic, inits = NULL, parametersBasic,
                     "BasicModel.jag", n.iter = 12000, n.burnin = 2000,
                     n.chains = 3, n.thin = 10)

jagsmodBasic

 # ---- Diagnostics ----

 # ---- Calculate predicted values ----
# Get mean parameter values from posterior
new_parms_mean <- jagsmodBasic$BUGSoutput$mean 

# 2.5 % credible interval
new_parms_mean_lo <- list("b" = as.numeric(jagsmodBasic$BUGSoutput$summary[1:3, 3]),
                      "c" = as.numeric(jagsmodBasic$BUGSoutput$summary[4:7, 3]),
                      "deviance" = jagsmodBasic$BUGSoutput$summary[8, 3],
                      "eps" = as.numeric(jagsmodBasic$BUGSoutput$summary[9:71, 3]),
                      "g_intcp" = jagsmodBasic$BUGSoutput$summary[72, 3],
                      "log.lambda" = jagsmodBasic$BUGSoutput$summary[73, 3],
                      "sd.re" = jagsmodBasic$BUGSoutput$summary[74, 3],
                      "sdhat" = jagsmodBasic$BUGSoutput$summary[75, 3])

# 97.5 % credible interval
new_parms_mean_hi <- list("b" = as.numeric(jagsmodBasic$BUGSoutput$summary[1:3, 7]),
                       "c" = as.numeric(jagsmodBasic$BUGSoutput$summary[4:7, 7]),
                       "deviance" = jagsmodBasic$BUGSoutput$summary[8, 7],
                       "eps" = as.numeric(jagsmodBasic$BUGSoutput$summary[9:71, 7]),
                       "g_intcp" = jagsmodBasic$BUGSoutput$summary[72, 7],
                       "log.lambda" = jagsmodBasic$BUGSoutput$summary[73, 7],
                       "sd.re" = jagsmodBasic$BUGSoutput$summary[74, 7],
                       "sdhat" = jagsmodBasic$BUGSoutput$summary[75, 7])
  
  
# Set tau & intervals of days original data
t_1 <- as.vector(jagsdataBasic$tau1[1:nrow(data)])
t_2 <- as.vector(jagsdataBasic$tau2[1:nrow(data)])
ints <- map(1:nrow(data), ~ seq(t_1[.x], t_2[.x]))
  
# Set year intervals original data
t_1_y <- with(data, tapply(1:nrow(data), year, min))
t_2_y <- with(data, tapply(1:nrow(data), year, max))
ints_y <- map(1:length(t_1_y), ~ seq(t_1_y[.x], t_2_y[.x]))

# Set tau & intervals of days
new_t_1 <- as.vector(newjagsdataBasic$tau1[1:nrow(new_data)])
new_t_2 <- as.vector(newjagsdataBasic$tau2[1:nrow(new_data)])
new_ints <- map(1:nrow(new_data), ~ seq(new_t_1[.x], new_t_2[.x]))

# Set year intervals
new_t_1_y <- with(new_data, tapply(1:nrow(new_data), year, min))
new_t_2_y <- with(new_data, tapply(1:nrow(new_data), year, max))
new_ints_y <- map(1:length(new_t_1_y), ~ seq(new_t_1_y[.x], new_t_2_y[.x]))

# Helper functions to calculate predicted biomass
y <- function(x, model, cred) {
  data <- switch(model,
                 new = newjagsdataBasic,
                 old = jagsdataBasic)
  
  if (model == "new") {
    parms_mean <- switch(cred,
                         base = new_parms_mean,
                         low = new_parms_mean_lo,
                         high = new_parms_mean_hi)
  } else {parms_mean <- switch(cred,
                               base = parms_mean,
                               low = parms_mean_lo,
                               high = parms_mean_hi)}
    
  parms_mean$g_intcp + parms_mean$log.lambda * data$year[x] +
    parms_mean$c[1] * data$daynr[x] + parms_mean$c[2] *
    data$daynr2[x] + parms_mean$c[3] * data$daynr[x] *
    data$year[x] + parms_mean$c[4] * data$daynr2[x] *
    data$year[x] + parms_mean$b[data$loctype[x]] +
    parms_mean$eps[data$plot[x]]
}

z <- function(x, ...) exp(y(x, ...))

# Calculate predicted biomass
# Old model
z_val <- unlist(map(1:nrow(data), ~ mean(z(unlist(ints[.x]), "old", "base"))))
z_val_y <- unlist(map(1:length(ints_y), ~ mean(z_val[unlist(ints_y[.x])])))

# 2.5 % credible interval
z_val_lo <- unlist(map(1:nrow(data), ~ mean(z(unlist(ints[.x]), "old", "low"))))
z_val_y_lo <- unlist(map(1:length(ints_y), ~ mean(z_val_lo[unlist(ints_y[.x])])))

# 97.5 % credible interval
z_val_hi <- unlist(map(1:nrow(data), ~ mean(z(unlist(ints[.x]), "old", "high"))))
z_val_y_hi <- unlist(map(1:length(ints_y), ~ mean(z_val_hi[unlist(ints_y[.x])])))

# New model
new_z_val <- unlist(map(1:nrow(new_data), ~ mean(z(unlist(new_ints[.x]), "new", "base"))))
new_z_val_y <- unlist(map(1:length(new_ints_y), ~ mean(new_z_val[unlist(new_ints_y[.x])])))

# 2.5 % credible interval
new_z_val_lo <- unlist(map(1:nrow(new_data), ~ mean(z(unlist(new_ints[.x]), "new", "low"))))
new_z_val_y_lo <- unlist(map(1:length(new_ints_y), ~ mean(new_z_val_lo[unlist(new_ints_y[.x])])))

# 97.5 % credible interval
new_z_val_hi <- unlist(map(1:nrow(new_data), ~ mean(z(unlist(new_ints[.x]), "new", "high"))))
new_z_val_y_hi <- unlist(map(1:length(new_ints_y), ~ mean(new_z_val_hi[unlist(new_ints_y[.x])])))

# Create initial dfs
m_bio_df <- data.frame("year" = 1:26) # original

new_m_bio_df <- data.frame("year" = 1:26) # new

# Fill dfs; set missing years to previous val
fill <- function(model, cred) {
  data <- c(rep(NA, 26))
  
  if (model == "new") {
    z_val_y <- switch(cred,
                      base = new_z_val_y,
                      low = new_z_val_y_lo,
                      high = new_z_val_y_hi)
    
  } else {
    z_val_y <- switch(cred,
                      base = z_val_y,
                      low = z_val_y_lo,
                      high = z_val_y_hi)
  }
  # Fill vals of years with data present
  data[1:7] <- z_val_y[1:7]
  data[9] <- z_val_y[8]
  data[11:13] <- z_val_y[9:11]
  data[15:26] <- z_val_y[12:23]
  
  # Fill vals of years with data absent
  data[8] <- data[7]
  data[10] <- data[9]
  data[14] <- data[13]
  return(data)
} # helper function

new_m_bio_df$m_bio <- fill("new", "base")
new_m_bio_df$m_bio_lo <- fill("new", "low")
new_m_bio_df$m_bio_hi <- fill("new", "high")

m_bio_df$m_bio <- fill("old", "base")
m_bio_df$m_bio_lo <- fill("old", "low")
m_bio_df$m_bio_hi <- fill("old", "high")
  
 # ---- Visualization ----
# Custom palette to match colors from paper
my_pal <- c("#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9",
            "#74add1", "#4575b4")

## Recreate Fig. 2a
tot_bm_data <- mutate(new_data, bm_p_day = biomass / (to.daynr - from.daynr)) %>%
  select(year, bm_p_day)

ggplot(tot_bm_data, aes(x = factor(year, levels = seq(1989, 2014)),
                        y = bm_p_day,
                        fill = year)) +
  geom_boxplot(outlier.alpha = .4, outlier.shape = 1) +
  scale_fill_gradient2(low = "#4575b4",
                        mid = "#fee090",
                        high = "#f46d43",
                        na.value = "grey50",
                        midpoint = 2005,
                        guide = "none") +
  geom_abline(intercept = parms_mean$g_intcp,
              slope = parms_mean$log.lambda) +
  geom_line(aes(y = m_bio, x = year), data = m_bio_df, color = "grey",
            inherit.aes = F) +
  geom_line(aes(y = m_bio_lo, x = year), data = m_bio_df, color = "grey",
            linetype = 2, inherit.aes = F) +
  geom_line(aes(y = m_bio_hi, x = year), data = m_bio_df, color = "grey",
            linetype = 2, inherit.aes = F) +
  geom_line(aes(y = m_bio, x = year), data = new_m_bio_df, color = "red",
            inherit.aes = F) +
  geom_line(aes(y = m_bio_lo, x = year), data = new_m_bio_df, color = "red",
            linetype = 2, inherit.aes = F) +
  geom_line(aes(y = m_bio_hi, x = year), data = new_m_bio_df, color = "red",
            linetype = 2, inherit.aes = F) +
  #geom_jitter(width = .1, alpha = .2) + # optional: datapoints
  scale_x_discrete(breaks = seq(1990, 2015, by = 5), drop = F) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.01, .02, .05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)) +
  theme_classic() +
  labs(y = "Biomass [g/d]", x = "Year")


## Recreate Fig. 2b
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
