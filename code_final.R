## Hallmann 2017

library(R2jags)
library(RCurl) # used for loading data from github
library(dplyr)
library(tidyr)
library(purrr)

setwd("C:/Users/Valentina/Desktop/code")
#-------------------------------------------------
# Load data an create the data subset
#-------------------------------------------------


# Read the data
url_1 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s1_data.csv")
url_2 <- getURL("https://raw.githubusercontent.com/leonthoma/Repro/master/Hallmann_s2_data.csv")


# ---- data sets to use ----
data <- read.csv(text = url_1, header = TRUE, sep = ",")
model.frame <- read.csv(text = url_2, header = TRUE, sep = ",")
model.frame$plot<- as.factor(model.frame$plot)

#---- Data wrangling ----

# Filter Biomass data (only first Dataset of every plot)
py <- select(data, c(plot, year)) %>%
  group_by(plot) %>%
  distinct(plot, year) %>%
  arrange(plot, year) # overview plot by year

new_plots <- group_by(py, plot) %>%
  slice(1) # get first sampling year of each plot

# Create new dfs
new_data <- semi_join(data, new_plots, by = c("plot", "year"))
new_model.frame <- semi_join(model.frame, new_plots, by = c("plot", "year"))

#-------------------------------------------------
# JAGS model
#-------------------------------------------------

# Jags model
## NEEDS TO RUN ONCE TO CREATE .JAGS FILE
## SET WD TO APROPRIATE DIRECTORY
{
  sink("BasicModel.jag")
  cat("model{
  ## Likelihood function for the latent expected daily biomass
  for (i in 1:n) {
  m_bio[i] ~ dnorm(sum(z[tau1[i]:tau2[i]]), sig_sq[i])
  sig_sq[i] <- 1/Var[i]
  Var[i] <- sum(vr[tau1[i]:tau2[i]])
  }

  ## Likelihood function for muHat, it's dependent function and variance
  for (i in 1:ndaily) {
  z[i] <- exp(y[i])
  y[i] <- g_intcp + log.lambda * year[i] + c[1] * daynr[i] + c[2] * daynr2[i] +
    c[3] * daynr[i] * year[i] + c[4] * daynr2[i] * year[i] + b[loctype[i]] +
    eps[plot[i]]
  vr[i] <- exp(2 * y[i] + lvar) * (exp(lvar) - 1)
  }

  ## Priors
  g_intcp ~ dnorm(0, .01)
  log.lambda ~ dnorm(0, .01)
  b[1] <- 0
  for( i in 2:3) {b[i] ~ dnorm(0, .01)}
  for( i in 1:4) {c[i] ~ dnorm(0, .01)}
  sdhat ~ dunif(0, 5)
  lvar <- pow(sdhat, 2)
  for (i in 1:nrandom) {
  eps[i] ~ dnorm(0, tau.re)
  }
  tau.re <- pow(sd.re, -2)
  sd.re ~ dunif(0, 1)
  }
  ")
  sink(NULL)
}

#-------------------------------------------------
# MODEL SET UP AND RUN 
#-------------------------------------------------

#----- MODEL SETUP ------
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

jagsdataBasic <- list(
  m_bio = data$biomass,
  tau1 = with(model.frame, tapply(1:nrow(model.frame), potID, min)),
  tau2 = with(model.frame, tapply(1:nrow(model.frame), potID, max)),
  plot = as.numeric(model.frame$plot),
  loctype = as.numeric(data$location.type[match(model.frame$potID,data$potID)]),
  daynr = as.numeric((model.frame$daynr-mean(data$mean.daynr)) / sd(data$mean.daynr)),
  daynr2 = as.numeric((model.frame$daynr-mean(data$mean.daynr)) / sd(data$mean.daynr))^2,
  year = model.frame$year,
  ndaily = nrow(model.frame),
  n = nrow(data),
  nrandom = max(as.numeric(model.frame$plot))
)

  parametersBasic <- c("g_intcp", "log.lambda", "b", "c", "eps", "sdhat", "sd.re")
  
  #----- MODEL RUN ------ 
 # newjagsmodBasic <- jags(newjagsdataBasic, inits = NULL, parameters = parametersBasic,
 #                     "BasicModel.jag", n.iter = 24000, n.burnin = 4000,
 #                     n.chains = 3, n.thin = 10)
 # 
 # newjagsmodBasic
 # 
 # jagsmodBasic <- jags(jagsdataBasic, inits = NULL, parametersBasic,
 #                      "BasicModel.jag", n.iter = 24000, n.burnin = 4000,
 #                      n.chains = 3, n.thin = 10)
 # 
 # 
 # jagsmodBasic

  #-------------------------------------------------
# LATENT DAILY BIOMASS CALCULATION FOR BOTH DATASETS
#-------------------------------------------------
  
settings<- function(x){
  if( x == "old"){
    print("old dataframe")
    model_used <<- jagsmodBasic
    jagsdata_used <<- jagsdataBasic
    data_used <<- data
    old_bm <- NA
  }
  if(x == "new" ){
    print("new dataframe")
    model_used <<- newjagsmodBasic
    jagsdata_used <<- newjagsdataBasic
    data_used <<- new_data
    new_bm <- NA
  }
  
  #extracting mean of every parameter
  parms_mean <- model_used$BUGSoutput$mean
  
  # 2.4 % credible interval
  parms_mean_lo <- list("b" = as.numeric(model_used$BUGSoutput$summary[1:3, 3]),
                        "c" = as.numeric(model_used$BUGSoutput$summary[4:7, 3]),
                        "deviance" = model_used$BUGSoutput$summary[8, 3],
                        "eps" = as.numeric(model_used$BUGSoutput$summary[9:71, 3]),
                        "g_intcp" = model_used$BUGSoutput$summary[72, 3],
                        "log.lambda" = model_used$BUGSoutput$summary[73, 3],
                        "sd.re" = model_used$BUGSoutput$summary[74, 3],
                        "sdhat" = model_used$BUGSoutput$summary[75, 3])
  
  # 97.5 % credible interval
  parms_mean_hi <- list("b" = as.numeric(model_used$BUGSoutput$summary[1:3, 7]),
                        "c" = as.numeric(model_used$BUGSoutput$summary[4:7, 7]),
                        "deviance" = model_used$BUGSoutput$summary[8, 7],
                        "eps" = as.numeric(model_used$BUGSoutput$summary[9:71, 7]),
                        "g_intcp" = model_used$BUGSoutput$summary[72, 7],
                        "log.lambda" = model_used$BUGSoutput$summary[73, 7],
                        "sd.re" = model_used$BUGSoutput$summary[74, 7],
                        "sdhat" = model_used$BUGSoutput$summary[75, 7])
  
  
  # Set tau & intervals of days
  t_1 <- as.vector(jagsdata_used$tau1[1:nrow(data_used)])
  t_2 <- as.vector(jagsdata_used$tau2[1:nrow(data_used)])
  ints <- map(1:nrow(data_used), ~ seq(t_1[.x], t_2[.x]))
  
  # Set year intervals
  t_1_y <- with(data_used, tapply(1:nrow(data_used), year, min))
  t_2_y <- with(data_used, tapply(1:nrow(data_used), year, max))
  ints_y <- map(1:length(t_1_y), ~ seq(t_1_y[.x], t_2_y[.x]))
  
  
  
  # Helper functions
  y <- function(x, type) {
    parms_mean <- switch(type,
                         base = parms_mean,
                         low = parms_mean_lo,
                         high = parms_mean_hi)
    
    parms_mean$g_intcp + parms_mean$log.lambda * jagsdata_used$year[x] +
      parms_mean$c[1] * jagsdata_used$daynr[x] + 
      parms_mean$c[2] * jagsdata_used$daynr2[x] + 
      parms_mean$c[3] * jagsdata_used$daynr[x] * jagsdata_used$year[x] + 
      parms_mean$c[4] * jagsdata_used$daynr2[x] *jagsdata_used$year[x] + 
      parms_mean$b[jagsdata_used$loctype[x]] +parms_mean$eps[jagsdata_used$plot[x]]
  }
  
  z <- function(x, ...) exp(y(x, ...))
  
  # Calculate predicted biomass
  #z_val_sum <- unlist(map(1:nrow(data_used), ~ sum(z(unlist(ints[.x]), "base"))))
  z_val <- unlist(map(1:nrow(data_used), ~ mean(z(unlist(ints[.x]), "base"))))
  z_val_y <- unlist(map(1:length(ints_y), ~ mean(z_val[unlist(ints_y[.x])])))
  
  # 2.5 % credible interval
  z_val_lo <- unlist(map(1:nrow(data_used), ~ mean(z(unlist(ints[.x]), "low"))))
  z_val_y_lo <- unlist(map(1:length(ints_y), ~ mean(z_val_lo[unlist(ints_y[.x])])))
  
  # 97.5 % credible interval
  z_val_hi <- unlist(map(1:nrow(data_used), ~ mean(z(unlist(ints[.x]), "high"))))
  z_val_y_hi <- unlist(map(1:length(ints_y), ~ mean(z_val_hi[unlist(ints_y[.x])])))
  
  
  m_bio_df <- data.frame("year" = 1:26, "m_bio" = NA, "m_bio_lo" = NA,
                         "m_bio_hi" = NA) # Create initial df
  
  # Fill df; set missing years to previous val
  fill <- function(type) {
    coln <- switch(type,
                   base = 2,
                   low = 3,
                   high = 4)
    
    z_val_y <- switch(type,
                      base = z_val_y,
                      low = z_val_y_lo,
                      high = z_val_y_hi)
    
    m_bio_df[1:7, coln] <<- z_val_y[1:7]
    m_bio_df[9, coln] <<- z_val_y[8]
    m_bio_df[11:13, coln] <<- z_val_y[9:11]
    m_bio_df[15:26, coln] <<- z_val_y[12:23]
    
    m_bio_df[8, coln] <<- m_bio_df[7, coln]
    m_bio_df[10, coln] <<- m_bio_df[9, coln]
    m_bio_df[14, coln] <<- m_bio_df[13, coln]
    
  } # helper function
  
  fill("base")
  fill("low")
  fill("high")
  
  if(x=="old"){
    old_bm <<- m_bio_df
    z_old<<- z_val
  }
  if(x=="new"){
    new_bm<<- m_bio_df
    z_new<<- z_val
  }
  
}

settings(x = "old")
settings(x = "new")

#-------------------------------------------------
# DIAGNOSTICS
#-------------------------------------------------

library(R2jags)
library(coda)
library(mcmcplots)

new_basic_mcmc<- as.mcmc(newjagsmodBasic)
old_basic_mcmc<- as.mcmc(jagsmodBasic)

new_basic_out<- summary(new_basic_mcmc)
old_basic_out<- summary(old_basic_mcmc)

#DIAGNOSTICS FOR ALL PARAMETERS
mcmcplot(new_basic_mcmc)
mcmcplot(old_basic_mcmc)


lambda_old<- old_basic_out$statistics[73,]
lambda_new<- new_basic_out$statistics[73,]

int_old<- old_basic_out$statistics[72,]
int_new<- new_basic_out$statistics[72,]
 

#-------------------------------------------------
# PLOTS 
#-------------------------------------------------

library(ggplot2)

# ----- DENSITY PLOTS TO CHECK THE DATA DISTRIBUTION-----
tot_bm_data <- mutate(new_data, bm_p_day = biomass / (to.daynr - from.daynr)) %>%
  select(year, bm_p_day)%>%
  mutate(model_data = z_new)


tot_bm_data_old <- mutate(data, bm_p_day = biomass / (to.daynr - from.daynr)) %>%
  select(year, bm_p_day)%>%
  mutate(model_data = z_old)

ggplot(tot_bm_data_old)+
  geom_density(aes(bm_p_day), col ="#1B49E7", size =1)+
  geom_density(aes(model_data), col= "#1B49E7", linetype = 2, size =1)


ggplot(tot_bm_data)+
  geom_density(aes(bm_p_day), col ="#E50235", size =1)+
  geom_density(aes(model_data), col= "#E50235", linetype = 2, size =1)


biomass_tot <-rbind((tot_bm_data_old%>%mutate(df_used = "original")),tot_bm_data%>%mutate(df_used = "fixed"))%>%
  gather(bm_p_day,model_data,key = "source",value = "biomass")

colores<- c("#486deb", "#e9315b")

ggplot(biomass_tot)+
  facet_wrap(~df_used)+
  geom_density(aes(biomass, col = source, fill = source), size = 1, alpha = 0.001)+
  scale_color_manual(values = colores)+
  scale_fill_manual(values =colores)




# ----- RECREATE THE PAPER PLOTS -----

# Custom palette to match colors from paper
my_pal <- c("#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9",
            "#74add1", "#4575b4")

biomass_df<- rbind((old_bm%>%mutate(df_used = "original")), (new_bm%>%mutate(df_used = "fixed")))

## Recreate Fig. 2a
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
  geom_line(aes(y = m_bio, x = year, col = df_used), data = biomass_df, inherit.aes = F, size =1 ) +
  geom_line(aes(y = m_bio_lo, x = year, col = df_used), data = biomass_df, linetype = 2 ,inherit.aes = F, size =1) +
  geom_line(aes(y = m_bio_hi, x = year, col = df_used), data = biomass_df, linetype = 2, inherit.aes = F, size =1) +
  scale_color_manual(name = "data frame used", values = colores)+
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

#------PLOTS WITH THE AGGRERATED BIOMASS-------

new_mydata<- new_data%>% 
  mutate(o_biomass = biomass / (to.daynr - from.daynr))%>%
  mutate(model_data = z_new)
new_mydata$exposure.d <- new_mydata$to.daynr - new_mydata$from.daynr
new_mydata$m_biomass<- new_mydata$model_data*new_mydata$exposure.d           

mydata <- data
mydata$exposure.d <- mydata$to.daynr - mydata$from.daynr
mydata<- mutate(mydata, model_data = z_old)%>%
  mutate(mydata, m_biomass = z_old*exposure.d)%>%
  mutate(mydata, o_biomass = biomass/exposure.d)


# DENSITY PLOTS FOR THE OLD AND NEW DATASET
ggplot(new_mydata)+
  geom_density(aes(biomass), col = "blue", size = 1)+
  geom_density(aes(m_biomass), col = "red", linetype =2, size = 1)

ggplot(mydata)+
  geom_density(aes(biomass), col = "blue", size = 1)+
  geom_density(aes(m_biomass), col = "red", linetype =2, size = 1)


# VARIATIONS ON EXPOSURE DAY FOR THE ORIGINAL DATASET
ggplot(mydata)+
  geom_boxplot(aes(year, exposure.d, group = year, fill = year))+
  scale_fill_gradient2(low = "#4575b4", mid = "#fee090", high = "#f46d43",midpoint = 2005, na.value = "grey50", guide = "none")+
  labs(y = "Exposure[d]", x = "Year", title = "Exposure time of traps")



#-------------------------------------------------
# DECAY RATE CALCULATION 
#-------------------------------------------------

# overall decline calculated with log.lambda
l_new<- c(lambda_new[1]-lambda_new[2],lambda_new[1],lambda_new[1]+lambda_new[2])
l_old<- c(lambda_old[1]-lambda_old[2],lambda_old[1],lambda_old[1]+lambda_old[2])

#decrease calculated using the geometric progresion formula
#n = 26 because we take the year 1 as inicial value (a0)
d_new<-((1+l_new)^26)-1
d_old<-((1+l_old)^26)-1

#only spring-summer period 
mydata_fil<- mydata%>% 
  filter(from.daynr >= 91 & to.daynr <= 303)%>%
  group_by(year)%>%
  summarise(min_bm = min(model_data), max_bm = max(model_data), mean_bm = mean(model_data))


new_mydata_fil<- new_mydata%>% 
  filter(from.daynr >= 91 & to.daynr <= 303)%>%
  group_by(year)%>%
  summarise(min_bm = min(model_data), max_bm = max(model_data), mean_bm = mean(model_data))


r_mean_old<- (mydata_fil$mean_bm[25] - mydata_fil$mean_bm[1])/mydata_fil$mean_bm[1]
#r_min<- (mydata_fil$min_bm[25] - mydata_fil$min_bm[1])/mydata_fil$min_bm[1]
#r_max<- (mydata_fil$max_bm[25] - mydata_fil$max_bm[1])/mydata_fil$max_bm[1]

r_mean_new<- (new_mydata_fil$mean_bm[23] - new_mydata_fil$mean_bm[1])/new_mydata_fil$mean_bm[1]
#r_min_n<- (new_mydata_fil$min_bm[23] - new_mydata_fil$min_bm[1])/new_mydata_fil$min_bm[1]
#r_max_n<- (new_mydata_fil$max_bm[23] - new_mydata_fil$max_bm[1])/new_mydata_fil$max_bm[1]

#only summer period
s_old<-mydata%>% 
  filter(from.daynr >= 172 & to.daynr <= 264)%>%
  group_by(year)%>%
  summarise(so_mean = mean(model_data), so_min = min(model_data), so_max = max(model_data))

y27_old<-s_old$so_mean[25]
y1_old<-s_old$so_mean[1]

(y27_old - y1_old)/y1_old

s_new<-new_mydata%>% 
  filter(from.daynr >= 172 & to.daynr <= 264)%>%
  group_by(year)%>%
  summarise(so_mean = mean(model_data), so_min = min(model_data), so_max = max(model_data))

y27_new<-s_new$so_mean[23]
y1_new<-s_new$so_mean[1]

(y27_new - y1_new)/y1_new

