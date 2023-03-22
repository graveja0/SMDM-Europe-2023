#############################################
## Rework mortality risk using Gompertz Model
#############################################

library(ggplot2); theme_set(theme_bw())
library(flexsurv) # For pgompertz
library(survival) # For fitting
library(here)

# Read in Raw Lifetable Data
deaths.male <- read.table(here::here("gompertz-mortality/data/mltper_1x1.txt"),skip=2,header=T,as.is=TRUE); head(deaths.male)
deaths.female <- read.table(here::here("gompertz-mortality/data/fltper_1x1.txt"),skip=2,header=T,as.is=TRUE); head(deaths.female)

# Make Some Edits
deaths.male[which(deaths.male$Age=="110+"),]$Age <- 110
deaths.female[which(deaths.female$Age=="110+"),]$Age <- 110
deaths.male$Age <- as.numeric(deaths.male$Age)
deaths.female$Age <- as.numeric(deaths.female$Age)

# Isolate each file to 2010 deaths only
deaths.male.2010 <- subset(deaths.male,Year==2010); head(deaths.male.2010)
deaths.female.2010 <- subset(deaths.female,Year==2010); head(deaths.female.2010)

# Generate the CDF ## CDF == 1-Survival = 1 - lx/10)
deaths.male.2010$cdf <- with(deaths.male.2010,1-lx/100000)
deaths.female.2010$cdf <- with(deaths.female.2010,1-lx/100000)

# Expand to full data frame
deaths.male.long <- data.frame(Age=rep(deaths.male.2010$Age, deaths.male.2010$dx))
deaths.male.long$Death <- 1
head(deaths.male.long)

deaths.female.long <- data.frame(Age=rep(deaths.female.2010$Age, deaths.female.2010$dx))
deaths.female.long$Death <- 1
head(deaths.female.long)

parameters = data.frame(
  Age = 0:110,
  male.shape=rep(NA,111),
  male.rate=rep(NA,111),
  female.shape=rep(NA,111),
  female.rate=rep(NA,111)
)

for(age in 0:109)
{
  print(age)
  surv.data <- with(deaths.male.long[deaths.male.long$Age >= age,], Surv(Age, Death, origin=age))
  surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")

  parameters$male.shape[age+1] <- surv.model$coefficients[1]
  parameters$male.rate[age+1] <- exp(surv.model$coefficients[2])

  surv.data <- with(deaths.female.long[deaths.female.long$Age >= age,], Surv(Age, Death, origin=age))
  surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")

  parameters$female.shape[age+1] <- surv.model$coefficients[1]
  parameters$female.rate[age+1] <- exp(surv.model$coefficients[2])
}

parameters <- parameters[1:110,]

library(tidyverse)
parameters %>%
  as_tibble() %>%
  mutate(a = male.shape, b = male.rate, t = 0) %>%
  mutate(r_death = (b/a) * (exp(a*t) - exp(a*(t-1)))) %>%
  write.csv(here::here("gompertz-mortality/data/gompertz-mortality.csv"), row.names=FALSE)

#  Let's check the post 40 fit.
surv.data <- deaths.male.long[deaths.male.long$Age >= 40,]

hist(surv.data$Age, freq=FALSE, main="Census Data from 2010", xlab="Male Deaths Starting at 40", breaks=20)
curve(dgompertz(x-40, parameters$male.shape[41], parameters$male.rate[41]) , add=TRUE, col='red')

#  Let's check the post 90 fit.
surv.data <- deaths.male.long[deaths.male.long$Age >= 90,]

hist(surv.data$Age, freq=FALSE, main="Census Data from 2010", xlab="Male Deaths Starting at 90", breaks=20)
curve(dgompertz(x-90, parameters$male.shape[91], parameters$male.rate[91]) , add=TRUE, col='red')
