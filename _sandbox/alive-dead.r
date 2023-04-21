###################
# Alive-Dead Model
###################

library(tidyverse)
library(demography)
library(MortalityLaws)
options("scipen" = 100, "digits" = 5)

# Functions
edit.na <- function(x, value) { x[is.na(x)] <- value; x}
alt_simp_coef <- function(i) c(17, 59, 43, 49, rep(48, i-8), 49, 43, 59, 17) / 48
cycle_adj      <- function(x,h) h*sum(alt_simp_coef(length(x)) * x)


#####################
# Mortality Modeling
#####################

mortality_year = 2019
radix = 100000

lt_usa_file <- "https://github.com/graveja0/SMDM-Europe-2023/raw/main/_learnr/smdm-europe-2023-cvd-model/www/usa-life-table.rds"

lt <- 
    readRDS(url(lt_usa_file)) %>% 
    demography::lifetable(.,series = "total", years = mortality_year) %>% 
    as_tibble() %>% 
    mutate_at(vars(lx,dx), function(x) x * radix) %>% 
    mutate(country = "USA") %>% 
    mutate(age = x)

lt_usa_file <- "https://github.com/graveja0/SMDM-Europe-2023/raw/main/_learnr/smdm-europe-2023-cvd-model/www/usa-life-table.rds"

lt <- 
    readRDS(url(lt_usa_file)) %>% 
    demography::lifetable(.,series = "total", years = mortality_year) %>% 
    as_tibble() %>% 
    mutate_at(vars(lx,dx), function(x) x * radix) %>% 
    mutate(country = "USA") %>% 
    mutate(age = x)

min_age = 50
max_age = 99

ages     <- lt$x[lt$x<=max_age & lt$x>=min_age]
deaths   <- lt$dx[lt$x<=max_age & lt$x>=min_age]
exposure <- lt$lx[lt$x<=max_age & lt$x>=min_age]

mort_fit <- MortalityLaw(
    x  = ages,
    Dx  = deaths,   # vector with death counts
    Ex  = exposure, # vector containing exposures
    law = "HP2",
    opt.method = "LF2")

######################
# Parameterize model
######################
params <-
    list(
        n_cycles = 100,
        time_step = 1,
        background_mortality = coef(mort_fit),
        payoff_qaly = c("alive" = 1, "dead" = 0),
        payoff_cost = c("alive" = 0, "dead" = 0)
    )

m_Pt_fn <- function(t, h = 1)
{
    lapply(t, function(tt){
        r_death<- HP2(min_age  + (tt)*h - 1, params$background_mortality)$hx
        p_death <- 1 - exp(-r_death * h)
        
        matrix(c(1-p_death,p_death,0,0), 
               byrow = TRUE,
               nrow=2,
               dimnames = list(c("alive","dead"),c("alive","dead")))
        
    })
}

params$m_P <- m_Pt_fn(1:params$n_cycles, h = params$time_step)

################
# Simulation
################

sim_alive_dead <- function(params) {
    tr_ <- t(c("alive" = 1, "dead" = 0))
    tr <- 
        do.call(rbind,lapply(params$m_P, function(tp) {
            tr_ <<- tr_ %*% tp
        }))
    tr <- rbind(t(c(1,0)),tr)
    return(tr)
}

tr_alive_dead <- sim_alive_dead(params)

###########
# Payoffs
###########

life_exp <- 
    cycle_adj(tr_alive_dead %*% params$payoff_qaly, 1)
life_exp


