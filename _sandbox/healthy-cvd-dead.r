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
max_age =  # Max age in life table
    99
min_age =  # Min age in life table
    0
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

# Cause-Deleted Mortality 

ihme_cvd <- 
    tibble::tribble(
        ~age_name,        ~val,
        1L, 0.038771524,
        5L, 0.038546046,
        10L, 0.044403585,
        15L, 0.033781126,
        20L, 0.035856165,
        25L, 0.053077797,
        30L, 0.086001439,
        35L, 0.130326551,
        40L, 0.184310334,
        45L,  0.21839762,
        50L, 0.243705394,
        55L, 0.256334637,
        60L,  0.26828001,
        65L, 0.272698709,
        70L,  0.28529754,
        75L, 0.310642009,
        0L, 0.016750489,
        80L, 0.353518012,
        85L, 0.399856716,
        90L, 0.447817792,
        95L, 0.495305502
    ) %>% 
    mutate(age_ihme = cut(age_name,unique(c(0,1,seq(0,95,5),105)),right=FALSE))  %>% 
    select(age_ihme,  pct_cvd = val) 

lt_ <-
    lt %>% 
    mutate(age_ihme = cut(age,unique(c(0,1,seq(0,95,5),105)),right=FALSE)) %>% 
    left_join(ihme_cvd,"age_ihme") %>%
    mutate(dx_i = round(dx * pct_cvd)) %>% 
    select(age_ihme,
           age,
           D = dx,  # Deaths
           Di = dx_i, # Cause-specific deaths
           lx = lx) %>% # Living
    mutate(a = ifelse(age_ihme == "[0,1)", 0.152, 0.5)) %>% 
    
    # The conditional probability of dying of a given cause given survival to 
    # the age group is easy to obtain, we just multiply the overall probability 
    # by the ratio of deaths of a given cause to all deaths.
    
    mutate(q = edit.na(1 - lead(lx)/lx, 1),
           qi = q * Di/D) %>% 
    
    # The unconditional counts of deaths of any cause and of a given cause 
    # are calculated multiplying by the number surviving to the start of each 
    # age group, which is lx. Recall that to die of cause i in the interval 
    # [x, x+n) one must survive all causes up to age x.
    
    mutate(d = lx * q, 
           di = lx * qi) %>% 
    
    # In preparation for the next part, note that if we had nmx and we were willing 
    # to assume that the hazard is constant in each age group we would have had a 
    # slightly different estimate of the survival function. Let us “back out” 
    # the rates from the probabilities:
    
    mutate(n = c(diff(age),NA), 
           m =  edit.na( q/(n - q * (n - a)), 1/tail(a,1))) %>%  # m[last] = 1/a[last]
    
    # With these rates we compute the cumulative hazard and survival as
    
    mutate(H = cumsum(n * m), 
           S = edit.na(exp(-lag(H)), 1)) %>%  # S[1] = 1
    
    # We compute cause-specific rates by dividing deaths of a given cause into person-years 
    # of exposure, which is equivalent to multiplying the overall rate by the ratio of 
    # deaths of a given cause to the total. Here we want deaths for causes other than 
    # neoplasms. I will use the subscript d for deleted:
    
    mutate(Rd = (D - Di)/D,
           md = m * Rd) %>% 
    
    # We compute the conditional probability of surviving an age group after 
    # deleting a cause as the overall probability raised to Rd, and then calculate 
    # the survival function as a cumulative product
    
    mutate(pd = (1 - q)^Rd,
           ld = 100000 * cumprod(c(1, pd[-length(pd)]))) %>% 
    
    # Then we construct a survival function in the usual way, but treating this 
    # hazard as if it was the only one operating:
    
    mutate(Hd = cumsum(n * md), 
           Sd = edit.na(exp(-lag(Hd)), 1)) %>% # Sd[1] = 1 
    mutate(Pd =  edit.na((Sd - lead(Sd))/md, tail(Sd/md, 1))) %>% 
    
    # Now do it for cause-specific death. 
    
    mutate(Ri = Di / D, 
           pi = (1 - qi)^Ri,
           li = 100000 * cumprod(c(1, pi[-length(pi)]))) %>% 
    mutate(mi = m - md)

###########################
# Cause-Specific Mortality 
###########################
ages_     <- lt_$age[lt_$age<=max_age & lt_$age>=min_age]
deaths_   <- lt_$d[lt_$age<=max_age & lt_$age>=min_age] - lt_$di[lt_$age<=max_age & lt_$age>=min_age]
exposure_  <- lt_$lx[lt_$age<=max_age & lt_$age>=min_age]

mort_fit_CVDdeleted <- MortalityLaw(
    x  = ages_,
    Dx  = deaths_,   # vector with death counts
    Ex  = exposure_, # vector containing exposures
    law = "HP2",
    opt.method = "LF2")

######################
# Parameterize model
######################
params <-
    list(
        min_age = 50,
        n_cycles = 100,
        time_step = 1,
        r_H_CVD = 0.1,
        u_H = 1, 
        u_CVD = 1,
        u_D = 0,
        u_DCVD = 0,
        c_H = 0, 
        c_CVD = 0,
        c_DCVD = 0,
        c_D = 0, 
        background_mortality = coef(mort_fit_CVDdeleted),
        cause_specific_mortality = approxfun(lt_$age, lt_$mi,rule = 2)
    )

params <- modifyList(params, list(
    payoff_qaly = c("healthy" = params$u_H, "cvd" = params$u_CVD, "cvddeath" = params$u_DCVD, "dead" = params$u_D),
    payoff_cost = c("healthy" = params$c_H, "cvd" = params$c_CVD, "cvddeath" = params$c_DCVD, "dead" = params$c_D)
))

m_Qt_markov <- function(t, h = params$time_step)
{
    lapply(t, function(tt){
        current_age <- params$min_age  + (tt)*h - 1; current_age
        r_death = HP2(current_age, params$background_mortality)$hx; r_death
        r_cvd <- params$cause_specific_mortality(current_age) ; r_cvd
        
        m_Q <- 
            matrix(c(
            0, params$r_H_CVD , 0, r_death ,
            0, 0, r_cvd , r_death ,
            0, 0, 0, 0,
            0, 0, 0, 0),
            nrow = 4,
            byrow= TRUE,
            dimnames = list(names(params$payoff_qaly),names(params$payoff_qaly))
        )
        diag(m_Q) <- -rowSums(m_Q)
        m_Q
        
    })
}

m_Qt_non_markov <- function(t,acc = c("accCVD"), tunnel = c("T1CVD"))
{
    markov     <- m_Qt_markov(t)
    non_markov <- list()
    for(i in t)
    {
        r_ <- nrow(markov[[i]])
        add_ <- length(acc) + length(tunnel) + 1 
        
        # Expand matrix
        non_markov[[i]] <- cbind(markov[[i]],     matrix(rep(0, r_ * add_), nrow=r_))
        non_markov[[i]] <- rbind(non_markov[[i]], matrix(rep(0, add_ * (r_+add_)), nrow=add_))
    }
    
    
    lapply(non_markov, function(m){
        # Put in State Names
        rownames(m) <- c(rownames(m)[-which(rownames(m)=="")], acc, tunnel, "N")
        colnames(m) <- c(colnames(m)[-which(colnames(m)=="")], acc, tunnel, "N")
        
        # Define Accumulator
        m["healthy", "accCVD"] <- m["healthy", "cvd"] 
        
        # Define Tunnel state entry
        m["healthy", "T1CVD"]  <- m["healthy", "cvd"]
        
        m  # Note: Tunnel states are not fully defined at this point.
    })
    
}

m_Pt_simple <- function(t,h = params$time_step) {
    m_P <- # Embed the matrices into the timesteps
        lapply(m_Qt_markov(t), function(m) expm::expm(m * h))

}

m_Pt_fn <- function(t,h = params$time_step) {
    m_P_ <- lapply(m_Qt_non_markov(t), function(m) expm::expm(m * h))
    
    lapply(m_P_, function(m)
    {
        # It is possible to exit tunnel to external risk of death
        m["T1CVD", "N"]  <- m["healthy", "cvddeath"] + m["healthy", "dead"]
        m["T1CVD", "T1CVD"] <- 0  # Cannot remain in tunnel
        
        # last state in the tunnel is a terminal state
        m["T1CVD", "N"]  <- 1
        
        # Note: At this point, the "N" state could be stripped as it was
        #       only required for the embedding, and serves no other purpose
        #       at this point
        
        states <- colnames(m)[-which(colnames(m)=="N")]
        m[states, states]
    })
}

params$m_P <- m_Pt_fn(1:params$n_cycles, h = params$time_step)

################
# Simulation
################

sim_alive_dead <- function(params) {
    tr_ <- t(c("healthy" = 1, "cvd"  = 0,   "cvddeath"   = 0,   "dead"   = 0, "accCVD"   = 0,  "T1CVD" = 0))
    tr <- 
        do.call(rbind,lapply(params$m_P, function(tp) {
            tr_ <<- tr_ %*% tp
        }))
    tr <- rbind(t(c(1,0,0,0,0,0)),tr)
    return(tr)
}

tr_alive_dead <- 
    sim_alive_dead(params)

###########
# Payoffs
###########

life_exp <- 
    cycle_adj(tr_alive_dead[,c("healthy","cvd","cvddeath","dead")] %*% params$payoff_qaly, 1)
life_exp


