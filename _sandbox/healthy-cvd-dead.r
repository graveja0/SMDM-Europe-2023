library(tidyverse)
library(demography)
library(MortalityLaws)
library(directlabels)
library(ggsci)
library(hrbrthemes)
library(MASS)
library(mgcv)
library(patchwork)
select <- dplyr::select
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
cohort_starting_age =  # Min age in life table
    25
mortality_year = 2021
radix = 100000

lt_usa_file <- "https://github.com/graveja0/SMDM-Europe-2023/raw/main/_learnr/smdm-europe-2023-cvd-model/www/usa-life-table.rds"

lt <- 
    readRDS(url(lt_usa_file)) %>% 
    #readRDS("./_sandbox/mortality/usa-life-table.rds") %>% 
    demography::lifetable(.,series = "total", years = mortality_year) %>% 
    as_tibble() %>% 
    mutate_at(vars(lx,dx), function(x) x * radix) %>% 
    mutate(country = "USA") %>% 
    mutate(age = x)

# Cause-Deleted Mortality 

ihme_cvd <-   # Source: https://vizhub.healthdata.org/gbd-results/
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

lt_ <-  # Source: https://grodri.github.io/demography/neoplasms
    lt %>% 
    mutate(age_ihme = cut(age,unique(c(0,1,seq(0,95,5),105)),right=FALSE)) %>% 
    left_join(ihme_cvd,"age_ihme") %>%
    mutate(dx_i = round(dx * pct_cvd)) %>% 
    select(age_ihme,
           age,
           D = dx,  # Deaths
           Di = dx_i, # Cause-specific deaths
           lx = lx) %>% # Living
    mutate(a = ifelse(age_ihme == "[0,1)", 0.152, 0.5))  %>% 
    mutate(age_interval = c(diff(age), NA)) %>% 
    select(age,lx,D,Di,a,age_interval) %>% 
    # Probability of death in interval =  deaths / total living at beginning of interval. 
    mutate(q = replace_na(1 - lead(lx) / lx, 1)) %>%  
    # Cause specific probability of death. 
    mutate(qi = q * Di / D) %>% 
    mutate(m = -log(1-q)/age_interval) %>% 
    # We compute cause-specific rates by dividing deaths of a given cause into person-years
    # of exposure, which is equivalent to multiplying the overall rate by the ratio of
    # deaths of a given cause to the total.
    mutate(Rd = (D - Di) / D) %>% 
    mutate(md = m * Rd) %>% 
    mutate(mi = m - md)

bx <- mutate(lt_, agem = age + age_interval/2, mi = m - md)[-nrow(lt_), ]

p_cd <- 
    bx %>% 
    select(agem,m,md,mi) %>% 
    gather(series,value,-agem) %>% 
    mutate(series = factor(series,levels = c("m","md","mi"),labels = c("All Cause", "Non-CVD","CVD"))) %>% 
    ggplot(aes(x = agem, y = value, colour = series)) + geom_line() + scale_y_log10() + 
    ggsci::scale_color_aaas() + 
    geom_dl(method = "smart.grid",aes(label=series)) + 
    scale_x_continuous(expand = c(0.5,0)) + 
    labs(x = "Age", y = "Mortality Rate") + 
    geom_point(data = lt_ %>% mutate(series = "Life Table") %>% filter(age<100), aes(x = age, y = m), alpha =0.2, colour = "darkblue") + 
    geom_point(data = lt_ %>% mutate(series = "Life Table") %>% filter(age<100), aes(x = age, y = md), alpha =0.2, colour = "red") + 
    geom_point(data = lt_ %>% mutate(series = "Life Table") %>% filter(age<100), aes(x = age, y = m - md), alpha =0.2, colour = "darkgreen") +
    theme_ipsum_tw(base_family = "Arial") +
    theme(legend.position = "none"); p_cd

###########################
# Cause-Specific Mortality 
###########################
ages_     <- lt_$age[lt_$age<=max_age & lt_$age>=cohort_starting_age]
deaths_   <- lt_$D[lt_$age<=max_age & lt_$age>=cohort_starting_age] - lt_$Di[lt_$age<=max_age & lt_$age>=cohort_starting_age]
exposure_  <- lt_$lx[lt_$age<=max_age & lt_$age>=cohort_starting_age]

mort_fit_CVDdeleted <- MortalityLaw(
    x  = ages_,
    Dx  = deaths_,   # vector with death counts
    Ex  = exposure_, # vector containing exposures
    law = "HP2",
    opt.method = "LF2")

plot(mort_fit_CVDdeleted)


# library(splines)
# cvd_fit <- 
#     glm(mi ~ bs(age,knots = c(1,5,10,20,30)), data = lt_ %>% filter(age < max_age),family = gaussian(link="log"),start = c(0,0,0,0,0,0,0,0,0))
# 
# plot(lt_$age,log(lt_$mi),col='red')    
# points(cohort_starting_age:max_age,predict(cvd_fit,newdata = list(age = cohort_starting_age:max_age),type="link"))
# 

cvd_fit <-
    gam(mi ~ s(age), family = gaussian(link="log"), data = lt_ %>% filter(age < max_age))
plot(lt_$age,log(lt_$mi),col='red')
points(cohort_starting_age:max_age,predict(cvd_fit,newdata = tibble(age = cohort_starting_age:max_age),type="link"))


######################
# Parameterize model
######################
params <-
    list(
        cohort_starting_age = cohort_starting_age,
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
        cause_specific_mortality = approxfun(lt_$age, lt_$mi,rule = 2),
        cause_specific_mortality_gam = cvd_fit
    )

params <- modifyList(params, list(
    payoff_qaly = c("healthy" = params$u_H, "cvd" = params$u_CVD, "cvddeath" = params$u_DCVD, "dead" = params$u_D),
    payoff_cost = c("healthy" = params$c_H, "cvd" = params$c_CVD, "cvddeath" = params$c_DCVD, "dead" = params$c_D)
))

m_Qt_markov <- function(t, h = params$time_step)
{
    lapply(t, function(tt){
        current_age <- params$cohort_starting_age  + (tt)*h - 1; current_age
        r_death = HP2(current_age, params$background_mortality)$hx; r_death
        r_cvd <- params$cause_specific_mortality(current_age) ; r_cvd
            # predict(params$cause_specific_mortality_gam,newdata = tibble(age = current_age),type="link") %>% 
            # exp(.)
 
        
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

m_Qt_non_markov <- function(t,acc = c("accCVD"), tunnel = c("t_CVD_1","t_CVD_2"))
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
        m["healthy", "t_CVD_1"]  <- m["healthy", "cvd"]
        
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
        m["t_CVD_1", "N"]  <- m["healthy", "cvddeath"] + m["healthy", "dead"]
        m["t_CVD_1", "t_CVD_1"] <- 0  # Cannot remain in tunnel
        m["t_CVD_1", "t_CVD_2"] <- 1 - m["t_CVD_1", "N"]  # The tunnel is everything else 
        
        # last state in the tunnel is a terminal state
        m["t_CVD_2", "t_CVD_2"] <- 0
        m["t_CVD_2", "N"]  <- 1
        
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

sim_cvd <- function(params) {
    tr_ <- t(c("healthy" = 1, "cvd"  = 0,   "cvddeath"   = 0,   "dead"   = 0, "accCVD"   = 0,  "t_CVD_1" = 0, "t_CVD_2" = 0))
    tr <- 
        do.call(rbind,lapply(params$m_P, function(tp) {
            tr_ <<- tr_ %*% tp
        }))
    tr <- rbind(t(c(1,0,0,0,0,0,0)),tr)
    return(tr)
}

tr_cvd <- 
    sim_cvd(params)

###########
# Payoffs
###########

life_exp <- 
    cycle_adj(tr_cvd[,c("healthy","cvd","cvddeath","dead")] %*% params$payoff_qaly, 1)
life_exp

lt %>% filter(age==cohort_starting_age) %>% pull(ex)

lt_markov <- 
    tr_cvd %>% 
    as_tibble() %>%
    mutate(alive = healthy + cvd) %>% 
    mutate(lx = radix * alive) %>%
    mutate(q = edit.na(1 - lead(lx)/lx, 1)) %>%
    mutate(age = params$cohort_starting_age + (row_number()-1)/1) %>% 
    inner_join(lt %>% filter(age>=cohort_starting_age & age <max_age),"age") %>% 
    select(age, q, qx) %>% 
    gather(source,value,-age) %>% 
    mutate(source = factor(source,levels = c("q","qx"), labels = c("Markov","Life Table")))

p1 <- 
    lt_markov %>% 
    ggplot(aes(x = age, y = value, colour = source)) + #scale_y_log10() + 
    geom_point() +
    hrbrthemes::theme_ipsum_pub(base_family = "Arial") + 
    ggsci::scale_colour_aaas(name="") + 
    theme(legend.position = "top") + 
    labs(x = "Age" , y = "Death Rate")


p2 <- 
    lt_markov %>% 
    ggplot(aes(x = age, y = value, colour = source)) + scale_y_log10() + 
    geom_point() +
    hrbrthemes::theme_ipsum_pub(base_family = "Arial") + 
    ggsci::scale_colour_aaas(name="") + 
    theme(legend.position = "top") + 
    labs(x = "Age" , y = "log(Death Rate)")

p1 + p2
