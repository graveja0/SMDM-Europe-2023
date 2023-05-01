library(learnr)
library(tidyverse)
library(demography)
library(MortalityLaws)
library(directlabels)
library(ggsci)
library(hrbrthemes)
library(MASS)
library(mgcv)
library(patchwork)
library(knitr)
library(kableExtra)
library(here)
library(ggsci)
library(expm)
library(glue)
select =dplyr::select
options("scipen" = 100, "digits" = 5)

alt_simp_coef <- function(i) c(17, 59, 43, 49, rep(48, i-8), 49, 43, 59, 17) / 48
cycle_adj      <- function(x,h) h*sum(alt_simp_coef(length(x)) * x)


params = list(
    t_names = c("without_drug", "with_drug"),
    n_treatments =2,
    s_names  =c("Asymptomatic_disease", "Progressive_disease", "DeadCause","Dead"),
    n_states =4,
    a_names = c("accProgressive_disease"),
    tunnel_names = c("trProgressive","trDeadCause"),
    
    n_cohort =1000,
    cycle =1,
    
    n_cycles = 46,
    Initial_age =55,
    effect =0.5, # at 0.048 the two approaches will straddle a WTP of 100k
    
    cAsymp =500,
    cDeath =1000,
    cDrug =1000,
    cProg =3000,
    uAsymp =0.95,
    uProg =0.75,
    oDr = 0, #0.06,
    cDr = 0, #0.06,
    tpDcm =0.15,
    tpProg =0.01,
    tpDn =0.0379 
    
)

t = 1:params$n_cycles
tt = params$n_cycles
h = params$cycle
traditional = FALSE

m_Pt <- function(t,h = params$cycle, traditional = FALSE) {
    lapply(t, function(tt){
        
        current_age <- params$Initial_age  + (tt)*h - 1
        cycle = (tt)*h #- 1

        tpDn_lookup <-
            c("(34,44]" = 0.0017,
              "(44,54]" = 0.0044,
              "(54,64]" = 0.0138,
              "(64,74]" = 0.0379,
              "(74,84]" = 0.0912,
              "(84,100]" = 0.1958)
        
        age_grp <- cut(current_age, breaks = c(34,44,54,64,74,84,100))
        tpDn <- tpDn_lookup[age_grp]
        r_death = -log(1-tpDn)
        
        
        #r_death = HP2(current_age, params$background_mortality)$hx
        tpProg_ <- params$tpProg * (cycle )
        tpDcm_ <- params$tpDcm 
        effect_ <- params$effect
        n_states_ <- params$n_states
        s_names_ <- params$s_names
        t_names_ <- params$t_names
        n_treatments_ <- params$n_treatments
        
        #########################
        # Markovian Rate Matrix
        #########################
        m_R_ <- 
            array(data = c(0, 0, 0,0,
                                   -log(1-tpProg_), 0, 0,0,
                                   0,-log(1-tpDcm_),0,0,
                           r_death,  r_death, 0,0,
                           
                                   0, 0, 0,0,
                                   -log(1-(tpProg_*(1-effect_))), 0, 0,0,
                                   0,-log(1-tpDcm_),0,0,
                           r_death,  r_death,0, 0),
                          dim = c(n_states_, n_states_, n_treatments_),
                          dimnames = list(from = s_names_,
                                          to = s_names_,
                                          t_names_))

        m_R_ <- apply(m_R_,3,function(x){
            diag(x) <- -rowSums(x)
            x
        },simplify=FALSE)
        
        ###############################
        # Add Non-Markovian Components
        ###############################
        
        m_R <- lapply(m_R_,function(x){
            r_ <- nrow(x)
            add_ <- length(params$a_names)   + length(params$tunnel_names) + 1
            
            new_row_names <- c(rownames(x), params$a_names,params$tunnel_names,"XX")
            new_col_names <- c(colnames(x), params$a_names,params$tunnel_names,"XX")
            # Expand matrix
            x <- cbind(x,     matrix(rep(0, r_ * add_), nrow=r_))
            x  <- rbind(x, matrix(rep(0, add_ * (r_+add_)), nrow=add_))
            
            rownames(x) <- new_row_names
            colnames(x) <- new_col_names
            
            x["Asymptomatic_disease","accProgressive_disease"] = x["Asymptomatic_disease","Progressive_disease"]
            #x["Asymptomatic_disease","trDead"] = x["Asymptomatic_disease","Dead"]
            x["Asymptomatic_disease","trProgressive"] = x["Asymptomatic_disease","Progressive_disease"]
            # In the original model, everyone in progressive disease who dies of other causes also gets the death costs
            x["Progressive_disease","trDeadCause"] = x["Progressive_disease","DeadCause"] 
            x
        })
        
        m_P_ = lapply(m_R,expm)
        
        if (traditional) {
            m_P_ <- 
                array(data = c(0, 0, 0,0,
                               tpProg_, 0, 0,0,
                               0,tpDcm_+tpDn,0,0,
                               tpDn ,  0, 0,0,
                               
                               0, 0, 0,0,
                               tpProg_*(1-effect_), 0, 0,0,
                               0,tpDcm_+tpDn,0,0,
                               tpDn,  0,0, 0),
                      dim = c(n_states_, n_states_, n_treatments_),
                      dimnames = list(from = s_names_,
                                      to = s_names_,
                                      t_names_))
            diag(m_P_[,,1]) <- 1 - rowSums(m_P_[,,1])
            diag(m_P_[,,2]) <- 1 - rowSums(m_P_[,,2])
            m_P <- list()
            m_P[[1]] <- m_P_[,,1]
            m_P[[2]] <- m_P_[,,2]
            
        } else {
            # Now make final edits to tunnel states
            m_P <- lapply(m_P_,function(x){
                # leave after getting recorded
                x["trDeadCause","trDeadCause"] = 0
                # x["trDead","trDead"] = 0
                x["trProgressive","trProgressive"] = 0
                x[-which(rownames(x)=="XX"),-which(colnames(x)=="XX")]    
            }) 
        }
        

        
        return(m_P)
    })
}

m_P <- m_Pt(1:(params$n_cycles-1))
m_P_trad <- m_Pt(1:(params$n_cycles-1),traditional = TRUE)

# Reproduce Original Results

u_payoff_traditional <- with(params,{
    array(c("Asymptomatic_disease" = uAsymp, "Progressive_disease" = uProg, "DeadCause" = 0, "Dead" = 0, 
            "Asymptomatic_disease" = uAsymp, "Progressive_disease" = uProg,  "DeadCause" = 0, "Dead" = 0 ),
          dim = c(1, n_states, n_treatments),
          dimnames = list(from = "cost",
                          to = c(s_names),
                          t_names))
})

c_payoff_traditional <- with(params,{
    array(c("Asymptomatic_disease" = cAsymp, "Progressive_disease" = cProg, "DeadCause" = 0 , "Dead" = 0, 
            "Asymptomatic_disease" = cAsymp+cDrug, "Progressive_disease" = cProg, "DeadCause" = 0 , "Dead" = 0 ),
          dim = c(1, n_states, n_treatments),
          dimnames = list(from = "cost",
                          to = c(s_names),
                          t_names))
})

c_transition_payoff_traditional <- with(params,{
    array(c("Asymptomatic_disease" = 0, "Progressive_disease" = 0, "DeadCause" = cDeath , "Dead" = 0, 
            "Asymptomatic_disease" = 0, "Progressive_disease" = 0, "DeadCause" = cDeath , "Dead" = 0 ),
          dim = c(1, n_states, n_treatments),
          dimnames = list(from = "cost",
                          to = c(s_names),
                          t_names))
})


sim_cohort_traditional <- function(params) {
    
    m_P <- m_Pt(1:(params$n_cycles-1),traditional = TRUE)
    
    tr <- 
        m_P %>% transpose() %>% 
        map(~({
            tr_ <- t(c("Asymptomatic_disease" = 1000, "Progressive_disease" = 0, "DeadCause" = 0, "Dead" = 0))
            do.call(rbind,lapply(.x, function(tp) {
                tr_ <<- tr_ %*% tp
            }))
        }))
    tr <- 
        tr %>% 
        map(~({
            .x <- rbind(t(c(1000,0,0,0)),.x)
        }))
    
    return(tr)
}

markov_trace <- 
    markov_trace_traditional <- sim_cohort_traditional(params)

state_costs_without <- markov_trace[[1]] %*% c_payoff_traditional[,,1]
transition_costs_without <- matrix(c(0,diff(markov_trace[[1]][,"DeadCause"]))) * params$cDeath
costs_without <- state_costs_without + transition_costs_without
qalys_without <- markov_trace[[1]] %*% u_payoff_traditional[,,1]

state_costs_with <- markov_trace[[2]] %*% c_payoff_traditional[,,2]
transition_costs_with <- matrix(c(0,diff(markov_trace[[2]][,"DeadCause"]))) * params$cDeath
costs_with <- state_costs_with + transition_costs_with
qalys_with <- markov_trace[[2]] %*% u_payoff_traditional[,,2]

cycle_adj <- alt_simp_coef(params$n_cycles)
discounting <- 1/(1+params$oDr)^(0:(params$n_cycles-1))

tot_costs_without <- sum(as.vector(costs_without[-1]) * discounting[-1])
tot_costs_with <- sum(as.vector(costs_with[-1])  * discounting[-1])

tot_qalys_without <- sum(as.vector(qalys_without) *  discounting)
tot_qalys_with <- sum(as.vector(qalys_with) *  discounting)

#inc_cost <- (tot_costs_with - tot_costs_without) 
inc_qaly <- (tot_qalys_with - tot_qalys_without); inc_qaly
inc_cost <- (tot_costs_with - tot_costs_without); inc_cost

icer = inc_cost / inc_qaly ; icer

c_payoff <- with(params,{
    array(c("Asymptomatic_disease" = cAsymp, "Progressive_disease" = cProg, "DeadCause" = 0 , "Dead" = 0, "accProgressive_disease" = 0,"trProgressive" = 0, "trDeadCause" = cDeath,
            "Asymptomatic_disease" = cAsymp+cDrug, "Progressive_disease" = cProg, "DeadCause" = 0 , "Dead" = 0, "accProgressive_disease" = 0,"trProgressive" = 0, "trDeadCause" = cDeath),
          dim = c(1, n_states+length(a_names)+length(tunnel_names), n_treatments),
          dimnames = list(from = "cost",
                          to = c(s_names,a_names,tunnel_names),
                          t_names))
    
    })

u_payoff <- with(params,{
    array(c("Asymptomatic_disease" = uAsymp, "Progressive_disease" = uProg, "DeadCause" = 0 , "Dead" = 0, "accProgressive_disease" = 0,"trProgressive" = 0, "trDeadCause" = 0,
            "Asymptomatic_disease" = uAsymp, "Progressive_disease" = uProg, "DeadCause" = 0 , "Dead" = 0, "accProgressive_disease" = 0,"trProgressive" = 0,"trDeadCause" = 0),
          dim = c(1, n_states+length(a_names)+length(tunnel_names), n_treatments),
          dimnames = list(from = "cost",
                          to = c(s_names,a_names,tunnel_names),
                          t_names))
    
})


sim_cohort<- function(params) {
    
    m_P <- m_Pt(1:(params$n_cycles-1),traditional = FALSE)
    
    tr <- 
        m_P %>% transpose() %>% 
        map(~({
            tr_ <- t(c("Asymptomatic_disease" = 1000, "Progressive_disease" = 0, "DeadCause" = 0, "Dead" = 0,
                       "accProgressive_disease" = 0, "trProgressive"=0, "trDeadCause" = 0))
            do.call(rbind,lapply(.x, function(tp) {
                tr_ <<- tr_ %*% tp
            }))
        }))
    tr <- 
        tr %>% 
        map(~({
            .x <- rbind(t(c(1000,0,0,0,0,0,0)),.x)
        }))
    
    return(tr)
}

m_P <- m_Pt(1:params$n_cycles)
markov_trace <- sim_cohort(params)

costs_without <- markov_trace[[1]] %*% c_payoff[,,1]
qalys_without <- markov_trace[[1]] %*% u_payoff[,,1]

costs_with <- markov_trace[[2]] %*% c_payoff[,,2]
qalys_with <- markov_trace[[2]] %*% u_payoff[,,2]

cycle_adj <- alt_simp_coef(params$n_cycles)
discounting <- 1/(1+params$oDr)^(0:(params$n_cycles-1))

tot_costs_without <- sum(as.vector(costs_without[-1]) * discounting[-1])
tot_costs_with <- sum(as.vector(costs_with[-1])  * discounting[-1])

tot_qalys_without <- sum(as.vector(qalys_without) * discounting)
tot_qalys_with <- sum(as.vector(qalys_with) * discounting)

inc_qaly2 <- (tot_qalys_with - tot_qalys_without); inc_qaly
inc_cost2 <- (tot_costs_with - tot_costs_without); inc_cost

icer2 = inc_cost2 / inc_qaly2 ; icer2

icer2/icer

c("traditional" = icer, "embedded" = icer2)




# In round 1, 9.931 people progress, but only 9.1 remain.
# In round 2, an additional 19.392 progress, meaning taht 28.493 ever enter the state with progressed disease, or develop it during the state.
# But by the end of the round, 3.09506 people die out of the progressive state, leaving 25.3982  by the end of the round. 
# If costs only matter at the end of the round, then 25.3982 is the right number to use. 


