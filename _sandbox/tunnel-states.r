##################
# 0. Setup
##################

source("slides/manifest.r")

##################
# 1. Parameterize
##################

params = 
    list(
        t_names = c("natural_history"),           # Strategy names. 
        n_treatments = 1,                         # Number of treatments
        s_names  = c("Healthy", "CVD", "Dead"),   # State names
        n_states = 3,                             # Number of states
        n_cohort = 1,                             # Cohort size
        n_cycles = 100,                           # Number of cycles in model.  
        cycle = 1,                                # Cycle length
        initial_age = 55,                         # Cohort starting age
        r_H_CVD = 0.15,                           # Rate of healthy -> CVD
        hr_CVD = 10,                              # Hazard Ratio: CVD Death
        r_H_D = 0.01                              # Rate of healthy -> dead
    )

# Markovian Rate Matrix

params$mR <- 
    with(params,{
        r_CVD_D <- hr_CVD * r_H_D
        R_ <- 
            array(data = c(0, 0, 0,  
                           r_H_CVD,0, 0,
                           r_H_D,r_H_D + r_CVD_D,0,
                           r_H_D, 0,0),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))
        R <- apply(R_,3, simplify=FALSE,function(x) {
            diag(x) = -rowSums(x)
            x * cycle
        })
        R    
    })
R_ <- params$mR[[1]]

R_

# Add tunnel and receiving buckets

R <- 
    cbind(R_,"tunCVDy1" = c(0,0,0), "tunCVDy2" = c(0,0,0),"NULL" = c(0,0,0)) %>% 
    rbind(.,"tunCVDy1" = c(0,0,0,0,0,0), "tunCVDy2" = c(0,0,0,0,0,0), "NULL" = c(0,0,0,0,0,0))  
R

P = expm(R)

# First step is to enter the tunnel
P["Healthy","tunCVDy1"] = P["Healthy","CVD"]
# It is possible to exit the tunnel state to death. 
P["tunCVDy1","Dead"] = P["CVD","Dead"]
# Cannot remain in the tunnel
P["tunCVDy1","tunCVDy1"] = 0
P["tunCVDy1","tunCVDy2"] = 1 - P["CVD","Dead"]
# It is possible to exit the second tunnel state to death
P["tunCVDy2","Dead"] = P["CVD","Dead"]
# Cannot remain in the tunnel
P["tunCVDy2","tunCVDy2"] = 0
# in the next year, you enter into steady-state CVD 
P["tunCVDy2","CVD"] = 1 - P["CVD","Dead"]
# Don't allow entry into steady state until after 2 years
P["Healthy","CVD"] = 0

P[c("Healthy","tunCVDy1","tunCVDy2","CVD","Dead"),c("Healthy","tunCVDy1","tunCVDy2","CVD","Dead")]



