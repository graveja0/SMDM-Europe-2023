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
select <- dplyr::select
options("scipen" = 100, "digits" = 5)

gbd_lt = read_rds(here("_gbd_lifetables/life-tables-2019.rds"))


# Original HIV Model

m_P <- matrix(c(0.721, 0.202, 0.067, 0.010,
                0.000, 0.581, 0.407, 0.012,
                0.000, 0.000, 0.750, 0.250,
                0.000, 0.000, 0.000, 1.000), 
              nrow=4, byrow=TRUE,
              dimnames=list(c("A", "B", "C", "D"),
                            c("A", "B", "C", "D")))
m_P  

# Eigenvector decomposition to get the continuous generator

V  <- eigen(m_P)$vectors
iV <- solve(V)
Ap <- iV %*% m_P %*% V
lAp <- diag(log(diag(Ap)), nrow(Ap), ncol(Ap))
R  <- V %*% lAp %*% iV
R 
R[abs(R) < 1e-6 ] <- 0
rownames(R) <- c("A", "B", "C", "D")
colnames(R) <- c("A", "B", "C", "D")
R

location = 1
run_sim <- function(location) {
    
    mort_coef <- 
        gbd_lt %>% 
        filter(location_id==location) %>% 
        pull(hp_mort) %>% 
        pluck(1) %>% 
        as.vector()
    
    params_hiv <- 
        list(
            R = R,
            dmca =	 1701 , # Direct medical costs associated with state A
            dmcb =	 1774 , # Direct medical costs associated with state B
            dmcc =	 6948 , # Direct medical costs associated with state C
            ccca =	 1055 , # Community care costs associated with state A
            cccb =	 1278 , # Community care costs associated with state B
            cccc =	 2059 , # Community care costs associated with state C
            cAZT =	 2278 , # Zidovudine drug cost
            cLam =	 2087 , # Lamivudine drug cost
            RR	  = 0.509, # Treatment effect (RR)
            cDR	 =  0.06,   # Annual discount rate - costs (%)
            oDR	 =  0       # Annual discount rate - benefits (%)
        )
    
    starting_age = 40
    
    markov_rate_matrices <- function(t)
    {
        lapply(t, function(tt) {
            current_age <- tt + starting_age - 1
            
            r_death <- HP2(current_age, mort_coef)$hx
            
            x <- params_hiv$R
            x["A","D"] <- r_death
            
            diag(x) <- 0
            diag(x) <- -rowSums(x)
            x
        })
    }
    
    transition_prob <- function(t,h=1,combo = FALSE)
    {
        # Embed the matrices into the timesteps
        tp <- lapply(markov_rate_matrices(t), function(m) expm(m*h))
        
        if (combo) {
            m_RR =            matrix(c(1, params_hiv$RR, params_hiv$RR, params_hiv$RR,
                                       1, 1, params_hiv$RR, params_hiv$RR,
                                       1, 1, 1, params_hiv$RR,
                                       1, 1, 1, 1), 
                                     nrow=4, byrow=TRUE,
                                     dimnames=list(c("A", "B", "C", "D")))
            tp[[1]] = tp[[1]] * m_RR
            diag(tp[[1]]) <- rep(0,4)
            diag(tp[[1]]) <- 1 - rowSums(tp[[1]])
            tp[[2]] = tp[[2]] * m_RR
            diag(tp[[2]]) <- rep(0,4)
            diag(tp[[2]]) <- 1 - rowSums(tp[[2]])
            
        }

        tp
        
    }
    
    Y <- t(c(A = 1, B = 0, C = 0, D = 0))
    
    res_mono <- do.call(rbind, lapply(transition_prob(1:20), function(tp) {
        Y <<- Y %*% tp
    }))
    
    res_mono <- rbind(
        c(A = 1, B = 0, C = 0, D = 0), 
        res_mono
    )
    
    
    Y <- t(c(A = 1, B = 0, C = 0, D = 0))
    res_combo <- do.call(rbind, lapply(transition_prob(1:20, combo = TRUE), function(tp) {
        Y <<- Y %*% tp
    }))
    res_combo <- rbind(
        c(A = 1, B = 0, C = 0, D = 0), 
        res_combo
    )

    
    payoff_health  = matrix(c(1,1,1,0),byrow=TRUE,ncol=1,dimnames=list(c("A","B","C","D"),c("qaly")))
    
    qaly_mono <- as.matrix(res_mono[-1,]) %*% payoff_health
    qaly_disc_mono <- qaly_mono / (1 + params_hiv$oDR)^c(1:(nrow(res_mono)-1))
    
    qaly_combo <- as.matrix(res_combo[-1,]) %*% payoff_health
    qaly_disc_combo <- qaly_combo / (1 + params_hiv$oDR)^c(1:(nrow(res_combo)-1))
    
    
    payoff_cost_mono <- with(params_hiv,matrix(c(dmca + ccca + cAZT,
                                                 dmcb + cccb + cAZT,
                                                 dmcc + cccc + cAZT,
                                                 0)))
    cost_mono <- as.matrix(res_mono[-1,]) %*% payoff_cost_mono %>% round(.,0)
    cost_disc_mono <- cost_mono / (1 + params_hiv$cDR)^c(1:(nrow(res_mono)-1))
    
    
    payoff_cost_combo <- with(params_hiv,matrix(c(dmca + ccca + cLam + cAZT,
                                                 dmcb + cccb + cLam + cAZT,
                                                 dmcc + cccc + cLam + cAZT,
                                                 0)))
    
    payoff_cost_combo2 <- with(params_hiv,matrix(c(dmca + ccca + cAZT,
                                                  dmcb + cccb + cAZT,
                                                  dmcc + cccc +  cAZT,
                                                  0)))
    
    
    cost_combo1 <- as.matrix(res_combo[2:3,]) %*% payoff_cost_combo 
    cost_combo2 <- as.matrix(res_combo[-(1:3),]) %*% payoff_cost_combo2 
    cost_combo <- rbind(cost_combo1,cost_combo2)
    
    cost_disc_combo <- cost_combo / (1 + params_hiv$cDR)^c(1:(nrow(res_combo)-1))
    
    
    data.frame(location = location, qaly_mono = sum(qaly_disc_mono) , cost_mono = sum(cost_disc_mono),
               qaly_combo = sum(qaly_disc_combo) , cost_combo = sum(cost_disc_combo))
}

res <- 
    gbd_lt %>% 
    mutate(res = map(location_id,~(run_sim(.x)))) %>% 
    unnest(cols = c(res)) %>% 
    mutate(inc_cost = cost_combo - cost_mono,
           inc_qaly = qaly_combo - qaly_mono) %>% 
    mutate(icer = inc_cost / inc_qaly)

res %>% 
    ggplot(aes(x = icer)) + geom_histogram() +
    hrbrthemes::theme_ipsum() +
    labs(x = "ICER", y = "Count")


###########################
# Model Mortality in GBD
###########################
# The below code only needs to be run once to patch together the various life tables. 
    
# cbd_life_tables <- c(
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_28_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_5_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_6_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_7_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_8_0.zip",    
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_9_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_10_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_11_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_12_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_13_0.zip",    
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_14_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_15_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_16_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_17_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_18_0.zip",    
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_19_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_20_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_30_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_31_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_32_0.zip", 
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_33_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_44_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_45_0.zip",
#     "https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_LIFE_TABLES_1950_2019_ID_148_0.zip"
# )

# cbd_life_tables %>% map(~({
#     cat(glue::glue("File: {.x}\n"))
#     tmp <- tempfile()
#     tmpdir <- tempdir()
#     download.file(.x, tmp)
#     nn <- gsub("https://ghdx.healthdata.org/sites/default/files/record-attached-files/","",.x) %>% gsub(".ZIP|.zip","",.)
#     unzip(tmp)
#     ff <- list.files(here(),pattern = ".CSV|.csv")
#     files <- ff[grep("WSHOCK",ff)] %>% map(~(read.csv(.x) %>% as_tibble() %>% mutate(file = .x))) %>% pluck(1)
#     write_rds(files,file = glue::glue("_gbd_lifetables/{nn}"))
#     file.remove(ff)   
# }))

# radix = 100000

# gbd_lt_ <-
#     list.files(here("_gbd_lifetables")) %>%
#     map_df(~({
#         read_rds(here(glue("_gbd_lifetables/{.x}")))
#     }))
# 
# age_group_lut <- c("28" = 0, "5" = 5, "6" = 10, "7" = 15, "8" = 20,
#                    "9" = 25, "10" = 30, "11" = 35, "12"= 40, "13" = 45,
#                    "14" = 50, "15" = 55, "16" = 60, "17"=65, "18"=70,
#                    "19" = 75, "20" = 80, "30"= 85,"31" = 90,"32" = 95, "33" = 100,
#                    "44"=105,"45" = 110,"148" = 111)
# 
# gbd_lt <-
#     gbd_lt_ %>%
#     filter(year_id==2019) %>%
#     group_by(location_name,location_id) %>%
#     nest()  %>%
#     mutate(life_table = map(data,~({
#         .x %>% select(sex_name,age_group_id,age_group_name,measure_name,val) %>%
#             mutate(age = age_group_lut[paste0(age_group_id)]) %>%
#             spread(measure_name,val) %>%
#             janitor::clean_names() %>%
#             rename(ex = life_expectancy,
#                    qx = probability_of_death) %>%
#             arrange(sex_name,age) %>%
#             select(sex_name,age,ex,qx) %>%
#             gather(measure,value,-sex_name,-age) %>%
#             mutate(measure = paste0(measure,gsub("_both","",paste0("_",sex_name)))) %>%
#             select(-sex_name)  %>%
#             spread(measure,value) %>% 
#             mutate(p = 1 - qx,
#                    p_male = 1 - qx_male,
#                    p_female = 1 - qx_female) %>% 
#             mutate(lx = radix * cumprod(c(1,p[-nrow(.)])),
#                    lx_male = radix * cumprod(c(1,p_male[-nrow(.)])),
#                    lx_female = radix * cumprod(c(1,p_female[-nrow(.)])))  
#         
#     })))
# gbd_lt <- 
#     gbd_lt %>% 
#     mutate(hp_mort = map(life_table,~({
#         ages     <- .x$age
#         deaths   <- .x$lx * .x$qx
#         exposure <- .x$lx
#         
#         mort_fit <- MortalityLaw(
#             x  = ages,
#             Dx  = deaths,   # vector with death counts
#             Ex  = exposure, # vector containing exposures
#             law = "HP2",
#             opt.method = "LF2")
#         coef(mort_fit) %>% data.frame() %>% t() 
#     })))


# 
# gbd_lt %>%
#    select(-data) %>% 
#     write_rds(here("_gbd_lifetables/life-tables-2019.rds"))


