library(tidyverse)
library(here)
library(glue)
hmd <- 
    tibble::tribble(
    ~code, ~country,
    "AUS",                              "Australia",
    "AUT",                               "Austria",
    "BLR",                               "Belarus",
    "BEL",                               "Belgium",
    "BGR",                              "Bulgaria",
    "CAN",                                "Canada",
    "CHL",                                 "Chile",
    "CZE",                        "Czech Republic",
    "DNK",                               "Denmark",
    "EST",                               "Estonia",
    "FIN",                               "Finland",
    "FRATNP",             "France (total population)",
    "FRACNP",          "France (civilian population)",
    "DEUTNP",            "Germany (total population)",
    "DEUTE",                        "Germany (east)",
    "DEUTW",                        "Germany (west)",
    "GRC",                                "Greece",
    "HUN",                               "Hungary",
    "ISL",                               "Iceland",
    "IRL",                               "Ireland",
    "ISR",                                "Israel",
    "ITA",                                 "Italy",
    "JPN",                                 "Japan",
    "LVA",                                "Latvia",
    "LTU",                             "Lithuania",
    "LUX",                            "Luxembourg",
    "NLD",                           "Netherlands",
    "NZL_NP",        "New Zealand (total population)",
    "NZL_MA",        "New Zealand (Maori population)",
    "NZL_NM",    "New Zealand (non-Maori population)",
    "NOR",                                "Norway",
    "POL",                                "Poland",
    "PRT",                              "Portugal",
    "RUS",                                "Russia",
    "SVK",                              "Slovakia",
    "SVN",                              "Slovenia",
    "ESP",                                 "Spain",
    "SWE",                                "Sweden",
    "CHE",                           "Switzerland",
    "TWN",                                "Taiwan",
    "GBR_NP",                        "United Kingdom",
    "GBRTENW",    "England & Wales (total population)",
    "GBRCENW", "England & Wales (civilian population)",
    "GBR_SCO",                              "Scotland",
    "GBR_NIR",                      "Northern Ireland",
    "USA",                                "U.S.A.",
    "UKR",                               "Ukraine"
) 


countries <- c("AUS"   ,   "AUT"    ,  "BLR"   ,  "BEL"   ,  "BGR"     ,"CAN"  ,   "CHL"  ,    "HRV"  ,    "CZE"  ,   "DNK"     , "EST"    , "FIN"  ,   "FRATNP" , "DEUTNP" , "GRC"  ,  "HKG"  ,    "HUN" ,  
 "ISL"  ,   "IRL"    ,   "ISR"  ,   "ITA"  ,   "JPN"    , "LVA" ,    "LTU" ,     "LUX" ,     "NLD" ,    "NZL_NP" ,  "NOR"   ,  "POL" ,    "PRT"   ,  "KOR"   ,  "RUS"  ,  "SVK"  ,    "SVN" ,  
 "ESP"  ,   "SWE"    ,   "CHE"  ,   "TWN"  ,   "GBR_NP" , "USA" ,    "UKR"     )



countries %>% map(~({
    if (!file.exists(here(glue("_hmd_lifetables/life-tables/{.x}_{yy}.rds")))) {
        cat(glue("Reading {.x}\n"))
        tmp <- demography::hmd.mx(.x,username = hmd_u, password = hmd_p, paste0(.x))
        yy = max(tmp$year)
        lt <- demography::lifetable(tmp,series = "total", years = yy) 
        saveRDS(lt,here(glue("_hmd_lifetables/life-tables/{.x}_{yy}.rds")))
    }
}))

source("~/auth-aws.r")
.x  = "AUS"
mort_year = 2019

