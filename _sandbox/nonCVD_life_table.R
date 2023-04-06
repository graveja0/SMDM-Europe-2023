library(tidyverse)
library(readxl)

# get the US life tables from CDC(2020)
# https://www.cdc.gov/nchs/data/nvsr/nvsr71/nvsr71-01.pdf

flt2020_raw = read_xlsx("life table/flt2020.xlsx", skip = 2, col_names = TRUE)
colnames(flt2020_raw)[1] = "age"
flt2020 = flt2020_raw %>%
  drop_na() %>%
  rowwise %>%
  mutate(Age = ifelse(str_detect(age, '–'), 
                      as.numeric(str_sub(age,1,unlist(gregexpr('–', age))[1]-1)), 
                      as.numeric(str_sub(age,1,unlist(gregexpr(' ', age))[1]-1))),
         rate = -log(1-qx)/1)

mlt2020_raw = read_xlsx("life table/flt2020.xlsx", skip = 2, col_names = TRUE)
colnames(mlt2020_raw)[1] = "age"
mlt2020 = mlt2020_raw %>%
  drop_na() %>%
  rowwise %>%
  mutate(Age = ifelse(str_detect(age, '–'), 
                      as.numeric(str_sub(age,1,unlist(gregexpr('–', age))[1]-1)), 
                      as.numeric(str_sub(age,1,unlist(gregexpr(' ', age))[1]-1))),
         rate = -log(1-qx)/1)

# get the CVD-cause death percents from GBD(2019)
# https://vizhub.healthdata.org/gbd-results/

CVD_death_raw = read.csv("life table/GBD_2019_CVD_death_percent.csv")
CVD_death = CVD_death_raw %>%
  select(sex,age,val,upper,lower)
CVD_death_f = CVD_death %>% filter(sex == "Female") %>% select(-sex)
CVD_death_m = CVD_death %>% filter(sex == "Male") %>% select(-sex)

# change the age group into a constant number and replicate the left parts

CVD_death_f1 = CVD_death_f %>% filter(str_detect(age, "<")) %>%
  mutate(start = 0, 
         reptimes = as.numeric(str_sub(age,
                            unlist(gregexpr('<', age))[1]+1,
                            unlist(gregexpr(' ', age))[1]-1)),
         end = start + reptimes - 1)

CVD_death_f2 = CVD_death_f %>% 
  filter(str_detect(age,"-") & !str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,str_length(age))),
         reptimes = end - start + 1)

CVD_death_f3 = CVD_death_f %>% filter(str_detect(age,"-") & str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,unlist(gregexpr(' ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_f4 = CVD_death_f %>% filter(str_detect(age, "[+] ")) %>%
  rowwise() %>%
  mutate(end = max(flt2020$Age),
         start = as.numeric(str_sub(age, 1, unlist(gregexpr('[+] ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_female = rbind(
  as.data.frame(lapply(CVD_death_f1, rep, CVD_death_f1$reptimes)) %>%
  group_by(age) %>%
  mutate(Age = start - 1 + 1:n()) %>%
  ungroup(),
  as.data.frame(lapply(CVD_death_f2, rep, CVD_death_f2$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_f3, rep, CVD_death_f3$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_f4, rep, CVD_death_f4$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup()) %>%
  arrange(Age) %>%
  select(Age, val, upper, lower)

maxage = min(max(flt2020$Age), max(CVD_death_female$Age))
minage_cause = min(CVD_death_female$Age)

## dx
flt2020_nonCVD = rbind(ifelse(minage_cause>1, flt2020[1:minage_cause, c("Age","dx")], NA),
                       data.frame(Age = minage_cause:maxage, 
                                  dx = flt2020[(minage_cause+1):(maxage+1), "dx"] * (1-CVD_death_female[, "val"])),
                       ifelse(maxage < max(flt2020$Age), flt2020[(maxage+2):nrow(flt2020),c("Age","dx")], NA)) %>%
  drop_na()

flt2020_CVD = data.frame(Age = minage_cause:maxage,
                         dx = flt2020[(minage_cause+1):(maxage+1), "dx"] * CVD_death_female[, "val"])

sum(flt2020_CVD$dx) + sum(flt2020_nonCVD$dx)

# directly using the percent and qx
for (i in 1:nrow(flt2020)) {
  if (i %in% (minage_cause+1):(maxage+1)) {
    flt2020_nonCVD[i, "qx"] = flt2020[i, "qx"] * (1-CVD_death_female[i - minage_cause, "val"])
  } else {
    flt2020_nonCVD[i, "qx"] = flt2020[i, "qx"]
  }
}

write.csv(flt2020_nonCVD, "life table/flt2020-CVD-deleted.csv", row.names = FALSE)

## male

CVD_death_m1 = CVD_death_m %>% filter(str_detect(age, "<")) %>%
  mutate(start = 0, 
         reptimes = as.numeric(str_sub(age,
                                       unlist(gregexpr('<', age))[1]+1,
                                       unlist(gregexpr(' ', age))[1]-1)),
         end = start + reptimes - 1)

CVD_death_m2 = CVD_death_m %>% 
  filter(str_detect(age,"-") & !str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,str_length(age))),
         reptimes = end - start + 1)

CVD_death_m3 = CVD_death_m %>% filter(str_detect(age,"-") & str_detect(age, " ")) %>%
  rowwise %>%
  mutate(start = as.numeric(str_sub(age,1,unlist(gregexpr('-', age))[1]-1)),
         end = as.numeric(str_sub(age,unlist(gregexpr('-', age))[1]+1,unlist(gregexpr(' ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_m4 = CVD_death_m %>% filter(str_detect(age, "[+] ")) %>%
  rowwise() %>%
  mutate(end = max(mlt2020$Age),
         start = as.numeric(str_sub(age, 1, unlist(gregexpr('[+] ', age))[1]-1)),
         reptimes = end - start + 1)

CVD_death_male = rbind(
  as.data.frame(lapply(CVD_death_m1, rep, CVD_death_m1$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_m2, rep, CVD_death_m2$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_m3, rep, CVD_death_m3$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup(),
  as.data.frame(lapply(CVD_death_m4, rep, CVD_death_m4$reptimes)) %>%
    group_by(age) %>%
    mutate(Age = start - 1 + 1:n()) %>%
    ungroup()) %>%
  arrange(Age) %>%
  select(Age, val, upper, lower)

maxage = min(max(mlt2020$Age), max(CVD_death_male$Age))
minage_cause = min(CVD_death_male$Age)

## dx
mlt2020_nonCVD = rbind(ifelse(minage_cause>1, mlt2020[1:minage_cause, c("Age","dx")], NA),
                         data.frame(Age = minage_cause:maxage, 
                                    dx = mlt2020[(minage_cause+1):(maxage+1), "dx"] * (1-CVD_death_male[, "val"])),
                         ifelse(maxage < max(mlt2020$Age), mlt2020[(maxage+2):nrow(mlt2020),c("Age","dx")], NA)) %>%
  drop_na()

mlt2020_CVD = data.frame(Age = minage_cause:maxage,
                           dx = mlt2020[(minage_cause+1):(maxage+1), "dx"] * CVD_death_male[, "val"])

sum(mlt2020_CVD$dx) + sum(mlt2020_nonCVD$dx)

# directly using the percent and qx
for (i in 1:nrow(mlt2020)) {
  if (i %in% (minage_cause+1):(maxage+1)) {
    mlt2020_nonCVD[i, "qx"] = mlt2020[i, "qx"] * (1-CVD_death_male[i - minage_cause, "val"])
  } else {
    mlt2020_nonCVD[i, "qx"] = mlt2020[i, "qx"]
  }
}

write.csv(mlt2020_nonCVD, "life table/mlt2020-CVD-deleted.csv", row.names = FALSE)

# validation
ggplot(flt2020_nonCVD, aes(x = Age, y = dx)) +
  geom_point(color = "blue")+
  geom_point(data = flt2020_CVD, aes(x = Age, y = dx), color = "red") +
  ggtitle("Female Number dying between ages x and x + 1 (blue: non-CVD cause, red: CVD cause)")

ggplot(flt2020_nonCVD, aes(x = Age, y = cumsum(dx))) +
  geom_point(color = "blue")+
  geom_point(data = flt2020_CVD, aes(x = Age, y = cumsum(dx)), color = "red") +
  ggtitle("Female Number dying between ages 0 and x + 1 (blue: non-CVD cause, red: CVD cause)")


