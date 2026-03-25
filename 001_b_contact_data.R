
# read dta file
library(haven)
library(tidyverse)

# 2018 household ----------------------------------------------------------

# read family dta file
data_family <- read_dta("../data/CLDS/household2018.dta")

# select columns: contains birthyear
data_family_clean <- data_family |> 
     filter(PROV2018 == 43) |> 
     select(FID2018, main, contains("birthyear"), contains("out")) |> 
     mutate(FID2018 = as.character(FID2018),
            # set default value for main is 1
            main = if_else(is.na(main), 1, main)) |> 
     pivot_longer(cols = !c(FID2018, main), names_to = "variable", values_to = "value") |> 
     # drop na values
     filter(!is.na(value)) |> 
     # separate variable into two columns by _
     separate(variable, into = c("variable", "number"), sep = "_") |> 
     pivot_wider(names_from = variable, values_from = value) |> 
     # drop out family members
     filter(is.na(hhmemberoutreason)) |> 
     # select columns
     select(FID2018, main, number, birthyear) |> 
     # calculate age
     mutate(age = 2018 - birthyear,
            number = as.numeric(number))

# find object age
data_family_main <- data_family_clean[data_family_clean$number == data_family_clean$main,c('FID2018', 'age')]

# add main_age to data_family_clean
data_family_2018 <- data_family_clean |> 
     filter(number != main) |> 
     left_join(data_family_main, by = 'FID2018') |> 
     rename(age = age.x,
            age_main = age.y) |> 
     select(age, age_main, FID = FID2018)

remove(data_family_main, data_family)

# 2016 household ----------------------------------------------------------

# read family dta file
data_family <- read_dta("../data/CLDS/household2016.dta")

# select columns: contains birthyear
data_family_clean <- data_family |> 
     filter(PROV2016 == 43) |> 
     select(FID2016, contains("Fmbirthy"), contains("Fmoutcause")) |> 
     mutate(FID2016 = as.character(FID2016),
            # set default value for main is 1
            main = 1) |> 
     pivot_longer(cols = !c(FID2016, main), names_to = "variable", values_to = "value") |> 
     # drop na values
     filter(!is.na(value)) |> 
     # separate variable into two columns by _
     separate(variable, into = c("variable", "number"), sep = "_") |> 
     pivot_wider(names_from = variable, values_from = value) |> 
     # drop out family members
     filter(is.na(Fmoutcause)) |> 
     # select columns
     select(FID2016, main, number, birthyear = Fmbirthy) |> 
     # calculate age
     mutate(age = 2016 - birthyear,
            number = as.numeric(number))

# find object age
data_family_main <- data_family_clean[data_family_clean$number == data_family_clean$main,c('FID2016', 'age')]

# add main_age to data_family_clean
data_family_2016 <- data_family_clean |> 
     filter(number != main) |> 
     left_join(data_family_main, by = 'FID2016') |> 
     rename(age = age.x,
            age_main = age.y) |> 
     filter(!is.na(age_main)) |> 
     select(age, age_main, FID = FID2016)

remove(data_family_main, data_family_clean, data_family)

# 2014 household ----------------------------------------------------------

# read family dta file
data_family <- read_dta("../data/CLDS/CLDS2014-household-160520-STATA.dta")

# select columns: contains birthyear
data_family_clean <- data_family |> 
     filter(PROVINCE == 430000) |> 
     select(FID2014, fmnum, contains('fmage'), contains('fmoutcause')) |> 
     # drop columns contains oth
     select(-contains('oth')) |>
     mutate(FID2014 = as.character(FID2014),
            # set default value for main is 1
            main = '1') |> 
     pivot_longer(cols = !c(FID2014, main, fmnum), names_to = "variable", values_to = "value") |>
     # drop na values
     filter(!is.na(value)) |>
     # separate variable into two columns by _
     separate(variable, into = c("variable", "number"), sep = "_") |>
     pivot_wider(names_from = variable, values_from = value) |>
     # drop out family members
     filter(is.na(fmoutcause)) |>
     # select columns
     select(FID2014, main, number, age = fmage)

# find object age
data_family_main <- data_family_clean[data_family_clean$number == data_family_clean$main,c('FID2014', 'age')]

# add main_age to data_family_clean
data_family_2014 <- data_family_clean |> 
     filter(number != main) |> 
     left_join(data_family_main, by = 'FID2014') |> 
     rename(age = age.x,
            age_main = age.y) |> 
     select(age, age_main, FID = FID2014)

rm(data_family_main, data_family_clean, data_family)

# 2012 household ----------------------------------------------------------

# read family dta file
data_family <- read_dta("../data/CLDS/household2012.dta")

# select columns: contains birthyear
data_family_clean <- data_family |> 
     filter(PROVINCE == 430000) |> 
     # select columns end iwth *_3_1
     select(PROVINCE, FID, matches(".*_3_1$")) |> 
     mutate(FID = as.character(FID),
            # set default value for main is 1
            main = 'F1_1_101_3_1') |> 
     pivot_longer(cols = !c(PROVINCE, FID, main), names_to = "number", values_to = "value") |> 
     # drop na values
     filter(!is.na(value) & value > 1900 & value < 2014) |>
     # select columns
     select(FID, main, number, birthyear = value) |> 
     # calculate age
     mutate(age = 2012 - birthyear)

# find object age
data_family_main <- data_family_clean[data_family_clean$number == data_family_clean$main,c('FID', 'age')]

# add main_age to data_family_clean
data_family_2012 <- data_family_clean |> 
     filter(number != main) |> 
     left_join(data_family_main, by = 'FID') |> 
     rename(age = age.x,
            age_main = age.y) |> 
     select(age, age_main, FID)

rm(data_family_main, data_family_clean, data_family)

# combine data ----------------------------------------------------------

data_family <- bind_rows(data_family_2018, data_family_2016, data_family_2014, data_family_2012)

write.csv(data_family,
          "../data/contact_CLDS_family.csv",
          row.names = FALSE)
