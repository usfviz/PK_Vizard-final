library(dplyr)

#raw data files downloaded from: 
#https://www.fda.gov/drugs/guidancecomplianceregulatoryinformation/surveillance/adversedrugeffects/ucm082193.htm
drug <- read.csv('DRUG16Q4.txt', sep='$')
reac <- read.csv('REAC16Q4.txt', sep = '$')
demo <- read.csv('DEMO16Q4.txt', sep = '$')

drug <- drug[drug$prod_ai != "", c('caseid', 'drug_seq', 'role_cod', 'drugname', 'prod_ai',
                                   'route')]
reac <- reac[, c('caseid', 'pt')]
reac <- reac[!duplicated(reac), ]

demo <- demo[, c('caseid', 'fda_dt', 'age', 'age_cod', 'age_grp', 'sex', 'e_sub', 
                 'occp_cod', 'reporter_country', 'occr_country')]

standardize_age <- function(age, df) {
    age[df$age_cod == 'HR'] <- df[df$age_cod == 'HR', 'age'] / (365.0 * 24)
    age[df$age_cod == 'DY'] <- df[df$age_cod == 'DY', 'age'] / (365.0)
    age[df$age_cod == 'WK'] <- df[df$age_cod == 'WK', 'age'] / (52.0)
    age[df$age_cod == 'MON'] <- df[df$age_cod == 'MON', 'age'] / (12.0)
    age[df$age_cod == 'DEC'] <- df[(df$age_cod == 'DEC') & (df$age < 12), 'age'] * 10
    age[age < 0] <- 0
    return(age)
}

standard_age <- standardize_age(demo$age, demo)
demo$age <- standard_age

drugs_count <- drug[!duplicated(drug), ] %>%
    group_by(prod_ai) %>%
    summarize(count = n()) %>% 
    arrange(-count)


reacs_count <- reac %>% 
    group_by(pt) %>%  
    summarize(count = n()) %>% 
    arrange(-count)

write.table(demo, 'demo_1May.txt', sep="$")
write.table(drug, 'drug_1May.txt', sep="$")
write.table(reac, 'reac_1May.txt', sep="$")
write.table(drugs_count, 'drugCount_1May.txt', sep="$")
write.table(reacs_count, 'reacCount_1May.txt', sep="$")
