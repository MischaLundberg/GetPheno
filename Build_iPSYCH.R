####################################
####    Comparison to iPsych    #### 
####################################
#### This script was created by ####
#### Johanne Østerby Sørensen   ####
#### johanneoesterby@gmail.com  ####
####----------------------------####
####       and modified by      ####
####       Mischa Lundberg      ####
#### mischa.lundberg@gmail.com  ####
####################################

# Clear environment
rm (list = ls())

# Load necessary libraries
library (data.table)

# set directory, when using singularity
#setwd ( "/work dir" )

# set directory when using modules
setwd ( "" )

# Read in datatable
#scores <- fread("/dpibp/shared/aSchorkLab/FGRS/scores/PAFGRS_11psych_240501.csv", header = TRUE, sep = "," )
PSYK_LPR <- fread("/dpibp/shared/johost/depression age_project/getpheno/output/250218 getpheno mdd aao_adhd_SUD_with_LPR", header = TRUE)
#ipsyk_def <- fread("/dpibp/shared/johost/depression age_project/getpheno/output/250205_ipsych definitions", header = TRUE)
comorb <- fread("/dpibp/shared/johost/depression_age_project/getpheno/output/250206_getphone mdd_aao_results", header = TRUE)

# make more simple datatable ipsyk definitions
# c <- ipsyk_def[, c(1,8,11, 14,17,20,23,26,29,32,35,38,41,44,47,50,51,52,53,54, 57,63) ]
# c[, statd := as.Date(statd)]
# c[, birthdate := as.Date(birthdate)]

# make more simple datatable
c < - comorb[, c (1,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,69,71,72,70,75,81,82)]
c[, statd := as.Date(statd)]
c[, birthdate := as.Date(birthdate)]

# restricted cohort
# born in DK
# Danish municipalities                     10 - 900
# Danish courts                             1101 - 1199
# Danish state offices, part 1              1301 - 1315
# Danish state offices, part 2              1317 - 1325
# Undisclosed place in Denmark              2401 - 2599
# Danish churches                           4601 - 4688
# Danish church districts                   6001 - 6903
# Danish parishes                           7001 - 9348
# Denmark (country)                         5100    
# Partially undisclosed place in Denmark    4998

# keep rows where fkode is in this list
c1 <- c[fkode %in% c(10:900,
                    1101:1199,
                    1301:1315,
                    1317:1325,
                    2401:2599,
                    4601:4688,
                    6001:6903,
                    7001:9348,
                    5100,
                    4998)]

# born to mother with a cpr number?
dk_mom <- c1[!is.na(pnr_m),]
c1 <- dk_mom

# selected birth window
c2 <- c1[birthdate > "1981-04-30" & birthdate < "2009-01-01",]

# alive and residing in DK on their 1 year birthday
c3 <- c2[c2$statd >= (c2$birthdate + 365) |
        c2$stat == "1" | 
        c2$stat == "3", ]

# exclude twins
# Remove all rows where mother-birthdate combination appears more than once
# fromLast = TRUE in duplicated() checks for duplicates starting from the end of the dataset instead of the beginning.
c4 <- c3[!duplicated(c3, by = c("pnr_m", "birthdate")) &
        !duplicated(c3, by = c("pnr_m", "birthdate"), fromLast = TRUE)]

# subset to only depression cases
c5 <- c4[c4$MDD == "Case",]

# get age at onset
c5[, MDD_earliest_date := as.Date(MDD_earliest_date)]
c5$age_at_MDD <- as.numeric(c5$MDD_earliest date - c5$birthdate) / 365.25

# exclude diagnosis <10y to make it comparable to iPsych
с6 <- c5[c5$age_at_MDD >= 10,]

# year of diagnosis before 2016
c6$year_of_first_Dx_MDD <- as.numeric(format(as.Date(c6$MDD_earliest_date), "%Y"))
c7 <- c6[year_of_first_Dx_MDD < 2016,]

# make AAO groups
c7$child_onset <- ifelse(c7$age_at_MDD < 14, 1, 0)
c7$adolescent_onset <- ifelse(c7$age_at_MDD >= 14 &
                              c7$age_at_MDD < 18, 1, 0)
c7$adult_onset <- ifelse(c7$age_at_MDD >= 18, 1, 0)

# add onset group column
c7$onset group <- factor(ifelse(c7$child_onset == 1, "child",
                         ifelse(c7$adolescent_onset == 1, "adolescent",
                         "adult")))

## import GetPheno with PSYK + LPR to get cases of ADHD and SUD
psyklpr <- PSYK_LPR[, c(1,9,13,23,27)] 
comorb_dt <- merge(c7, psyklpr, by = "pnr", all.x = TRUE)
names(comorb_dt) <- gsub("-y$", "", names(comorb_dt))

write(comorb_dt, "250219_MDD_AAO_comorbidities data_ipsych_restrictions.csv") 

##########################################################################
## Check no. of other cases to compare to iPsych with ipsych-definitions
# get year of diagnosis
earliest date <- c("skiz_earliest _date",
                    "skizospek _earliest_date",
                    "bipol_earliest_date",
                    "affek_earliest_date",
                    "autism_earliest_date",
                    "adhd_earliest_date",
                    "anorek_earliest_date")


# Loop through and create new year columns
for(i in earliest_date) {
    new_col_name <- paste(i, "_year")
    c4[[new_col_name]] <- as.numeric(format(as.Date(c4[[i]]), "%Y"))
}

skiz <- c4[skiz_earliest_date_year < 2016,] 
skizospek <- c4[skizospek_earliest_date_year < 2016,] 
bipol <- c4[bipol_earliest_date_year < 2016,] 
affek <- c4faffek_earliest date_year < 2016,]
autism <- c4[autism_earliest_date_year < 2016,] 
adhd <- c4|adhd_earliest_date_year < 2016,] 
anorek <- c4[anorek_earliest_date_year < 2016,]

phenotypes = "skiz\tF20
skizospek\tF20,F21,F22,F23,F24,F25,F28,F29
bipol\tF30,F31
affek\F30,F31,F32,F33,F34,F35,F36,F37,F38,F39
autism\tF84.0,F84.1,F84.5,F84.8,F84.9
adhd\tF90.0
anorek\tF50.0,F50.1"
