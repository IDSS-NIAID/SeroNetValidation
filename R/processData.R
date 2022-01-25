# processData.R
# should be called from SPIKEanalysis.Rmd


########
# LLOQ #
########

# read in and process LLOQ data
lloq <- read_excel(f, sheet = 'LLOQ', na = c('', 'Sample did not dilute down properly-can not use data')) %>%
    rename(acon = `Result (Calculated Con.)`) %>%
    
    # flag out-of-range values
    mutate(range = acon == 'Range?',
           
           # set acon to 0 where it is out of range on the low side
           # set acon to Inf where it is out of range on the high side
           # use optical density to differentiate
           acon = case_when(range & OD <  1 ~ 0,
                            range & OD >= 1 ~ Inf,
                            TRUE            ~ suppressWarnings(as.numeric(acon))),
           
           # Sample_ID: LLOQ<sample>-<replicate>, where
           #            sample in 1:8
           #            replicate in 1:4
           # replicate is redundant with `Replicate_per_sample`, so we drop it here
           Sample_ID = map_chr(Sample_ID, ~ strsplit(.x, '-')[[1]][1])) %>%
    
    # drop out of range samples
    filter(!range) %>%
    
    # base concentration
    group_by(Sample_ID, Day, Analyst) %>%
    mutate(min_dil_factor = min(Dil_Factor),
           base_con = geo_mean(acon[Dil_Factor == min(Dil_Factor)])) %>%
    ungroup()

# update summary table 1
tables <- summary_table_update(tables, lloq, 'LLOQ',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# set dilution factor
dilution_factor <- 50

# calculate statistics for each group
lloq_sum <- group_by(lloq, Assay, Sample_ID, Day, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
    summarize(theo_con = unique(base_con * min_dil_factor / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon),
              std = geo_sd(acon),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std)) %>%
    
    ungroup()

tables <- get_summary_table1(tables, 
                             test_name = 'lloq',
                             stats = lloq_sum,
                             minmax = min,
                             which.minmax = which.min,
                             pctl = 0.975,
                             dilution_factor = 50) # base dilution factor - this is what was used before, so I'm assuming it is the same. Might be 150?

                          
########
# ULOQ #
########

# read in and process LLOQ data
uloq <- read_excel(f, sheet = 'ULOQ') %>%
    rename(acon = `Result (Calculated Con.)`) %>%
    
    # flag out-of-range values
    mutate(range = acon == 'Range?',
           
           # set acon to 0 where it is out of range on the low side
           # set acon to Inf where it is out of range on the high side
           # use optical density to differentiate
           acon = case_when(range & OD <  1 ~ 0,
                            range & OD >= 1 ~ Inf,
                            TRUE            ~ suppressWarnings(as.numeric(acon))),
           
           # Sample_ID: LLOQ<sample>-<replicate>, where
           #            sample in 1:8
           #            replicate in 1:4
           # replicate is redundant with `Replicate_per_sample`, so we drop it here
           Sample_ID = map_chr(Sample_ID, ~ strsplit(.x, '-')[[1]][1])) %>%
    
    # candidate base concentrations (can't just go with the bottom one on this one
    group_by(Sample_ID, Day, Analyst, Dil_Factor) %>%
    mutate(base_con = geo_mean(acon[!range])) %>%
    ungroup() %>%
    
    # pick base concentration is closest to 10
    group_by(Sample_ID, Day, Analyst) %>%
    mutate(base_dil = Dil_Factor[which.min(abs(10 - base_con))],
           base_con = base_con[which.min(abs(10 - base_con))]) %>%
    ungroup()

tables <- summary_table_update(tables, uloq, 'ULOQ',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# set dilution factor - checking on this?
#dilution_factor <- 36450

# calculate statistics for each group
uloq_sum <- group_by(uloq, Assay, Sample_ID, Day, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
    summarize(theo_con = unique(base_con * base_dil / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon),
              std = geo_sd(acon),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std)) %>%
    
    ungroup()

tables <- get_summary_table1(tables, 
                             test_name = 'uloq',
                             stats = uloq_sum,
                             minmax = max,
                             which.minmax = which.max,
                             pctl = 0.025,
                             dilution_factor = 50) # base dilution factor - this is what was used before, so I'm assuming it is the same. Might be 150?


#############
# Linearity #
#############

lin <- read_excel(f, sheet = 'LINEARITY')
names(lin)[names(lin) == ifelse('linearity_acon' %in% names(diff), diff$linearity_acon, 'Result (Calculated Con.)')] <- 'acon'

    # flag out-of-range values
lin <- mutate(lin, range = acon == 'Range?',
           
           # set acon to 0 where it is out of range on the low side
           # set acon to Inf where it is out of range on the high side
           # use optical density to differentiate
           acon = case_when(range & OD <  1 ~ 0,
                            range & OD >= 1 ~ Inf,
                            TRUE            ~ suppressWarnings(as.numeric(acon))),
           
           # Sample_ID: LLOQ<sample>-<replicate>, where
           #            sample in 1:8
           #            replicate in 1:4
           # replicate is redundant with `Replicate_per_sample`, so we drop it here
           Sample_ID = map_chr(Sample_ID, ~ strsplit(.x, '-')[[1]][1])) %>%
    
    # candidate base concentrations (can't just go with the bottom one on this one
    group_by(Sample_ID, Day, Analyst, Dil_Factor) %>%
    mutate(base_con = geo_mean(acon[!range])) %>%
    ungroup() %>%
    
    # pick base concentration is closest to 10
    group_by(Sample_ID, Day, Analyst) %>%
    mutate(base_dil = Dil_Factor[which.min(abs(10 - base_con))],
           base_con = base_con[which.min(abs(10 - base_con))]) %>%
    ungroup()

# update summary table 1
tables <- summary_table_update(tables, lin, 'Linearity',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# calculate statistics for each group
lin_sum <- group_by(lin, Assay, Sample_ID, Day, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
    summarize(theo_con = unique(base_con * base_dil / Dil_Factor), # just need one per group
              
              # geometric mean
              xbar = geo_mean(acon),
              std = geo_sd(acon),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std)) %>%
    
    ungroup()

# calculate linearity summaries by Assay, Sample_ID
tables$lin <- filter(lin_sum, delta < 50 & rsd < 30) %>%
    group_by(Assay, Sample_ID) %>%
    summarize(n = length(xbar),
              r = ifelse(n < 3, NA, cor.test(theo_con, xbar)$estimate)) %>%
    ungroup() %>%
    
    group_by(Assay) %>%
    summarize(`min Accept` = ifelse(length(n) == filter(tables$table1, Experiment == 'Linearity')$`Samples (n)`,
                                    min(n), 0),
              `max Accept` = max(n),
              `min r` = min(r),
              `max r` = max(r)) %>%
    ungroup()


############
# Cutpoint #
############

cutpt <- read_excel(f, sheet = 'CUTPOINT')
names(cutpt)[names(cutpt) == ifelse('cutpoint_acon' %in% names(diff), diff$cutpoint_acon, 'Result (AU/mL)')] <- 'acon'

# update summary table 1
tables <- summary_table_update(tables, cutpt, 'Cutpoint', 'Calculate Geometric Mean and 95% CI', '')
    
cutpt_sum <- group_by(cutpt, Assay, Sample_ID) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              pctl_95 = exp(qnorm(0.95, log(xbar), log(std)))) %>%
    ungroup()

tables$cutpt <- group_by(cutpt, Assay) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              `95th Percentile` = exp(qnorm(0.95, log(xbar), log(std)))) %>%
    ungroup() %>%
    select(-xbar, -std)


###############
# Sensitivity #
###############

sens <- read_excel(f, sheet = 'LLOQ_Challenge', na = c('', 'Plate failed'))
names(sens)[names(sens) == ifelse('sens_acon' %in% names(diff), diff$sens_acon, 'Result (AU/mL)')] <- 'acon'

sens <- filter(sens, !is.na(acon)) %>%
    mutate(Dil_Factor = 2^(as.numeric(substr(Sample_ID, 7, 7)) - 1)) # calculate dilution factor based on communication with lab folks

# update summary table 1
tables <- summary_table_update(tables, sens, 'Sensitivity (LLOQ Challenge)',
                               'Calculate Geometric Mean, RSD, and Percent Error for each sample concentration level',
                               'Samples at LLOQ and higher must pass these criteria to accept LLOQ: Percent Error ≤ 50%; RSD ≤ 30%')

# set theoretical concentration at median of all observed results
sens_sum <- group_by(sens, Assay) %>%
    mutate(base_con = geo_mean(acon[Sample_ID == 'LLOQ_C1'])) %>%
    ungroup() %>%
    
    group_by(Assay, Sample_ID) %>%
    summarize(theo_con = unique(base_con / Dil_Factor),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std),
              delta = pct_err(xbar, theo_con)) %>%
    ungroup()

tables$sens <- filter(sens_sum, delta <= 50 & rsd <= 30) %>%
    group_by(Assay) %>%
    summarize(`Pct Error` = delta[which.min(theo_con)],
              RSD = rsd[which.min(theo_con)],
              `Concentration (AU/mL)` = xbar[which.min(theo_con)]) %>%
    ungroup()

##########################
# Accuracy and Precision #
##########################

#### Accuracy ####
acc <- read_excel(f, sheet = 'ACCURACY')
names(acc)[names(acc) == ifelse('acc_acon' %in% names(diff), diff$acc_acon, 'Result (AU/mL)')] <- 'acon'

acc <- filter(acc, !is.na(acon)) %>%
    
    # capture dilution series
    mutate(dil_series = {strsplit(Sample_ID, '_', fixed = TRUE) %>%
                         map_chr(~ gsub(.x[2], pattern = 'L', replacement = '', fixed = TRUE)) %>%
                         as.numeric()}) %>%
    
    # set theoretical concentration at what is observed at ACC_01, dilution series is 1/3 each time
    group_by(Assay, Analyst) %>%
    mutate(theo_con = mean(acon[dil_series == 1], na.rm = TRUE)) %>%
    ungroup() %>%
    
    mutate(theo_con = theo_con / (3^(dil_series - 1)))

# update summary table 1
tables <- summary_table_update(tables, acc, 'Accuracy',
                               'Calculate Geometric Mean and Percent Error', 
                               'Percent Error ≤ 25%')


tables$acc <- 
    # filter such that delta ≤ 25%
    mutate(acc, 
           delta = pct_err(acon, theo_con)) %>%
    #filter(delta <= 25) %>%
    
    # convert analyst initials to integers in {1,2,...}
    group_by(Assay, Sample_ID) %>%
    mutate(Analyst = as.numeric(as.factor(Analyst))) %>%
    ungroup() %>%
    
    # calculate mean, delta, rsd by Assay,sample_ID/Analyst
    group_by(Assay, Sample_ID, Analyst) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              sd = geo_sd(acon, na.rm = TRUE),
              RSD = rsd(xbar, sd),
              `Pct Error` = pct_err(xbar, unique(theo_con))) %>%
    ungroup() %>%
    rename(`Geometric Mean (AU/mL)` = xbar) %>%
    select(-sd)

# other tables
acc_by_day_analyst <- group_by(acc, Assay, Sample_ID, Analyst, Day) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              delta = pct_err(xbar, unique(theo_con))) %>%
    ungroup()

acc_by_analyst <- group_by(acc, Assay, Sample_ID, Analyst) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              delta = pct_err(xbar, unique(theo_con))) %>%
    ungroup()

acc_by_sample <- group_by(acc, Assay, Sample_ID) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              delta = pct_err(xbar, unique(theo_con))) %>%
    ungroup()

# put things together for comparison of within / between group pct error
acc_by_day <- full_join(acc_by_day_analyst,
                        acc_by_analyst,
                        by = c("Assay", "Sample_ID", "Analyst", "xbar", "delta")) %>%
    arrange(Assay, Sample_ID, Analyst, Day) %>%
    mutate(Day = ifelse(is.na(Day), 'Between Days', Day))

acc_by_analyst <- full_join(acc_by_analyst,
                            acc_by_sample,
                            by = c("Assay", "Sample_ID", "xbar", "delta")) %>%
    arrange(Assay, Sample_ID, Analyst) %>%
    mutate(Analyst = ifelse(is.na(Analyst), 'Between Analysts', Analyst))


#### Precision ####
prec <- read_excel(f, sheet = ifelse('prec_sheet' %in% names(diff), diff$prec_sheet, 'PRECISION'))
names(prec)[names(prec) == ifelse('prec_acon' %in% names(diff), diff$prec_acon, 'Result (AU/mL)')] <- 'acon'
prec <- filter(prec, !is.na(acon))

# update summary table 1
tables <- summary_table_update(tables, prec, 'Precision',
                               'Calculate Geometric Mean and RSD for Intra-plate, Inter-plate, and Inter-Analyst', 
                               'RSD ≤ 25% for Intra-plate, Inter-plate, and Inter- Analyst')

tables$prec <- 
    # filter such that delta ≤ 25%
    group_by(prec, Assay, Sample_ID, Analyst, Day) %>%
    mutate(xbar = geo_mean(acon, na.rm = TRUE),
           std = geo_sd(acon, na.rm = TRUE),
           rsd = rsd(xbar, std)) %>%
    ungroup() %>%
    filter(rsd <= 25) %>%
    
    # start with inter-analyst CV
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_inter_analyst = rsd(xbar, std)) %>%
    ungroup() %>%
    
    # inter-day CV
    group_by(Assay, Sample_ID, Day) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_inter_day = rsd(xbar, std)) %>%
    ungroup() %>%
    
    # intra-day CV
    group_by(Assay, Sample_ID) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_intra_day = rsd(xbar, std)) %>%
    ungroup() %>%
    
    # summarize
    group_by(Assay) %>%
    summarize(`Intra-Day CV` = mean(rsd_intra_day),
              `Inter-Day CV` = mean(rsd_inter_day),
              `Inter-Analyst CV` = mean(rsd_inter_analyst)) %>%
    ungroup()

# Other tables
prec_within_day_by_analyst_sample <- group_by(prec, Assay, Sample_ID, Analyst, Day) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std)) %>%
    ungroup()

prec_within_analyst_by_sample <- group_by(prec, Assay, Sample_ID, Analyst) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std)) %>%
    ungroup()

prec_by_sample <- group_by(prec, Assay, Sample_ID) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std)) %>%
    ungroup()

# put things together for comparison of within / between group RSD
prec_by_day <- full_join(prec_within_day_by_analyst_sample,
                         prec_within_analyst_by_sample,
                         by = c("Assay", "Sample_ID", "Analyst", "xbar", "std", "rsd")) %>%

    mutate(Day = ifelse(is.na(Day), 'Between Days', Day)) %>%
    arrange(Assay, Sample_ID, Analyst)

prec_by_analyst <- full_join(prec_within_analyst_by_sample,
                             prec_by_sample,
                             by = c("Assay", "Sample_ID", "xbar", "std", "rsd")) %>%
    arrange(Assay, Sample_ID, Analyst) %>%
    mutate(Analyst = ifelse(is.na(Analyst), 'Between Analysts', Analyst))


###############
# Specificity #
###############

spec <- read_excel(f, sheet = 'Specificity') %>%
    rename(acon = `Result (AU/mL)`) %>%
    filter(!is.na(acon)) %>%
    
    mutate(group = map_chr(Sample_ID, ~ strsplit(.x, 'PC')[[1]][2]))
    
# update summary table 1
tables <- summary_table_update(tables, spec, 'Specificity',
                               'Calculate Geometric Mean and RSD', 
                               '1) ≥90% of sample concentration must be inhibited with ... type specific VLPs. 2) Non-spike type specific inhibition must be ≤ 25%.')

    # expected concentration
spec_sum <- group_by(spec, Assay, group) %>%
    mutate(theo_con = geo_mean(acon[`Spiked Antigen` == 'Unspiked'])) %>%
    ungroup() %>%
    
    # don't need these anymore
    filter(`Spiked Antigen` != 'Unspiked') %>%

    # percent error
    group_by(Assay, Sample_ID, theo_con) %>%
    summarize(xbar = geo_mean(acon),
              delta = pct_err(xbar, unique(theo_con))) %>%
    ungroup()
    

tables$spec <- group_by(spec_sum, Assay, Sample_ID) %>%
    summarize(`Pct Error` = mean(delta))

##############
# Carry-over #
##############

covr <- read_excel(f, sheet = 'CARRYOVER')
names(covr)[names(covr) == ifelse('covr_acon' %in% names(diff), diff$covr_acon, 'Result (AU/mL)')] <- 'acon'
covr <- filter(covr, !is.na(acon))

# update summary table 1
tables <- summary_table_update(tables, covr, 'Carry-over',
                               'Calculate Geometric Mean and RSD', 
                               'RSD ≤ 25% for Low and High sample. The Neg_Serum sample must be ≤ LLOQ.')

covr_sum <- get_rsd_by_analyst(covr)

tables$covr <- group_by(covr_sum, Assay, Sample_ID) %>%
    summarize(`Pct Error` = mean(rsd)) %>%
    ungroup()


#############
# Stability #
#############

# row for stability updates
i <- which(tables$table1$Experiment == "Stability (Freeze/Thaw, Matrix Components, Sample Integrity)")

##### Freeze Thaw Data #####

freez <- read_excel(f, sheet = 'STABILITY_FREEZE-THAW')
names(freez)[names(freez) == ifelse('freez_acon' %in% names(diff), diff$freez_acon, 'Result (AU/mL)')] <- 'acon'
freez <- filter(freez, !is.na(acon)) %>%
    
    # figure out which are high concentration (the rest are low)
    mutate(high = substr(Sample_ID, 1, 4) == 'High')

# update summary table 1 (freeze/thaw)
tab1 <- summary_table_update(tables, freez, 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)',
                             'Calculate Geometric Mean and Percent Error', 
                             'Percent Error ≤ 25%')$table1[i,]

    # theoretical concentration
freez_sum <- group_by(freez, Assay, high) %>%
    mutate(theo_con = median(acon)) %>%
    ungroup() %>%
    
    filter(is.na(Treatment)) %>%
    
    get_pct_err_by_analyst()

                
tables$freez <- tibble(Assay = unique(freez_sum$Assay),  # remove reference  
                       `Freeze/Thaw Pct Error` = get.stability.max(filter(freez_sum, !grepl('1x', Sample_ID))))


##### Matrix Components #####

# lipemia
lip <- read_excel(f, sheet = 'STABILITY_LIPEMIA', na = c('', 'Flagged for repeat'))
names(lip)[names(lip) == ifelse('lip_acon' %in% names(diff), diff$lip_acon, 'Result (AU/mL)')] <- 'acon'
lip <- filter(lip, !is.na(acon)) %>%
    
    # figure out concentration level
    mutate(lvl = substr(Sample_ID, 1, 1))

# update summary table 1 (lipemia)
tab1 <- summary_table_update(tables, lip, 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)',
                             'Calculate Geometric Mean and Percent Error', 
                             'Percent Error ≤ 25%')$table1[i,] %>%
    bind_rows(tab1)

    # theoretical concentration
lip_sum <- group_by(lip, Assay, lvl) %>%
    mutate(theo_con = median(acon[is.na(Treatment)])) %>%
    ungroup() %>%
    
    filter(!is.na(Treatment)) %>%
    
    get_pct_err_by_analyst()

# Billirubin
bili <- read_excel(f, sheet = 'STABILITY_BILIRUBIN', na = c('', 'Flagged for repeat'))
names(bili)[names(bili) == ifelse('bili_acon' %in% names(diff), diff$bili_acon, 'Result (AU/mL)')] <- 'acon'
bili <- filter(bili, !is.na(acon)) %>%
    
    # figure out concentration level
    mutate(lvl = substr(Sample_ID, 1, 1))

# update summary table 1 (lipemia)
tab1 <- summary_table_update(tables, bili, 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)',
                             'Calculate Geometric Mean and Percent Error', 
                             'Percent Error ≤ 25%')$table1[i,] %>%
    bind_rows(tab1)

    # theoretical concentration
bili_sum <- group_by(bili, Assay, lvl) %>%
    mutate(theo_con = median(acon[is.na(Treatment)])) %>%
    ungroup() %>%
    
    filter(!is.na(Treatment)) %>%
    
    get_pct_err_by_analyst()

# hemolysis
hemo <- read_excel(f, sheet = 'STABILITY_HEMOLYSIS', na = c('', 'Flagged for repeat'))
names(hemo)[names(hemo) == ifelse('hemo_acon' %in% names(diff), diff$hemo_acon, 'Result (AU/mL)')] <- 'acon'
hemo <- filter(hemo, !is.na(acon)) %>%
    
    # figure out concentration level
    mutate(lvl = substr(Sample_ID, 1, 1))

# update summary table 1 (hemolysis)
tab1 <- summary_table_update(tables, hemo, 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)',
                             'Calculate Geometric Mean and Percent Error', 
                             'Percent Error ≤ 25%')$table1[i,] %>%
    bind_rows(tab1)

    # theoretical concentration
hemo_sum <- group_by(hemo, Assay, lvl) %>%
    mutate(theo_con = median(acon[is.na(Treatment)])) %>%
    ungroup() %>%
    
    filter(!is.na(Treatment)) %>%
    
    get_pct_err_by_analyst()

# additional updates for table1

tables$matrx <- tibble(Assay = unique(lip_sum$Assay),
                       `Hemoglobin Pct Error` = get.stability.max(hemo_sum),
                       `Lipid Pct Error` = get.stability.max(lip_sum),
                       `Bilirubin Percent Error` = get.stability.max(bili_sum)) %>%
    filter(`Hemoglobin Pct Error` != '≤-Inf' | `Lipid Pct Error` != '≤-Inf' | `Bilirubin Percent Error` != '≤-Inf')


##### Heat Treatment #####

heat <- read_excel(f, sheet = 'STABILITY_HEAT', na = c('', 'No heat', 'No Heat')) # sometimes they leave treatment blank, other times they label it "No h/Heat" - easier just to have it be missing
names(heat)[names(heat) == ifelse('heat_acon' %in% names(diff), diff$heat_acon, 'Result (AU/mL)')] <- 'acon'
heat <- filter(heat, !is.na(acon)) %>%
    
    # fix non-standard labels (compared to other stability tabs)
    mutate(Sample_ID = paste0(gsub('Heat ', '', Sample_ID, fixed = TRUE),
                              ifelse(is.na(Treatment), '', ' Heat')),
           
           # figure out concentration level
           lvl = substr(Sample_ID, 1, 1))

# update summary table 1 (heat)
tab1 <- summary_table_update(tables, heat, 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)',
                             'Calculate Geometric Mean and Percent Error', 
                             'Percent Error ≤ 25%')$table1[i,] %>%
    bind_rows(tab1)

    
    # theoretical concentration
heat_sum <- group_by(heat, Assay, lvl) %>%
    mutate(theo_con = median(acon[is.na(Treatment)])) %>%
    ungroup() %>%
    
    filter(!is.na(Treatment)) %>%
    
    get_pct_err_by_analyst()

tables$heat <- tibble(Assay = unique(heat_sum$Assay),       # remove reference  
                      `Heat Pct Error` = get.stability.max(filter(heat_sum, !grepl('LOW', Sample_ID))))

# final updates of Table 1
tables$table1$`Samples (n)` %<>% as.character
tables$table1$Replicates %<>% as.character
tables$table1[i,] <- tibble(Experiment = unique(tab1$Experiment),
                            `Samples (n)` = paste(range(tab1$`Samples (n)`), collapse = '-'),
                            `Replicates` = paste(range(tab1$`Replicates`), collapse = '-'),
                            `Analyst(s)` = unique(tab1$`Analyst(s)`),
                            `Days (Runs)` = unique(tab1$`Days (Runs)`),
                            Methodology = unique(tab1$Methodology),
                            `Acceptance Criteria` = unique(tab1$`Acceptance Criteria`))

##########################
# Lot-to-Lot Comparisons #
##########################

conj <- read_excel(f, sheet = 'Lot-to-Lot Conjugate') %>%
    rename(lot = `Lot Status: Old or New`)
names(conj)[names(conj) == ifelse('conj_acon' %in% names(diff), diff$conj_acon, 'Result (AU/mL)')] <- 'acon'
conj <- filter(conj, !is.na(acon) & !is.na(Assay)) %>%
    mutate(lot = tolower(lot))

# update summary table 1
tables <- summary_table_update(tables, conj, 'Stability (Critical Reagent Lot Change)',
                               'Calculate Geometric Mean and Percent Error', 
                               'Percent Error ≤ 25%')

    # expected concentration
conj_sum <- group_by(conj, Assay, Sample_ID) %>%
    mutate(theo_con = geo_mean(acon[lot == 'old'])) %>%
    ungroup() %>%
    
    filter(lot == 'new') %>%
    
    group_by(Assay, Sample_ID, theo_con) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              delta = pct_err(xbar, theo_con)) %>%
    ungroup()


antigen_sum <- read_excel(f, sheet = 'Lot-to-Lot Antigen', na = c('', 'Undefined (>10,000)', '>70000')) %>%
    rename(lot = `Lot Status: Old or New`)
names(antigen_sum)[names(antigen_sum) == ifelse('antigen_acon' %in% names(diff), diff$antigen_acon, 'Result (AU/mL)')] <- 'acon'
antigen_sum <- filter(antigen_sum, !is.na(acon) & !is.na(Assay)) %>%
    mutate(lot = tolower(lot)) %>%
    
    # expected concentration
    group_by(Assay, Sample_ID) %>%
    mutate(theo_con = geo_mean(acon[lot == 'old'])) %>%
    ungroup() %>%
    
    filter(lot == 'new') %>%
    
    group_by(Assay, Sample_ID, theo_con) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              delta = pct_err(xbar, theo_con)) %>%
    ungroup()

tables$l2l <- group_by(conj_sum, Assay) %>%
    summarize(`Conjugate Passing` = sum(delta <= 25),
              `Conjugate Total` = length(delta),
              `Conjugate w/ Error ≤ 25%` = paste0(round(`Conjugate Passing` / `Conjugate Total`, 3) * 100, "%")) %>%
    ungroup()

tables$l2l <- group_by(antigen_sum, Assay) %>%
    summarize(`Antigen Passing` = sum(delta <= 25),
              `Antigen Total` = length(delta),
              `Antigen w/ Error ≤ 25%` = paste0(round(`Antigen Passing` / `Antigen Total`, 3) * 100, "%")) %>%
    ungroup() %>%
    right_join(tables$l2l, by = 'Assay')
