# processData.R
# should be called from SPIKEanalysis.Rmd


########
# LLOQ #
########

# read in and process LLOQ data
lloq <- read_excel(f, sheet = 'LLOQ', na = c('', 'Sample did not dilute down properly-can not use data')) %>%
    rename(acon = `Result (Calculated Con.)`) %>%
  
    # drop out of range samples
    filter(acon != 'Range?' & !is.na(acon)) %>%
    
    # base concentration
    group_by(Assay, Sample_ID, Day, Analyst) %>%
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
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(base_con * min_dil_factor / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    
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

# read in and process ULOQ data
uloq <- read_excel(f, sheet = 'ULOQ') %>%
    rename(acon = `Result (Calculated Con.)`) %>%
    
    # candidate base concentrations (can't just go with the bottom one on this one)
    group_by(Assay, Sample_ID, Day, Analyst, Dil_Factor) %>%
    mutate(base_con = geo_mean(acon, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # pick base concentration closest to 10
    group_by(Assay, Sample_ID, Day, Analyst) %>%
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
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(base_con * base_dil / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    
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

lin <- read_excel(f, sheet = 'LINEARITY', na = c('', 'Range?'))
names(lin)[names(lin) == diff$linearity_acon] <- 'acon'

lin <- group_by(lin, Sample_ID, Day, Analyst, Dil_Factor) %>%
    # candidate base concentrations (can't just go with the bottom one on this one)
    mutate(base_con = geo_mean(acon, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # pick base concentration is closest to 10
    group_by(Assay, Sample_ID, Day, Analyst) %>%
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
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(base_con * base_dil / Dil_Factor), # just need one per group
              
              # geometric mean
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    
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
              `min r` = min(r, na.rm = TRUE),
              `max r` = max(r, na.rm = TRUE)) %>%
    ungroup()


############
# Cutpoint #
############

cutpt <- read_excel(f, sheet = 'CUTPOINT')
names(cutpt)[names(cutpt) == diff$cutpoint_acon] <- 'acon'

# update summary table 1
tables <- summary_table_update(tables, cutpt, 'Cutpoint', 'Calculate Geometric Mean and 95th Percentile', '')
    
cutpt_sum <- group_by(cutpt, Assay, Sample_ID) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              pctl_95 = qtl_limits(xbar, std, qtl = 0.95, log_scale = TRUE)[2]) %>%
    ungroup()

tables$cutpt <- group_by(cutpt, Assay) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              `95th Percentile` = qtl_limits(xbar, std, qtl = 0.95, log_scale = TRUE)[2]) %>%
    ungroup() %>%
    select(-xbar, -std)


###############
# Sensitivity #
###############

sens <- read_excel(f, sheet = 'LLOQ_Challenge', na = c('', 'Plate failed'))
names(sens)[names(sens) == diff$sens_acon] <- 'acon'

sens <- filter(sens, !is.na(acon)) %>%
    mutate(Dil_Factor = 2^(as.numeric(substr(Sample_ID, 7, 7)) - 1)) # calculate dilution factor based on communication with lab folks

# update summary table 1
tables <- summary_table_update(tables, sens, 'Sensitivity (LLOQ Challenge)',
                               'Calculate Geometric Mean, RSD, and Percent Error for each sample concentration level',
                               'Samples at LLOQ and higher must pass these criteria to accept LLOQ: Percent Error ≤ 50%; RSD ≤ 30%')

# set theoretical concentration at median of all observed results
sens_sum <- group_by(sens, Assay) %>%
    mutate(base_con = geo_mean(acon[Sample_ID == 'LLOQ_C1'], na.rm = TRUE)) %>%
    ungroup() %>%
    
    group_by(Assay, Sample_ID) %>%
    summarize(theo_con = unique(base_con / Dil_Factor),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std, log_scale = TRUE),
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
acc <- read_excel(f, sheet = 'ACCURACY') %>%
  filter(!is.na(Assay))
names(acc)[names(acc) == diff$acc_acon] <- 'acon'

# pre-validation data to be used for expected concentration for accuracy
assay_translation <- c(`CoV2 S` = 'CoV2_S', `CoV2 N` = 'CoV2_N', `CoV1 S` = 'CoV1_S',
                       `MERS S` = 'MERS_S', `OC43 S` = 'OC43_S', `229E S` = '229E_S',
                        HKU1    = 'HKU1_S', `NL63 S` = 'NL63_S', `Ragon RBD` = 'R_RBD',
                       `Ragon RBD UK` = 'R_RBD_UK', `Ragon RBD E484K` = 'R_RBD_E484',
                       `Mount Sinai RBD` = 'M_RBD', `Mount Sinai RBD UK` = 'M_RBD_UK',
                       `Mount Sinai RBD SA` = 'M_RBD_SA',
                       `Mount Sinai RBD E484K` = 'M_RBD_E484K')
pre_validation <- read_excel(f2, sheet = 'ACC_Tracking', skip = 4)[,-c(43:46)]

for(i in 2:nrow(pre_validation))
{
  # propagate assay names down (used merged cells in excel)
  if(is.na(pre_validation$Assay[i]))
    pre_validation$Assay[i] <- pre_validation$Assay[i - 1]
}

# add expected concentration for each Assay and Sample_ID
acc <- pivot_longer(pre_validation, cols = 3:42, values_to = 'acon') %>%
  mutate(Sample_ID = gsub('ACC-', '', Samples, fixed = TRUE),
         Assay = assay_translation[Assay]) %>%
  select(-name, -Samples) %>%
  group_by(Assay, Sample_ID) %>%
  summarize(theo_con = geo_mean(acon)) %>%
  ungroup() %>%
  
  right_join(acc, by = c("Assay", "Sample_ID"))
  

# update summary table 1
tables <- summary_table_update(tables, acc, 'Accuracy',
                               'Calculate Geometric Mean and Percent Error', 
                               'Percent Error ≤ 25%')


tables$acc <- acc %>%
    # convert analyst initials to integers in {1,2,...}
    group_by(Assay, Sample_ID) %>%
    mutate(Analyst = as.numeric(as.factor(Analyst))) %>%
    ungroup() %>%
    
    # calculate mean, delta, rsd by Assay,sample_ID/Analyst
    group_by(Assay, Sample_ID, Analyst) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              sd = geo_sd(acon, na.rm = TRUE),
              RSD = rsd(xbar, sd, log_scale = TRUE),
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
prec <- read_excel(f, sheet = diff$prec_sheet)
names(prec)[names(prec) == diff$prec_acon] <- 'acon'
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
           rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup() %>%
    filter(rsd <= 25) %>%
    
    # start with inter-analyst CV
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_inter_analyst = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup() %>%
    
    # inter-day CV
    group_by(Assay, Sample_ID, Day) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_inter_day = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup() %>%
    
    # intra-day CV
    group_by(Assay, Sample_ID) %>%
    mutate(xbar = geo_mean(acon),
           std = geo_sd(acon),
           rsd_intra_day = rsd(xbar, std, log_scale = TRUE)) %>%
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
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup()

prec_within_analyst_by_sample <- group_by(prec, Assay, Sample_ID, Analyst) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup()

prec_by_sample <- group_by(prec, Assay, Sample_ID) %>%
    summarize(xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
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


##############
# Carry-over #
##############

covr <- read_excel(f, sheet = 'CARRYOVER')
names(covr)[names(covr) == diff$covr_acon] <- 'acon'
names(covr)[names(covr) == diff$covr_sample_id] <- 'Sample_ID'
names(covr)[names(covr) == diff$covr_assay] <- 'Assay'
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

stability <- read_excel(f, sheet = 'STABILITY') %>%
  
  # drop 1X, 5X, and 10X from Treatment so we can identify samples
  mutate(test = gsub('\\d+X\\s', '', Treatment))

names(stability)[names(stability) == diff$stab_acon] <- 'acon'


stability_sum <- filter(stability, !is.na(acon)) %>%
  
  mutate(lvl = substr(Sample_ID, 1, 1), # figure out concentration level (low/high)
         trt = Sample_ID,               # the sample id is really the treatment here
         Sample_ID = Treatment) %>%     # when grouping, we really want to group by "Treatment" in the spreadsheet - swap labels for use in get_pct_err_by_analyst()

  # theoretical concentration
  group_by(Assay, lvl) %>%
  mutate(theo_con = median(acon[Sample_ID == ifelse(grepl('Freeze/Thaw', Sample_ID), '1X Freeze/Thaw', 'Control')])) %>%
  ungroup() %>%

  get_pct_err_by_analyst() %>%
  
  filter(! Sample_ID %in% c('Control', '1X Freeze/Thaw'))

tables$stability <- mutate(stability_sum, Treatment = Sample_ID) %>%
  group_by(Assay, Treatment) %>%
  summarize(`Pct Error` = paste0('≤', round(max(delta), digits = 1)))


# summary table information
i <- which(tables$table1$Experiment == 'Stability (Freeze/Thaw, Matrix Components, Sample Integrity)')

tables$table1$`Samples (n)`[i] <- group_by(stability, test) %>% 
  summarize(n = length(unique(Sample_ID))) %>% 
  select(n) %>% 
  range() %>% 
  unique() %>%
  paste(collapse = '-')

tables$table1$Replicates[i] <- group_by(stability, test) %>% 
  summarize(n = length(unique(Replicate_per_sample))) %>% 
  select(n) %>% 
  range() %>% 
  unique() %>%
  paste(collapse = '-')

tables$table1$`Analyst(s)`[i] <- group_by(stability, test) %>% 
  summarize(n = length(unique(Analyst))) %>% 
  select(n) %>% 
  range() %>% 
  unique() %>%
  paste(collapse = '-')

tables$table1$`Days (Runs)`[i] <- group_by(stability, test) %>% 
  summarize(n = length(unique(Day))) %>% 
  select(n) %>% 
  range() %>% 
  unique() %>%
  paste(collapse = '-')

tables$table1$Methodology[i] <- 'Calculate Geometric Mean and Percent Error'
tables$table1$`Acceptance Criteria`[i] <- 'Percent Error ≤ 25%'



##########################
# Lot-to-Lot Comparisons #
##########################

conj <- read_excel(f, sheet = 'Lot-to-Lot Conjugate') %>%
    rename(lot = `Lot Status: Old or New`)
names(conj)[names(conj) == diff$conj_acon] <- 'acon'
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
names(antigen_sum)[names(antigen_sum) == diff$antigen_acon] <- 'acon'
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


#######################
# Single vs Multiplex #
#######################

svm <- read_excel(f, sheet = 'Single vs Multiplex') %>%
  mutate(multiplex = grepl('multiplex', `Result File`))
names(svm)[names(svm) == diff$svm_acon] <- 'acon'

# update summary table 1
tables <- summary_table_update(tables, svm, 'Single vs Multiplex',
                               'Calculate Geometric Mean and Percent Error', 
                               'Percent Error ≤ 25%')

# comparison of the two assays
svm_sum <- svm %>%
  group_by(Assay, Sample_ID, Analyst) %>%
  summarize(single    = geo_mean(acon[!multiplex]),
            multiplex = geo_mean(acon[ multiplex]),
            delta     = pct_err(multiplex, single)) %>%
  ungroup()

tables$svm <- svm_sum %>%
  group_by(Assay) %>%
  summarize(`Total samples` = length(delta),
            `Samples passing` = sum(delta <= 25),
            `% passing` = paste0(round(`Samples passing` / `Total samples`, 3) * 100, '%')) %>%
  ungroup()
