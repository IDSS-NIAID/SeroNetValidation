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
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(min_dil_factor = min(Dil_Factor),
           base_con = geo_mean(acon[Dil_Factor == min(Dil_Factor)])) %>%
    ungroup()

# update summary table 1
tables <- summary_table_update(tables, lloq, 'LLOQ',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# calculate statistics for each group
lloq_sum <- group_by(lloq, Assay, Sample_ID, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(base_con * min_dil_factor / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              n = sum(!is.na(acon)),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    
    ungroup()

# estimation of threshold for use in plots below (see `get_summary_table1` for more details)
lloq_thresh <- lloq_sum %>%
  filter(delta <= 50 & rsd <= 30) %>%
  group_by(Sample_ID, Assay, Analyst) %>%
  summarize(statxbar = min(xbar, na.rm = TRUE),
            statsd = std[which.min(xbar)],
            dilFact = Dil_Factor[which.min(xbar)],
            n = n[which.min(xbar)]) %>%
  ungroup() %>%
  group_by(Assay) %>%
  summarize(lloq = map_dbl(1, ~ exp(metamean(n, log(statxbar), log(statsd))$TE.fixed))) %>%
  ungroup()

# create lloq table
tables <- get_summary_table1(tables, 
                             test_name = 'lloq',
                             stats = lloq_sum,
                             minmax = min,
                             which.minmax = which.min,
                             pctl = 0.975)

                          
########
# ULOQ #
########

# read in and process ULOQ data
uloq <- read_excel(f, sheet = 'ULOQ') %>%
    rename(acon = `Result (Calculated Con.)`) %>%
    
    # candidate base concentrations (can't just go with the bottom one on this one)
    group_by(Assay, Sample_ID, Analyst, Dil_Factor) %>%
    mutate(base_con = geo_mean(acon, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # pick base concentration closest to 10
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(base_dil = Dil_Factor[which.min(abs(10 - base_con))],
           base_con = base_con[which.min(abs(10 - base_con))]) %>%
    ungroup()

tables <- summary_table_update(tables, uloq, 'ULOQ',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# calculate statistics for each group
uloq_sum <- group_by(uloq, Assay, Sample_ID, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(base_con * base_dil / Dil_Factor), # just need one per group
              
              # mean
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              n = sum(!is.na(acon)),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    
    ungroup()

# estimation of threshold for use in plots below (see `get_summary_table1` for more details)
uloq_thresh <- uloq_sum %>%
  filter(delta <= 50 & rsd <= 30) %>%
  group_by(Sample_ID, Assay, Analyst) %>%
  summarize(statxbar = max(xbar, na.rm = TRUE),
            statsd = std[which.max(xbar)],
            dilFact = Dil_Factor[which.max(xbar)],
            n = n[which.max(xbar)]) %>%
  ungroup() %>%
  group_by(Assay) %>%
  summarize(uloq = map_dbl(1, ~ exp(metamean(n, log(statxbar), log(statsd))$TE.fixed))) %>%
  ungroup()

# create uloq table
tables <- get_summary_table1(tables, 
                             test_name = 'uloq',
                             stats = uloq_sum,
                             minmax = max,
                             which.minmax = which.max,
                             pctl = 0.025)


#############
# Linearity #
#############

inconsistent_names <- c('CoV1  S'            = 'CoV1 S',
                        'CoV1 S'             = 'CoV1 S',
                        
                        'Mt.Sinai RBD'       = 'M RBD',
                        'M RBD'              = 'M RBD',
                        
                        'Mt.Sinai RBD E484K' = 'M RBD E484K',
                        'M RBD E484K'        = 'M RBD E484K',
                        
                        'Mt.Sinai RBD SA'    = 'M RBD SA',
                        'M RBD SA'           = 'M RBD SA',
                        
                        'Mt.Sinai RBD UK'    = 'M RBD UK',
                        'M RBD UK'           = 'M RBD UK',
                        
                        '229E S'             = '229E S',
                        
                        'CoV2 N'             = 'CoV2 N',
                        
                        'CoV2 S'             = 'CoV2 S',
                        
                        'HKU1 S'             = 'HKU1 S',
                        
                        'MERS S'             = 'MERS S',
                        
                        'OC43 S'             = 'OC43 S',
                        
                        'NL63 S'             = 'NL63 S',
                        
                        'Ragon RBD'          = 'Ragon RBD',
                        
                        'Ragon RBD UK'       = 'Ragon RBD UK',
                        
                        'RagonRBD E484K'     = 'RagonRBD E484K',
                        'Ragon RBD E484K'    = 'RagonRBD E484K')

full_assay_name <- c('CoV1 S'         = 'CoV-1 S',
                     'M RBD'          = 'RBD WT (M)',
                     'M RBD E484K'    = 'RBD mutant E484K (M)',
                     'M RBD SA'       = 'RBD triple mutant (M)',
                     'M RBD UK'       = 'RBD mutant N501Y (M)',
                     '229E S'         = '229E CoV S',
                     'CoV2 N'         = 'CoV-2 N',
                     'CoV2 S'         = 'CoV-2 S',
                     'HKU1 S'         = 'HKU1 CoV S',
                     'MERS S'         = 'MERS CoV S',
                     'OC43 S'         = 'OC43 CoV S',
                     'NL63 S'         = 'NL63 CoV S',
                     'Ragon RBD'      = 'RBD WT (R)',
                     'Ragon RBD UK'   = 'RBD mutant N501Y (R)',
                     'RagonRBD E484K' = 'RBD mutant E484K (R)')
                       

lin <- read_excel(f, sheet = 'LINEARITY', na = c('', 'Range?')) %>%
  mutate(Assay = inconsistent_names[Assay])
names(lin)[names(lin) == diff$linearity_acon] <- 'acon'

lin <- group_by(lin, Sample_ID, Analyst, Dil_Factor) %>%
    # candidate base concentrations (can't just go with the bottom one on this one)
    mutate(base_con = geo_mean(acon, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # pick base concentration is closest to 10
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(base_dil = Dil_Factor[which.min(abs(10 - base_con))],
           base_con = base_con[which.min(abs(10 - base_con))],
           theo_con = base_con * base_dil / Dil_Factor,
           
           # added for debugging purposes
           delta_indiv = pct_err(acon, theo_con),
           drop = delta_indiv >= 50) %>%
    ungroup()

# update summary table 1
tables <- summary_table_update(tables, lin, 'Linearity',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               'Percent Error ≤ 50%; RSD ≤ 30%')

# calculate statistics for each group
lin_sum <- group_by(lin, Assay, Sample_ID, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
      # base_con * min_dil_factor /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor                   # divide by dilution factor to calculate expected concentration for each dilution
    summarize(theo_con = unique(theo_con), # just need one per group
              
              # sample statistics
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              min_acon = min(acon, na.rm = TRUE),
              max_acon = max(acon, na.rm = TRUE),
              
              # delta (% error)    
              delta = pct_err(xbar, theo_con),
              
              # Relative Standard Deviation
              rsd = rsd(xbar, std, log_scale = TRUE),
              
              # grab the number of observed replicates and the number with missing data
              nobs = sum(!is.na(acon)),
              nmissing = sum(is.na(acon)),
              
              # id used for filtering below
              id = unique(paste(Assay, Sample_ID, Analyst, Dil_Factor))) %>%
    
    ungroup()

##### plots #####

# start by filtering out parts of lin that we aren't using
lin <- lin %>%
  mutate(id = paste(Assay, Sample_ID, Analyst, Dil_Factor),
         keep = id %in% filter(lin_sum, delta < 50 & rsd < 30)$id,
         lacon = log(acon),
         ltheo_con = log(theo_con)) %>%
  
  # get upper and lower bounds
  group_by(Assay) %>%
  mutate(lloq = map_dbl(Assay, ~ filter(lloq_thresh, Assay == .x)$lloq),
         uloq = map_dbl(Assay, ~ filter(uloq_thresh, Assay == .x)$uloq)) %>%
  ungroup() %>%

  # filter out observations we don't want to keep
  filter(keep | (lloq <= acon & acon <= uloq))

# now calculate slopes
lin_assay_sum <- lin %>%
  group_by(Assay) %>%
  summarize(model_good = map(unique(Assay), ~ 
                               {
                                 tmp <- filter(lin, Assay == .x & keep)
                                 if(length(tmp$lacon) > 3)
                                 {
                                   return(lme(lacon ~ ltheo_con,
                                              random = ~ 1 | Sample_ID,
                                              data = tmp,
                                              na.action = na.omit))
                                  }else{
                                    return(NULL)
                                  }
            #                    }),
            # model_bad = map(unique(Assay), ~ 
            #                    {
            #                      tmp <- filter(lin, Assay == .x & !keep)
            #                      if(length(tmp$lacon) > 3)
            #                      {
            #                        return(lme(lacon ~ ltheo_con,
            #                                   random = ~ 1 | Sample_ID,
            #                                   data = tmp,
            #                                   na.action = na.omit))
            #                      }else{
            #                        return(NULL)
            #                      }
                               }),
            
            # number of replicates included in the calculation
            n_good = sum(keep)) %>%
  ungroup() %>%
  
  mutate(intercept = map_dbl(model_good, ~ 
                               {
                                 if(is.null(.x)) 
                                   return(NA)
                                 .x$coefficients$fixed['(Intercept)']
                               }),
         slope = map_dbl(model_good, ~ 
                           {
                             if(is.null(.x)) 
                               return(NA)
                             .x$coefficients$fixed['ltheo_con']
                           }))

ranges <- lin_sum %>%
  filter(delta < 50 & rsd < 30 & theo_con > 0 & !is.na(xbar)) %>%
  group_by(Assay) %>%
  summarize(amin = min(c(theo_con, xbar)),
            amax = max(c(theo_con, xbar))) %>%
  ungroup() %>%
  summarize(amin = min(amin),
            amax = max(amax),
            lmin = log10(amin),
            lmax = log10(amax))

axis_labs <- tibble(breaks = c(seq(from =    0.001, to =    0.009, length = 9),
                               seq(from =    0.01 , to =    0.09 , length = 9),
                               seq(from =    0.1  , to =    0.9  , length = 9),
                               seq(from =    1    , to =    9    , length = 9),
                               seq(from =   10    , to =   90    , length = 9),
                               seq(from =  100    , to =  900    , length = 9),
                                          1000                                ),
                    labels = c(   '0.001', rep('', 8),
                                  '0.01' , rep('', 8),
                                  '0.1'  , rep('', 8),
                                  '1'    , rep('', 8),
                                 '10'    , rep('', 8),
                                '100'    , rep('', 8),
                               '1000')) %>%
  filter(breaks > ranges$amin & breaks < ranges$amax)

xrange <- range(log10(lin_sum$theo_con))
yrange <- range(log10(lin_sum$xbar), na.rm = TRUE) 

figures$linearity_OvE_concentration <- 
  map(unique(lin_sum$Assay), ~ filter(lin_sum, delta < 50 & rsd < 30 &
                                      !is.na(xbar) & theo_con > 0 &
                                      Assay == .x) %>%
        ggplot(aes(theo_con, xbar)) +
        
        # plot points on log10 scale
        scale_x_log10(limits = c(ranges$amin, ranges$amax), 
                      breaks = axis_labs$breaks, labels = axis_labs$labels) +
        scale_y_log10(limits = c(ranges$amin, ranges$amax), 
                      breaks = axis_labs$breaks, labels = axis_labs$labels) +
        
        # add trend line
        geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
        geom_abline(slope = 1, intercept = 0) +

        # layer data on top
        geom_point(color = 'blue') +
        
        # add spread for each mean
        geom_errorbar(aes(ymin = min_acon, ymax = max_acon), width = 0, color = 'blue') +
        
        # labels
        annotate('text', x = 0, y = Inf, 
                 label = paste(' slope =', round(filter(lin_assay_sum, Assay == .x)$slope, 2)),
                 hjust = 0, vjust = 1) +
        ylab('Measured Concentration (AU/mL)') +
        xlab('Expected Concentration (AU/mL)') +
        ggtitle(full_assay_name[.x]))

names(figures$linearity_OvE_concentration) <- unique(lin_sum$Assay)

figures$linearity_OvE_concentration_bad <- 
  map(unique(lin_sum$Assay), ~ filter(lin_sum, 
                                      Assay == .x & 
                                      ((delta < 50 & rsd < 30) |
                                       (theo_con > filter(lloq_thresh, Assay == .x)$lloq &
                                        theo_con < filter(uloq_thresh, Assay == .x)$uloq))) %>%
        ggplot(aes(theo_con, xbar, color = delta < 50 & rsd < 30)) +
        
        # plot points on log10 scale
        scale_x_log10() +
        scale_y_log10(labels = scales::comma) +
        
        # add trend line
        geom_smooth(method = 'lm', se = TRUE) +
        geom_abline(slope = 1, intercept = 0) +
        
        # add cutoff lines
        # geom_vline(xintercept = filter(lloq_thresh, Assay == .x)$lloq, linetype = 2) +
        # geom_vline(xintercept = filter(uloq_thresh, Assay == .x)$uloq, linetype = 2) +

        # layer data on top
        geom_point() +
                
        # labels
        ylab('Titer (AU/mL)') +
        xlab('Theoretical Concentration') +
        ggtitle(.x))

names(figures$linearity_OvE_concentration_bad) <- unique(lin_sum$Assay)

#### save pdf figures ####
for(i in 1:length(figures$linearity_OvE_concentration))
{
  save_plot(paste0(root, '/figs/', names(figures$linearity_OvE_concentration)[i], '.pdf'),
            figures$linearity_OvE_concentration[[i]],
            base_asp = 1)
}

##### final table #####
# calculate linearity summaries by Assay, Sample_ID
tables$lin <- lin_sum %>%
  mutate(pass = delta < 50 & rsd < 30) %>%
  group_by(Assay, Sample_ID) %>%
  summarize(n = length(xbar[pass]),
            r = ifelse(n < 3, NA, cor.test(theo_con[pass], xbar[pass])$estimate)) %>%
  ungroup() %>%
  
  group_by(Assay) %>%
  summarize(`min Accept` = ifelse(length(n) == filter(tables$table1, Experiment == 'Linearity')$`Samples (n)`,
                                  min(n), 0),
            `max Accept` = max(n),
            `min r` = min(r, na.rm = TRUE),
            `max r` = max(r, na.rm = TRUE)) %>%
  ungroup() %>%
  right_join(dplyr::select(lin_assay_sum, Assay, slope), by = "Assay")


############
# Cutpoint #
############

cutpt <- read_excel(f, sheet = 'CUTPOINT')
names(cutpt)[names(cutpt) == diff$cutpoint_acon] <- 'acon'

# update summary table 1
tables <- summary_table_update(tables, cutpt, 'Cutpoint', 'Calculate Geometric Mean and 95% CI', '')
    
cutpt_sum <- group_by(cutpt, Assay, Sample_ID) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              upper_95_CI = qtl_limit(xbar, std, n = n, qtl = 0.975, log_scale = TRUE)) %>%
    ungroup()

tables$cutpt <- group_by(cutpt, Assay) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              `Upper 95% Confidence Bound` = qtl_limit(xbar, std, n = n, qtl = 0.95, log_scale = TRUE)) %>%
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
              `CV%` = rsd[which.min(theo_con)],
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
              `CV%` = rsd(xbar, sd, log_scale = TRUE),
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
prec <- filter(prec, !is.na(acon) &
                 !(Assay == 'CoV2 N' & Interpretation == 'Negative'))

# update summary table 1
tables <- summary_table_update(tables, prec, 'Precision',
                               'Calculate Geometric Mean and RSD for Intra-plate, Inter-plate, and Inter-Analyst', 
                               'RSD ≤ 25% for Intra-plate, Inter-plate, and Inter- Analyst')

tmp <- prec %>%
    # filter such that delta ≤ 25%
    group_by(Assay, Sample_ID, Analyst, Day) %>%
    mutate(xbar = geo_mean(acon, na.rm = TRUE),
           std = geo_sd(acon, na.rm = TRUE),
           rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup() %>%
    filter(rsd <= 25) %>%
    
    # start with inter-analyst CV
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(inter_analyst_xbar = geo_mean(acon),
           inter_analyst_std = geo_sd(acon),
           inter_analyst_n = sum(!is.na(acon))) %>%
    ungroup() %>%
    
    # inter-day CV
    group_by(Assay, Sample_ID, Day) %>%
    mutate(inter_day_xbar = geo_mean(acon),
           inter_day_std = geo_sd(acon),
           inter_day_n = sum(!is.na(acon))) %>%
    ungroup() %>%
    
    # intra-day CV
    group_by(Assay, Sample_ID) %>%
    mutate(intra_day_xbar = geo_mean(acon),
           intra_day_std = geo_sd(acon),
           intra_day_n = sum(!is.na(acon))) %>%
    ungroup()

# `tmp` has all sorts of duplicate rows for each different grouping, so we do it this way
tables$prec <- group_by(tmp, Assay) %>% 
  summarize(`Intra-Day CV` = NA,
            `Inter-Day CV` = NA,
            `Inter-Analyst CV` = NA) %>% 
  ungroup()
    
for(i in 1:nrow(tables$prec))
{
  # Intra-Day CV
    # means between High, Mid, and Low are pretty different, but standard deviations should be similar
    # Since we only care about standard deviation when log_scale = TRUE, we are OK here
    # we want the standard deviation -> use the weights in `obj` to estimate
  tables$prec[['Intra-Day CV']][i] <- tmp %>%
    filter(Assay == tables$prec$Assay[i]) %>%                            # keep only this assay
    select(Sample_ID, intra_day_xbar, intra_day_std, intra_day_n) %>%    # keep only intra-day statistics
    unique() %>%                                                         # keep only unique rows
    with(sum(log(intra_day_std)^2 * intra_day_n) / sum(intra_day_n)) %>% # compute weighted sum of variances, weighting by n
    sqrt() %>%                                                           # convert to standard deviation
    exp() %>%                                                            # exponentiate, since that is what rsd() is expecting
    rsd(xbar = 1, log_scale = TRUE)                                      # pass to rsd() (xbar is ignored when log_scale = TRUE)

  # Inter-Day CV
    # see comments above for additional details
  tables$prec[['Inter-Day CV']][i] <- tmp %>%
    filter(Assay == tables$prec$Assay[i]) %>%
    select(Sample_ID, inter_day_xbar, inter_day_std, inter_day_n) %>%
    unique() %>%
    with(sum(log(inter_day_std)^2 * inter_day_n) / sum(inter_day_n)) %>%
    sqrt() %>%
    exp() %>%
    rsd(xbar = 1, log_scale = TRUE)
  
  # Inter-Analyst CV
    # see comments above for additional details
  tables$prec[['Inter-Analyst CV']][i] <- tmp %>%
    filter(Assay == tables$prec$Assay[i]) %>%
    select(Sample_ID, inter_analyst_xbar, inter_analyst_std, inter_analyst_n) %>%
    unique() %>%
    with(sum(log(inter_analyst_std)^2 * inter_analyst_n) / sum(inter_analyst_n)) %>%
    sqrt() %>%
    exp() %>%
    rsd(xbar = 1, log_scale = TRUE)
}

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

tables$covr <- covr_sum %>%
  group_by(Assay, Sample_ID) %>%
  # only need weighted average of standard deviations for this calculation
  # see comments on Intra-Day CV for details on calculation (this line is dense, sorry)
  summarize(`Pct Error` = rsd(xbar = 1, exp(sqrt(sum(log(sdg)^2 * n) / sum(n))), 
                              log_scale = TRUE)) %>%
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
conj <- filter(conj, !is.na(acon) & !is.na(Assay) &
                 !(Assay == 'CoV2 N' & Interpretation == 'Negative')) %>%
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
antigen_sum <- filter(antigen_sum, !is.na(acon) & !is.na(Assay) &
                        !(Assay == 'CoV2 N' & Interpretation == 'Negative')) %>%
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
