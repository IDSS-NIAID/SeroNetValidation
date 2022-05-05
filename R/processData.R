# processData.R
# should be called from report.Rmd

# pct error
base_pct_error <- 30
lloq_pct_error <- 50
uloq_pct_error <- 50
lin_pct_error <- 50
sens_pct_error <- 50
acc_pct_error <- 20
stab_pct_error <- 20

# rsd pct
base_rsd_pct <- 25
lloq_rsd_pct <- base_rsd_pct
uloq_rsd_pct <- base_rsd_pct
lin_rsd_pct <- 30
sens_rsd_pct <- base_rsd_pct
prec_rsd_pct <- 25
covr_rsd_pct <- 20

# other constants
cutpt_bound <- 0.975
spec_inhibit <- 90
spec_non_inhibit <- 25


########
# LLOQ #
########

# read in and process LLOQ data
lloq <- read_excel(f, sheet = 'LLOQ', na = c('', 'Sample did not dilute down properly-can not use data'))
names(lloq)[names(lloq) == diff$lloq_acon] <- 'acon'
names(lloq)[names(lloq) == diff$lloq_assay] <- 'Assay'

lloq <- lloq %>%
    # drop out of range samples and convert to numeric
    filter(acon != 'Range?' & !is.na(acon)) %>%
    mutate(acon = as.numeric(acon)) %>%
    
    # base concentration
    group_by(Assay, Sample_ID, Analyst) %>%
    mutate(base_dil = min(Dil_Factor),
           base_con = geo_mean(acon[Dil_Factor == min(Dil_Factor)])) %>%
    ungroup()

# update summary table 1
tables <- summary_table_update(tables, lloq, 'LLOQ',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution', 
                               paste0('Percent Error ≤ ', lloq_pct_error, '%; ',
                                      'RSD ≤ ', lloq_rsd_pct, '%'))

# calculate statistics for each group
lloq_sum <- group_by(lloq, Assay, Sample_ID, Analyst, Dil_Factor) %>%
    
    # theoretical concentration
      # base_con * base_dil /  # base concentration -> normalize to dilution factor of 1
      # Dil_Factor             # divide by dilution factor to calculate expected concentration for each dilution
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
lloq_thresh <- lloq_sum %>%
  filter(delta <= lloq_pct_error & rsd <= lloq_rsd_pct) %>%
  group_by(Sample_ID, Assay, Analyst) %>%
  summarize(statxbar = min(xbar,na.rm = TRUE),
            statsd = std[which.min(xbar)],
            dilFact = Dil_Factor[which.min(xbar)],
            n = n[which.min(xbar)]) %>%
  ungroup() %>%
  group_by(Assay) %>%
  summarize(lloq = map_dbl(1, ~ exp(metamean(n, log(statxbar), log(statsd))$TE.fixed))) %>%
  ungroup()

# if you see a warning message for this or the update of `tables` on the next line,
# (warning in metamean about non-positive standard deviations) then there is probably 
# one or more samples with only a single sample (n == 1). This is OK - it will get
# skipped.

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
uloq <- read_excel(f, sheet = 'ULOQ')
names(uloq)[names(uloq) == diff$uloq_acon] <- 'acon'
names(uloq)[names(uloq) == diff$uloq_assay] <- 'Assay'

uloq <- uloq %>%
    # drop out of range samples and convert to numeric
    filter(acon != 'Range?' & !is.na(acon)) %>%
    mutate(acon = as.numeric(acon)) %>%
  
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
                               paste0('Percent Error ≤ ', uloq_pct_error, '%; ',
                                      'RSD ≤ ', uloq_rsd_pct, '%'))

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
  filter(delta <= uloq_pct_error & rsd <= uloq_rsd_pct) %>%
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

lloq$lab <- 'lloq'
uloq$lab <- 'uloq'

lin <- full_join(lloq, uloq, c("Assay", "Sample_ID", "Analyst", "Dil_Factor", "Day",
                               "Plate_per_Day", "Sample_per_plate", "Replicate_per_sample",
                               "acon", "% of first dilution", "base_con", "lab", 
                               "base_dil")) %>%
  
  mutate(theo_con = base_con * base_dil / Dil_Factor)
  
# update summary table 1
tables <- summary_table_update(tables, lin, 'Linearity',
                               'Calculate Geometric Mean, RSD, and Percent Error for each dilution',
                               paste0('Percent Error ≤ ', lin_pct_error, '%; RSD ≤ ', lin_rsd_pct, '%'))

# calculate statistics for each group
lin_sum <- group_by(lin, Assay, Sample_ID, Analyst, Dil_Factor, lab) %>%

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
              id = unique(paste(Assay, Sample_ID, Analyst, Dil_Factor, lab))) %>%

    ungroup()

##### plots #####

# start by filtering out parts of lin that we aren't using
lin <- lin %>%
  mutate(id = paste(Assay, Sample_ID, Analyst, Dil_Factor, lab),
         keep = id %in% filter(lin_sum, delta <= lin_pct_error & rsd <= lin_rsd_pct)$id,
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
                                 tmp <- filter(lin, Assay == .x & keep & is.finite(lacon))
                                 if(length(tmp$lacon) > 3)
                                 {
                                   return(lme(lacon ~ ltheo_con,
                                              random = ~ 1 | Sample_ID,
                                              data = tmp,
                                              na.action = na.omit))
                                  }else{
                                    return(NULL)
                                  }
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


figures$linearity_OvE_concentration <-
  map(unique(lin_sum$Assay), ~ filter(lin_sum, delta <= lin_pct_error & rsd <= lin_rsd_pct & Assay == .x) %>%
        ggplot(aes(theo_con, xbar)) +

        # plot points on log10 scale
        geom_point() +
        scale_x_log10() +
        scale_y_log10() +

        # add spread for each mean
        geom_errorbar(aes(ymin = min_acon, ymax = max_acon), width = 0) +

        # add trend line
        geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
        geom_abline(slope = 1, intercept = 0) +

        # labels
        annotate('text', x = 0, y = Inf,
                 label = paste(' slope =', round(filter(lin_assay_sum, Assay == .x)$slope, 2)),
                 hjust = 0, vjust = 1) +
        ylab('Titer (AU/mL)') +
        xlab('Theoretical Concentration') +
        ggtitle(.x))

names(figures$linearity_OvE_concentration) <- unique(lin_sum$Assay)

figures$linearity_OvE_concentration_bad <-
  map(unique(lin_sum$Assay), ~ filter(lin_sum,
                                      Assay == .x &
                                      ((delta < lin_pct_error & rsd < lin_rsd_pct) |
                                       (theo_con > filter(lloq_thresh, Assay == .x)$lloq &
                                        theo_con < filter(uloq_thresh, Assay == .x)$uloq))) %>%
        ggplot(aes(theo_con, xbar, color = delta < lin_pct_error & rsd < lin_rsd_pct)) +

        # plot points on log10 scale
        geom_point() +
        scale_x_log10() +
        scale_y_log10(labels = scales::comma) +

        # add trend line
        geom_smooth(method = 'lm', se = TRUE) +
        geom_abline(slope = 1, intercept = 0) +

        # add cutoff lines
        # geom_vline(xintercept = filter(lloq_thresh, Assay == .x)$lloq, linetype = 2) +
        # geom_vline(xintercept = filter(uloq_thresh, Assay == .x)$uloq, linetype = 2) +

        # labels
        ylab('Titer (AU/mL)') +
        xlab('Theoretical Concentration') +
        ggtitle(.x))

names(figures$linearity_OvE_concentration_bad) <- unique(lin_sum$Assay)

##### final table #####
# calculate linearity summaries by Assay, Sample_ID
tables$lin <- lin_sum %>%
  mutate(pass = delta < lin_pct_error & rsd < lin_rsd_pct) %>%
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

cutpt <- read_excel(f, sheet = diff$cutpoint_sheet, na = c('', 'No value'))
names(cutpt)[names(cutpt) == diff$cutpoint_acon] <- 'acon'
names(cutpt)[names(cutpt) == diff$cutpoint_assay] <- 'Assay'

# update summary table 1
tables <- summary_table_update(tables, cutpt, 'Cutpoint', 'Calculate Geometric Mean and 95% CI', '')
    
cutpt_sum <- group_by(cutpt, Assay, Sample_ID) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              upper_95_CI = qtl_limit(xbar, std, n = n, qtl = cutpt_bound, log_scale = TRUE)) %>%
    ungroup()

tables$cutpt <- group_by(cutpt, Assay) %>%
    summarize(n = sum(!is.na(acon)),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              `Upper 95% Confidence Bound` = qtl_limit(xbar, std, n = n, qtl = cutpt_bound, log_scale = TRUE)) %>%
    ungroup() %>%
    select(-xbar, -std)


###############
# Sensitivity #
###############

#### LLOQ ####
sens <- read_excel(f, sheet = diff$sens_sheet1, na = c('', 'Plate failed', 'Range?'))
names(sens)[names(sens) == diff$sens_acon] <- 'acon'
names(sens)[names(sens) == diff$sens_assay] <- 'Assay'

sens <- filter(sens, !is.na(acon)) %>%
  mutate(lacon = log(acon),
         ldil_factor = log(Dil_Factor))

# update summary table 1
tables <- summary_table_update(tables, sens, 'Sensitivity (LLOQ Challenge)',
                               'Calculate Geometric Mean, RSD, and Percent Error for each sample concentration level',
                               paste0('Samples at LLOQ and higher must pass these criteria to accept LLOQ: Percent Error ≤ ', sens_pct_error, '%; RSD ≤ ', sens_rsd_pct, '%'))

# run linear model for lloq sensitivity to estimate expected acon
sens_lloq_models <- sens %>%
  filter(acon >= unlist(tables$lloq_sum[1,'Geometric Mean']),
         acon <= unlist(tables$uloq_sum[1,'Geometric Mean'])) %>%
  
  group_by(Assay, Sample_ID) %>%
  summarize(model = map(1, ~ lm(lacon ~ ldil_factor))) %>%
  ungroup()

# get predictions and summary statistics
sens_sum <- sens %>%
    
    group_by(Assay, Sample_ID, Dil_Factor) %>%
    summarize(theo_con = pmap_dbl(list(assay = Assay, id = Sample_ID, dil = log(Dil_Factor)),
                                  function(assay, id, dil)
                                  {
                                    with(filter(sens_lloq_models, Assay == assay & Sample_ID == id),
                                         predict(model[[1]], newdata = data.frame(ldil_factor = dil))) %>%
                                      exp()
                                  }),
              xbar = geo_mean(acon, na.rm = TRUE),
              std = geo_sd(acon, na.rm = TRUE),
              rsd = rsd(xbar, std, log_scale = TRUE),
              delta = pct_err(xbar, theo_con)) %>%
    ungroup()

tables$sens <- filter(sens_sum, delta <= sens_pct_error & rsd <= sens_rsd_pct) %>%
    group_by(Assay, Sample_ID) %>%
    summarize(`Pct Error` = delta[which.min(theo_con)],
              `CV%` = rsd[which.min(theo_con)],
              `Concentration (AU/mL)` = xbar[which.min(theo_con)]) %>%
    ungroup()

#### ULOQ ####
sens2 <- read_excel(f, sheet = diff$sens_sheet2, na = c('', 'Plate failed', 'Range?'))
names(sens2)[names(sens2) == diff$sens_acon] <- 'acon'
names(sens2)[names(sens2) == diff$sens_assay] <- 'Assay'

sens2 <- filter(sens2, !is.na(acon)) %>%
  mutate(lacon = log(acon),
         ldil_factor = log(Dil_Factor))

# update summary table 1
tables <- summary_table_update(tables, sens2, 'Sensitivity (ULOQ Challenge)',
                               'Calculate Geometric Mean, RSD, and Percent Error for each sample concentration level',
                               paste0('Samples at LLOQ and higher must pass these criteria to accept LLOQ: Percent Error ≤ ', sens_pct_error, '%; RSD ≤ ', sens_rsd_pct, '%'))

# run linear model for uloq sensitivity to estimate expected acon
sens_uloq_models <- sens2 %>%
  filter(acon >= unlist(tables$lloq_sum[1,'Geometric Mean']),
         acon <= unlist(tables$uloq_sum[1,'Geometric Mean'])) %>%
  
  group_by(Assay, Sample_ID) %>%
  summarize(model = map(1, ~ lm(lacon ~ ldil_factor))) %>%
  ungroup()

# get predictions and summary statistics
sens2_sum <- sens2 %>%
  
  group_by(Assay, Sample_ID, Dil_Factor) %>%
  summarize(theo_con = pmap_dbl(list(assay = Assay, id = Sample_ID, dil = log(Dil_Factor)),
                                function(assay, id, dil)
                                {
                                  with(filter(sens_uloq_models, Assay == assay & Sample_ID == id),
                                       predict(model[[1]], newdata = data.frame(ldil_factor = dil))) %>%
                                    exp()
                                }),
            xbar = geo_mean(acon, na.rm = TRUE),
            std = geo_sd(acon, na.rm = TRUE),
            rsd = rsd(xbar, std, log_scale = TRUE),
            delta = pct_err(xbar, theo_con)) %>%
  ungroup()

tables$sens2 <- filter(sens2_sum, delta <= sens_pct_error & rsd <= sens_rsd_pct) %>%
  group_by(Assay, Sample_ID) %>%
  summarize(`Pct Error` = delta[which.max(theo_con)],
            `CV%` = rsd[which.max(theo_con)],
            `Concentration (AU/mL)` = xbar[which.max(theo_con)]) %>%
  ungroup()

# debugging
if(FALSE)
{
  sens2_sum %>%
    filter(delta <= sens_pct_error & rsd <= sens_rsd_pct) %>% # this is the step where we loose everything
    group_by(Sample_ID) %>%
    summarize(max(xbar))
  
  sens2 %>%
    #filter(Sample_ID == 'RDP0977_C1') %>%
    #mutate(qc = ifelse(delta <= sens_pct_error & rsd <= sens_rsd_pct, 'keep', 'filter')) %>%
    mutate(Day = as.factor(Day)) %>%
    ggplot(aes(Dil_Factor, acon, color = Day)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_log10() +
    scale_y_log10() +
    #geom_smooth(method = 'lm') +
    scale_color_manual(values = 1:6) +
    facet_wrap(~ Sample_ID)
}


##########################
# Accuracy and Precision #
##########################

#### Accuracy ####
acc <- read_excel(f, sheet = diff$acc_sheet)
names(acc)[names(acc) == diff$acc_assay] <- 'Assay'
names(acc)[names(acc) == diff$acc_acon] <- 'acon'
names(acc)[names(acc) == diff$acc_theo] <- 'theo_con'

acc <- filter(acc, !is.na(Assay))

# update summary table 1
tables <- summary_table_update(tables, acc, 'Accuracy',
                               'Calculate Geometric Mean and Percent Error', 
                               paste0('Percent Error ≤ ', acc_pct_error, '%'))


tables$acc <- acc %>%
    # convert analyst initials to integers in {1,2,...}
    group_by(Assay, Sample_ID) %>%
    mutate(Analyst = as.character(as.numeric(as.factor(Analyst)))) %>%
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

tables$acc <- acc %>%
  group_by(Assay, Sample_ID) %>%
  summarize(Analyst = 'All',
            xbar = geo_mean(acon, na.rm = TRUE),
            sd = geo_sd(acon, na.rm = TRUE),
            `CV%` = rsd(xbar, sd, log_scale = TRUE),
            `Pct Error` = pct_err(xbar, unique(theo_con))) %>%
  ungroup() %>%
  rename(`Geometric Mean (AU/mL)` = xbar) %>%
  select(-sd) %>%
  full_join(tables$acc, by = c("Assay", "Sample_ID", "Analyst", 
                               "Geometric Mean (AU/mL)", 'CV%',
                               'Pct Error'))
  

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
names(prec)[names(prec) == diff$prec_assay] <- 'Assay'
prec <- mutate(prec, # I think this is a typo - should be '<8'
               acon = ifelse(acon %in% c('<8', '<7'), 4, as.numeric(acon))) %>%
  filter(!is.na(acon))

# update summary table 1
tables <- summary_table_update(tables, prec, 'Precision',
                               'Calculate Geometric Mean and RSD for Intra-plate, Inter-plate, and Inter-Analyst', 
                               paste0('RSD ≤ ', prec_rsd_pct, '% for Intra-plate, Inter-plate, and Inter- Analyst'))

tmp <- prec %>%
    # filter such that delta ≤ 25%
    group_by(Assay, Sample_ID, Analyst, Day) %>%
    mutate(xbar = geo_mean(acon, na.rm = TRUE),
           std = geo_sd(acon, na.rm = TRUE),
           rsd = rsd(xbar, std, log_scale = TRUE)) %>%
    ungroup() %>%
    filter(rsd <= prec_rsd_pct) %>%
    
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


###############
# Specificity #
###############

spec <- read_excel(f, sheet = diff$spec_sheet)
names(spec)[names(spec) == diff$spec_acon] <- 'acon'
names(spec)[names(spec) == diff$spec_assay] <- 'Assay'

spec <- spec %>%
  mutate(acon = ifelse(acon == '<8', 4, as.numeric(acon)),
         Treatment = strsplit(Sample_ID, '_', fixed = TRUE) %>%
           map_chr(~ .x[2]))

for(i in 1:nrow(spec))
{
  if(grepl('_control', spec$Sample_ID[i], ignore.case = TRUE))
  {
    spec$Sample_ID[i] <- gsub('_control', '', spec$Sample_ID[i], ignore.case = TRUE)
  }else{
    spec$Sample_ID[i] <- spec$Sample_ID[i - 1]
  }
}
           
# update summary table 1
tables <- summary_table_update(tables, spec, 'Specificity',
                               'Calculate Geometric Mean and RSD',
                               paste0('1) ≥', spec_inhibit, '% of sample concentration must be inhibited with type-specific VLPs. 2) Non-type-specific inhibition must be ≤ ', spec_non_inhibit, '%.'))

# expected concentration
spec_sum <- group_by(spec, Assay, Sample_ID) %>%
  mutate(theo_con = geo_mean(acon[Treatment == 'Control'])) %>%
  ungroup() %>%

  # don't need these anymore
  filter(Treatment != 'Control') %>%

  # percent error
  group_by(Assay, Sample_ID, Treatment, theo_con) %>%
  summarize(`geometric mean` = geo_mean(acon),
            `Pct Error` = pct_err(`geometric mean`, unique(theo_con))) %>%
  ungroup()


tables$spec <- spec_sum


##############
# Carry-over #
##############

covr <- read_excel(f, sheet = diff$covr_sheet)
names(covr)[names(covr) == diff$covr_acon] <- 'acon'
names(covr)[names(covr) == diff$covr_sample_id] <- 'Sample_ID'
names(covr)[names(covr) == diff$covr_assay] <- 'Assay'

covr <- covr %>%
  mutate(acon = ifelse(acon == '<8', 4, as.numeric(acon))) %>%
  filter(!is.na(acon)) %>%
  
  # some of these ended up with replicate number added on the end of the sample ID
  mutate(Sample_ID = strsplit(Sample_ID, '_', fixed = TRUE) %>%
           map_chr(~ .x[1]))

# update summary table 1
tables <- summary_table_update(tables, covr, 'Carry-over',
                               'Calculate Geometric Mean and RSD', 
                               paste0('RSD ≤ ', covr_rsd_pct, '% for Low and High sample. The Neg_Serum sample must be ≤ LLOQ.'))

covr_sum <- get_rsd_by_analyst(covr)

tables$covr <- covr_sum %>%
  group_by(Assay, Sample_ID) %>%
  # only need weighted average of standard deviations for this calculation
  # see comments on Intra-Day CV for details on calculation (this line is dense, sorry)
  summarize(`geometric mean` = meang,
            `Pct Error` = rsd(xbar = 1, exp(sqrt(sum(log(sdg)^2 * n) / sum(n))), 
                              log_scale = TRUE)) %>%
  ungroup()


#############
# Stability #
#############

if(length(diff$stab_sheet) == 1)
{
  stability <- read_excel(f, sheet = diff$stab_sheet) %>%
  
    # drop 1X, 5X, and 10X from Treatment so we can identify samples
    mutate(test = gsub('\\d+X\\s', '', Treatment))

  names(stability)[names(stability) == diff$stab_acon] <- 'acon'
  names(stability)[names(stability) == diff$stab_assay] <- 'Assay'
  
  stability <- filter(stability, !is.na(acon)) %>%
    
    mutate(lvl = substr(Sample_ID, 1, 1), # figure out concentration level (low/high)
           trt = Sample_ID,               # the sample id is really the treatment here
           Sample_ID = Treatment) %>%     # when grouping, we really want to group by "Treatment" in the spreadsheet - swap labels for use in get_pct_err_by_analyst()
    
    # theoretical concentration
    group_by(Assay, lvl) %>%
    mutate(theo_con = median(acon[Sample_ID == ifelse(grepl('Freeze/Thaw', Sample_ID), '1X Freeze/Thaw', 'Control')])) %>%
    ungroup() %>%
    
    filter(! Sample_ID %in% c('Control', '1X Freeze/Thaw')) %>%
    
    # now that everything else is sorted out, this is really where this belongs
    mutate(Treatment = Sample_ID)
  
}else{
  # combine all sheets into one analysis
  stability <- map_df(diff$stab_sheet, ~ 
                      {
                        # read in each sheet
                        tmp <- read_excel(f, sheet = .x)
                        names(tmp)[names(tmp) == diff$stab_acon] <- 'acon'
                        names(tmp)[names(tmp) == diff$stab_assay] <- 'Assay'
                        
                        tmp %>%
                          
                          # remove extra crap from Sample_ID and properly label Treatment
                          mutate(lvl = strsplit(Sample_ID, ' ', fixed = TRUE) %>%
                                   map_chr(~ .x[length(.x)]),
                                 Sample_ID = strsplit(Sample_ID, ' ', fixed = TRUE) %>%
                                   map_chr(~ .x[1]) %>%
                                   gsub(pattern = 'Heat_', replacement = '', fixed = TRUE),
                                 Treatment = ifelse(is.na(Treatment) | lvl == '(1X)',
                                                    'control', Treatment),
                                 
                                 # flag different stability experiments
                                 test = strsplit(.x, '_', fixed = TRUE) %>% 
                                   map_chr(~ .x[length(.x)])) %>%
                          
                          return()
                      }) %>%
  
      # theoretical concentration
      group_by(Sample_ID, test) %>%
      mutate(theo_con = geo_mean(acon[Treatment == 'control'], na.rm = TRUE)) %>%
      ungroup() %>%
      filter(Treatment != 'control')
}

stability_sum <- get_pct_err_by_analyst(stability)


tables$stability <- stability_sum %>%
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
tables$table1$`Acceptance Criteria`[i] <- paste0('Percent Error ≤ ', stab_pct_error, '%')



##########################
# Lot-to-Lot Comparisons #
##########################

##### Lot to Lot Conjugate #####

conj <- read_excel(f, sheet = diff$conj_sheet)
names(conj)[names(conj) == diff$conj_lot] <- 'lot'
names(conj)[names(conj) == diff$conj_acon] <- 'acon'
names(conj)[names(conj) == diff$conj_assay] <- 'Assay'

conj <- conj %>%
  mutate(acon = ifelse(acon == '<8', 4, as.numeric(acon)),
         lot = tolower(lot)) %>%
  filter(!is.na(acon) & !is.na(Assay))

# update summary table 1
tables <- summary_table_update(tables, conj, 'Stability (Critical Reagent Lot Change)',
                               'Calculate Geometric Mean and Percent Error', 
                               paste0('Percent Error ≤ ', stab_pct_error, '%'))

# double check that there are two lots per sample
if(!(group_by(conj, Assay, Sample_ID) %>%
     summarize(n = length(acon),
               e = length(unique(lot)) == 2) %>%
     ungroup() %>%
     select(e) %>%
     unlist() %>%
     all()))
  stop('Not all Samples are tested with exactly two lots')

# calculate difference between two lots
conj_sum <- conj %>%
  group_by(Assay, Sample_ID, lot) %>%
  summarize(acon = geo_mean(acon, na.rm = TRUE)) %>%
  ungroup() %>%
  
  group_by(Assay, Sample_ID) %>%
  summarize(lotA = acon[1],
            lotB = acon[2],
            delta = pct_err(lotA, lotB)) %>%
  ungroup()


##### Lot to Lot Antigen #####

antigen <- read_excel(f, sheet = diff$antigen_sheet)
names(antigen)[names(antigen) == diff$antigen_lot] <- 'lot'
names(antigen)[names(antigen) == diff$antigen_acon] <- 'acon'
names(antigen)[names(antigen) == diff$antigen_assay] <- 'Assay'

antigen <- antigen %>%
  mutate(acon = ifelse(acon == '<8', 4, as.numeric(acon)),
         lot = tolower(lot)) %>%
  filter(!is.na(acon) & !is.na(Assay))

# double check that there are two lots per sample
if(!(group_by(antigen, Assay, Sample_ID) %>%
     summarize(n = length(acon),
               e = length(unique(lot)) == 2) %>%
     ungroup() %>%
     select(e) %>%
     unlist() %>%
     all()))
  stop('Not all Samples are tested with exactly two lots')

# calculate difference between two lots
antigen_sum <- antigen %>%
  group_by(Assay, Sample_ID, lot) %>%
  summarize(acon = geo_mean(acon, na.rm = TRUE)) %>%
  ungroup() %>%
  
  group_by(Assay, Sample_ID) %>%
  summarize(lotA = acon[1],
            lotB = acon[2],
            delta = pct_err(lotA, lotB)) %>%
  ungroup()

tables$l2l <- filter(conj_sum, delta > 0) %>%
    group_by(Assay) %>%
    summarize(`Conjugate Passing` = sum(delta <= stab_pct_error),
              `Conjugate Total` = length(delta),
              fixme = paste0(round(`Conjugate Passing` / `Conjugate Total`, 3) * 100, "%")) %>%
    ungroup()
names(tables$l2l)[names(tables$l2l) == 'fixme'] <- paste0('Conjugate w/ Error ≤ ', stab_pct_error, '%')

tables$l2l <- filter(antigen_sum, delta > 0) %>%
    group_by(Assay) %>%
    summarize(`Antigen Passing` = sum(delta <= stab_pct_error),
              `Antigen Total` = length(delta),
              fixme = paste0(round(`Antigen Passing` / `Antigen Total`, 3) * 100, "%")) %>%
    ungroup() %>%
    right_join(tables$l2l, by = 'Assay')
names(tables$l2l)[names(tables$l2l) == 'fixme'] <- paste0('Antigen w/ Error ≤ ', stab_pct_error, '%')


#######################
# Single vs Multiplex #
#######################

# svm <- read_excel(f, sheet = 'Single vs Multiplex') %>%
#   mutate(multiplex = grepl('multiplex', `Result File`))
# names(svm)[names(svm) == diff$svm_acon] <- 'acon'
# 
# # update summary table 1
# tables <- summary_table_update(tables, svm, 'Single vs Multiplex',
#                                'Calculate Geometric Mean and Percent Error', 
#                                'Percent Error ≤ 25%')
# 
# # comparison of the two assays
# svm_sum <- svm %>%
#   group_by(Assay, Sample_ID, Analyst) %>%
#   summarize(single    = geo_mean(acon[!multiplex]),
#             multiplex = geo_mean(acon[ multiplex]),
#             delta     = pct_err(multiplex, single)) %>%
#   ungroup()
# 
# tables$svm <- svm_sum %>%
#   group_by(Assay) %>%
#   summarize(`Total samples` = length(delta),
#             `Samples passing` = sum(delta <= 25),
#             `% passing` = paste0(round(`Samples passing` / `Total samples`, 3) * 100, '%')) %>%
#   ungroup()
