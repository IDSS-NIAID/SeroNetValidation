# helper functions

require(magrittr)


#' geo_mean
#' geo_sd
#' Geometric mean and standard deviation of a vector of data
#' 
#' @param x a vector of values to calculate the geometric mean
#' @param ... other parameters passed to mean and sd
#' 
#' @return The geometric mean
geo_mean <- function(x, ...)
{
    lx <- log(x)
    
    # calculate mean of logs, removing infinite values
    mean(lx[is.finite(lx)], ...) %>%
        exp()
}

geo_sd <- function(x, ...)
{
    lx <- log(x)
    
    # calculate sd of logs, removing infinite values
    sd(lx[is.finite(lx)], ...) %>%
        exp()
}


#' pct_err
#' Percent error of a an estimate
#' 
#' @param xbar A numeric value specifying the mean.
#' @param e A numeric value specifying the expected mean.
#' @param log_scale Logical. When TRUE, convert xbar and e to log scale prior to calculation
#' 
#' @return The percent error, delta.
pct_err <- function(xbar, e, log_scale = FALSE)
{
    if(log_scale)
    {
        retval <- abs((log(xbar) - log(e)) / log(e)) * 100
    }else{
        retval <- abs((xbar - e) / e) * 100
    }
    
    retval <- ifelse(is.na(retval), Inf, retval)
    
    return(retval)
}


#' rsd
#' Relative standard deviation
#' 
#' @param xbar Numeric. Mean
#' @param sd Numeric. Standard deviation
#' @param log_scale Logical. When TRUE, convert xbar and sd to log scale prior to calculation
#' 
#' @return The relative standard deviation
rsd <- function(xbar, sd, log_scale = TRUE)
{
    if(log_scale)
    {
        retval <- 100 * log(sd) / log(xbar)
    }else{
        retval <- 100 * sd / xbar
    }
    
    retval <- ifelse(is.na(retval), Inf, retval)
    
    return(retval)
}

#' summary_table_info
#' Get information for summary table 1
#' 
#' @param tables current list of tables
#' @param dat Data from input sheet
#' @param row Experiment name to identify correct row
summary_table_update <- function(tables, dat, row, methodology, acceptance)
{
    # summary table information
    i <- which(tables$table1$Experiment == row)
    tables$table1$`Samples (n)`[i] <- length(unique(dat$Sample_ID))
    tables$table1$Replicates[i] <- length(unique(dat$Replicate_per_sample))
    tables$table1$`Analyst(s)`[i] <- length(unique(dat$Analyst))
    tables$table1$`Days (Runs)`[i] <- length(unique(dat$Day))
    tables$table1$Methodology[i] <- methodology
    tables$table1$`Acceptance Criteria`[i] <- acceptance
    
    return(tables)
}


#' get_summary_table1
#' Summary table for LLOQ, ULOQ
#' 
#' @param tables current list of tables
#' @param test_name Name of current test (e.g. LLOQ, ULOQ)
#' @param stats Table of summary statistics
#' @param minmax Function to calculate min/max of statistics in stats_sum
#' @param which.minmax Corresponding which.min/max function to go along with minmax
#' @param pctl Percentile to calculate (0.025 for lower, 0.975 for upper)
#' @param dilution_factor Base dilution factor
get_summary_table1 <- function(tables, test_name, stats, minmax, which.minmax, pctl, dilution_factor)
{
    # summary table by sample ID and assay
    # filter by acceptance criteria
    igg <- filter(stats, delta <= 50 & rsd <= 30 & grepl('IgG', Assay)) %>%
        
        # calculate min/max concentration that passed - grouped by sample ID and analyst
        group_by(Sample_ID, Assay, Analyst) %>%
        summarize(statxbar = minmax(xbar, na.rm = TRUE),
                  dilFact = Dil_Factor[which.minmax(xbar)]) %>%
        ungroup() %>%
        
        # drop Analsyt initials and convert to 1/2
        mutate(tmp = factor(Analyst),
               Analyst = as.numeric(tmp)) %>%
        select(-tmp)
    
    igm <- filter(stats, delta <= 50 & rsd <= 30 & grepl('IgM', Assay)) %>%
        
        # calculate min/max concentration that passed - grouped by sample ID and analyst
        group_by(Sample_ID, Assay, Analyst) %>%
        summarize(statxbar = minmax(xbar, na.rm = TRUE),
                  dilFact = Dil_Factor[which.minmax(xbar)]) %>%
        ungroup() %>%
        
        # drop Analsyt initials and convert to 1/2
        mutate(tmp = factor(Analyst),
               Analyst = as.numeric(tmp)) %>%
        select(-tmp)
    
    
    tables[[test_name]] <- 
        bind_rows(pivot_wider(igg, id_cols = c(Sample_ID, Assay),
                              names_from = Analyst, values_from = statxbar,
                              names_prefix = 'Analyst '), 
                  pivot_wider(igm, id_cols = c(Sample_ID, Assay), 
                              names_from = Analyst, values_from = statxbar,
                              names_prefix = 'Analyst '))
    
    tables[[paste0(test_name, '_sum')]] <- 
        tibble(Measure = c('Mean_g', 'Median', 'SD_g', 'Min', 'Max', 
                           paste(ifelse(pctl < .5, 'Lower', 'Upper'), '95% CI'),
                           'Dilution Factor', paste(toupper(test_name), '(AU/mL)')),
               `IgG Statistic` = with(igg, c(geo_mean(statxbar, na.rm = TRUE),
                                             median(statxbar, na.rm = TRUE),
                                             geo_sd(statxbar, na.rm = TRUE),
                                             min(statxbar, na.rm = TRUE),
                                             max(statxbar, na.rm = TRUE),
                                             geo_mean(statxbar, na.rm = TRUE) + qnorm(pctl)*geo_sd(statxbar),
                                             dilution_factor,
                                             dilution_factor * (geo_mean(statxbar) + qnorm(pctl)*geo_sd(statxbar)))),
               `IgM Statistic` = with(igm, c(geo_mean(statxbar, na.rm = TRUE),
                                             median(statxbar, na.rm = TRUE),
                                             geo_sd(statxbar, na.rm = TRUE),
                                             min(statxbar, na.rm = TRUE),
                                             max(statxbar, na.rm = TRUE),
                                             geo_mean(statxbar, na.rm = TRUE) + qnorm(pctl)*geo_sd(statxbar),
                                             dilution_factor,
                                             dilution_factor * (geo_mean(statxbar) + qnorm(pctl)*geo_sd(statxbar)))))

    return(tables)
}

#' get.stability.max
#' 
#' @param dat Tibble with stability data. Must contain Assay and delta variables.
get.stability.max <- function(dat)
{
    paste0('≤',
           with(dat, round(c(max(delta[grepl('IgG', Assay)]),
                             max(delta[grepl('IgM', Assay)])),
                           digits = 1)))
}


#' get_rsd_by_analyst
#' Get RSD measures by analyst and between analyst
#' 
#' @param dat Data frame with the following variables: Assay, Sample_ID, Analyst, theo_con, acon
#' 
#' @return A summary data frame with the following variables: Assay, Sample_ID, Analyst, theo_con, meang, delta
get_rsd_by_analyst <- function(dat)
{
    dat_by_analyst <- group_by(dat, Assay, Sample_ID, Analyst) %>%
        summarize(meang = geo_mean(acon, na.rm = TRUE),
                  sdg = geo_sd(acon, na.rm = TRUE),
                  rsd = rsd(meang, sdg)) %>%
        ungroup()
    
    dat_by_sample <- group_by(dat, Assay, Sample_ID) %>%
        summarize(meang = geo_mean(acon, na.rm = TRUE),
                  sdg = geo_sd(acon, na.rm = TRUE),
                  rsd = rsd(meang, sdg)) %>%
        ungroup()
    
    full_join(dat_by_analyst,
              dat_by_sample,
              c("Assay", "Sample_ID", "meang", "sdg", "rsd")) %>%
        arrange(Assay, Sample_ID, Analyst) %>%
        mutate(Analyst = ifelse(is.na(Analyst), 'Between Analysts', Analyst)) %>%
        
        return()
}


#' get_pct_err_by_analyst
#' Get % error measures by analyst and between analyst
#' 
#' @param dat Data frame with the following variables: Assay, Sample_ID, Analyst, theo_con, acon
#' 
#' @return A summary data frame with the following variables: Assay, Sample_ID, Analyst, theo_con, meang, delta
get_pct_err_by_analyst <- function(dat)
{
    dat_by_analyst <- group_by(dat, Assay, Sample_ID, Analyst, theo_con) %>%
        summarize(meang = mean(acon, na.rm = TRUE),
                  delta = pct_err(meang, unique(theo_con))) %>%
        ungroup()
    
    dat_by_sample <- group_by(dat, Assay, Sample_ID, theo_con) %>%
        summarize(meang = mean(acon, na.rm = TRUE),
                  delta = pct_err(meang, unique(theo_con))) %>%
        ungroup()
    
    full_join(dat_by_analyst,
              dat_by_sample,
              c("Assay", "Sample_ID", "theo_con", "meang", "delta")) %>%
        arrange(Assay, Sample_ID, Analyst) %>%
        mutate(Analyst = ifelse(is.na(Analyst), 'Between Analysts', Analyst)) %>%
        
        return()
}
