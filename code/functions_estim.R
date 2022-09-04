#' @param dat Study data
#' @param numer_mod_forms List of model formulae needed to construct models for 
#'                        numerator of (stabilized) weights.
#' @param denom_mod_forms Same as `numer_mod_forms` but for denominator of
#'                        weights.
#' @param type If "first", get weights to estimate effect of INITIATING 
#'             treatment. If "all", get weights to estimate effect of 
#'             treatment strategy.
get_ip_weights_tvt <- function(dat, numer_mod_forms, denom_mod_forms, type = c("first", "all")) {
    type <- match.arg(type)
    stopifnot(length(numer_mod_forms)==length(denom_mod_forms))

    num_tp <- length(numer_mod_forms)

    ## Get propensity score if using type = "first"
    get_ps_for_initiation_wts <- function(dat, tp, ps) {
        treatment_this_time <- paste0("D", tp)
        if (tp >= 2) {
            treatment_prev_time <- paste0("D", tp-1)
            ifelse(
                dat[[treatment_this_time]]==1 & dat[[treatment_prev_time]]==1,
                1,
                ifelse(dat[[treatment_this_time]]==1, ps, 1-ps)
            )
        } else {
            ifelse(dat[[treatment_this_time]]==1, ps, 1-ps)
        }
    }
    
    ## Stabilized weights
    wts_numer_over_time_stab <- lapply(seq_along(numer_mod_forms), function(i) {
        mod_form <- numer_mod_forms[[i]]
        mod <- glm(mod_form, data = dat, family = "binomial")
        ps <- predict(mod, type = "response")
        if (type=="first") {
            get_ps_for_initiation_wts(dat = dat, tp = i, ps = ps)
        } else {
            treatment_this_time <- paste0("D", i)
            ifelse(dat[[treatment_this_time]]==1, ps, 1-ps)
        }
    })

    ## Unstabilized weights
    wts_numer_over_time_unstab <- lapply(denom_mod_forms, function(mod_form) {
        rep(1, nrow(dat))
    })

    # Estimate denominator of weights
    ps_over_time <- lapply(seq_along(denom_mod_forms), function(i) {
        mod_form <- denom_mod_forms[[i]]
        mod <- glm(mod_form, data = dat, family = "binomial")
        predict(mod, type = "response")
    })
    wts_denom_over_time <- lapply(seq_along(ps_over_time), function(i) {
        ps <- ps_over_time[[i]]
        if (type=="first") {
            get_ps_for_initiation_wts(dat = dat, tp = i, ps = ps)
        } else {
            treatment_this_time <- paste0("D", i)
            ifelse(dat[[treatment_this_time]]==1, ps, 1-ps)
        }
    })

    compute_ate_weights <- function(wts_numer_over_time, wts_denom_over_time) {
        log_wts_over_time <- lapply(seq_along(wts_numer_over_time), function(i) {
            numer_this_time <- wts_numer_over_time[[i]]
            denom_this_time <- wts_denom_over_time[[i]]
            log(numer_this_time) - log(denom_this_time)
        })
        log_wts <- Reduce("+", log_wts_over_time)
        exp(log_wts)
    }

    compute_gtatt_weights <- function(dat, wts_ate_unstab, ps_over_time, g, t, control_group = c("never", "notyet"), type = c("first", "all")) {
        ## Mark which units can serve as valid controls
        if (control_group=="never") {
            dat <- dat %>%
                mutate(is_valid_control = group==0)
        } else if (control_group=="notyet") {
            dat <- dat %>%
                mutate(is_valid_control = group==0 | group > max(g,t))
        }
        ## What treatment strategy corresponds to being in group g?
        treat_strat <- c(rep(0, g-1), rep(1, num_tp-g+1))
        dat_all_this_treat_strat <- matrix(rep(treat_strat, times = nrow(dat)), nrow = nrow(dat), ncol = length(treat_strat), byrow = TRUE)
        colnames(dat_all_this_treat_strat) <- paste0("D", seq_along(treat_strat))
        dat_all_this_treat_strat <- as_tibble(dat_all_this_treat_strat)

        ## Get propensity to follow this treatment strategy
        log_ps_prod_terms <- lapply(seq_along(ps_over_time), function(i) {
            ps <- ps_over_time[[i]]
            if (type=="first") {
                get_ps_for_initiation_wts(dat = dat_all_this_treat_strat, tp = i, ps = ps) %>% log()
            } else {
                treatment_this_time <- treat_strat[i]
                ifelse(treatment_this_time==1, ps, 1-ps) %>% log()
            }
        })
        ps_this_strat <- exp(Reduce("+", log_ps_prod_terms))

        ## Create weights
        gtatt_wts <- wts_ate_unstab*ps_this_strat
        ## Normalize weights for controls to sum to sample size for this group
        gtatt_wts[!dat$is_valid_control] <- 0
        n_this_group <- sum(dat$group==g)
        gtatt_wts <- gtatt_wts*n_this_group/sum(gtatt_wts)
        ## Treated units (under this treatment strategy) get weight = 1
        gtatt_wts[dat$group==g] <- 1
        
        gtatt_wts
    }

    ## Compute ATE weights
    wts_ate_stab <- compute_ate_weights(wts_numer_over_time_stab, wts_denom_over_time)
    wts_ate_unstab <- compute_ate_weights(wts_numer_over_time_unstab, wts_denom_over_time)

    ## Compute group-time ATT weights for different g, t, and control group
    df_wts_gtatt <- tidyr::crossing(
        group = seq_len(num_tp),
        time = seq_len(num_tp),
        ctrl_group = c("never", "notyet")
    ) %>% filter(time >= group)
    wts_gtatt <- lapply(seq_len(nrow(df_wts_gtatt)), function(r) {
        compute_gtatt_weights(dat, wts_ate_unstab, ps_over_time, g = df_wts_gtatt$group[r], t = df_wts_gtatt$time[r], control_group = df_wts_gtatt$ctrl_group[r], type = type)
    })
    df_wts_gtatt$wts <- wts_gtatt

    list(ate_stab = wts_ate_stab, ate_unstab = wts_ate_unstab, gtatt = df_wts_gtatt)
}

ipw_msm_fit <- function(msm_mod_form, dat, wts) {
    design <- svydesign(ids = ~0, weights = wts, data = dat)
    msm_fit <- svyglm(msm_mod_form, data = dat, design = design)
    msm_fit
}
