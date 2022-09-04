library(tidyverse)
library(survey)
library(did)

source("functions_sim_data.R")
source("functions_estim.R")

## Number of simulation repeats
num_reps <- 100

## Number of units
samp_sizes <- c(300, 3000)

## Number of time periods with available data
num_time_periods <- 4

## Strength of treatment-confounder feedback (coefficient on D_t-1 in the structural equation for X_t)
# treat_conf_strengths <- c(0, -0.05, -0.5, 0.05, 0.5)

## Multipliers of the X -> D and U/X -> Y coefficients
multipliers <- 1:3

## Code running times for 1 rep of 1 multiplier
## 45.785 for samp_sizes <- c(300, 3000) ------------> 3.8 hours
## 108.816 for samp_sizes <- c(300, 3000, 6000) -----> 9.1 hours
## 344.178 for samp_sizes <- c(300, 3000, 30000) ----> 28.7 hours

extract_info_att_gt <- function(obj) {
    ## Group-time atts
    res_gt <- tibble(group = obj$group, time = obj$t, estimate = obj$att, se = obj$se)

    ## Aggregating group-time ATTs by event time
    agg_dyna <- aggte(obj, type = "dynamic")
    res_agg_dyna <- tibble(group = "dynamic", time = agg_dyna$egt, estimate = agg_dyna$att.egt, se = agg_dyna$se.egt)

    bind_rows(
        res_gt %>% mutate(group = as.character(group)),
        res_agg_dyna
    )
}

extract_info_msm <- function(obj, gtatt = FALSE) {
    if (gtatt) {
        tibble(coeff = names(coefficients(obj)), estimate = coefficients(obj), se = SE(obj)) %>% filter(coeff != "(Intercept)")
    } else {
        tibble(coeff = names(coefficients(obj)), estimate = coefficients(obj), se = SE(obj))
    }
}

estimate_effects <- function(sim_data) {
    ## Insert ID and group variables (group as in group-time ATT for use in
    ## the did package)
    num_times_treated <- sim_data %>% select(starts_with("D")) %>% as.matrix() %>% rowSums()
    treat_group <- num_time_periods - num_times_treated + 1
    treat_group[treat_group==(num_time_periods+1)] <- 0
    sim_data <- sim_data %>%
        mutate(id = seq_len(nrow(sim_data)), group = treat_group)


    ## Effect estimation in the TVT framework
    denom_mod_forms <- list(
        D1 ~ X1,
        D2 ~ D1 + X1 + X2,
        D3 ~ D1 + D2 + X1 + X2 + X3,
        D4 ~ D1 + D2 + D3 + X1 + X2 + X3 + X4
    )
    numer_mod_forms <- list(
        D1 ~ 1,
        D2 ~ D1,
        D3 ~ D1 + D2,
        D4 ~ D1 + D2 + D3
    )
    
    ## Get weights up through time 1, ..., up through time 4
    ipw_tvt_wts_t1 <- get_ip_weights_tvt(dat = sim_data, numer_mod_forms = numer_mod_forms[1], denom_mod_forms = denom_mod_forms[1], type = "first")
    ipw_tvt_wts_t2 <- get_ip_weights_tvt(dat = sim_data, numer_mod_forms = numer_mod_forms[1:2], denom_mod_forms = denom_mod_forms[1:2], type = "first")
    ipw_tvt_wts_t3 <- get_ip_weights_tvt(dat = sim_data, numer_mod_forms = numer_mod_forms[1:3], denom_mod_forms = denom_mod_forms[1:3], type = "first")
    ipw_tvt_wts_t4 <- get_ip_weights_tvt(dat = sim_data, numer_mod_forms = numer_mod_forms, denom_mod_forms = denom_mod_forms, type = "first")

    ## Using one set of weights for all time point models
    ipw_tvt_fit1 <- ipw_msm_fit(msm_mod_form = Y1 ~ D1, dat = sim_data, wts = ipw_tvt_wts_t4$ate_stab) %>% extract_info_msm()
    ipw_tvt_fit2 <- ipw_msm_fit(msm_mod_form = Y2 ~ D1+D2, dat = sim_data, wts = ipw_tvt_wts_t4$ate_stab) %>% extract_info_msm()
    ipw_tvt_fit3 <- ipw_msm_fit(msm_mod_form = Y3 ~ D1+D2+D3, dat = sim_data, wts = ipw_tvt_wts_t4$ate_stab) %>% extract_info_msm()
    ipw_tvt_fit4 <- ipw_msm_fit(msm_mod_form = Y4 ~ D1+D2+D3+D4, dat = sim_data, wts = ipw_tvt_wts_t4$ate_stab) %>% extract_info_msm()

    ## Using a separate set of weights for each time point
    ipw_tvt_sepweights_list <- list(t1 = ipw_tvt_wts_t1, t2 = ipw_tvt_wts_t2, t3 = ipw_tvt_wts_t3, t4 = ipw_tvt_wts_t4)
    ipw_tvt_fit1_sepweights <- ipw_msm_fit(msm_mod_form = Y1 ~ D1, dat = sim_data, wts = ipw_tvt_wts_t1$ate_stab) %>% extract_info_msm()
    ipw_tvt_fit2_sepweights <- ipw_msm_fit(msm_mod_form = Y2 ~ D1+D2, dat = sim_data, wts = ipw_tvt_wts_t2$ate_stab) %>% extract_info_msm()
    ipw_tvt_fit3_sepweights <- ipw_msm_fit(msm_mod_form = Y3 ~ D1+D2+D3, dat = sim_data, wts = ipw_tvt_wts_t3$ate_stab) %>% extract_info_msm()

    ## Using IPW TVT framework to estimate group-time ATTs (as in did package)
    estimate_gtatt <- function(dat, g, t, ctrl_g, which_wts = c("one_set", "separate"), change_score = FALSE) {
        if (change_score) {
            dat$change_score <- dat[[paste0("Y", t)]] - dat[[paste0("Y", g-1)]]
            mod_form <- paste0("change_score ~ D", g)
        } else {
            mod_form <- paste0("Y", t, " ~ D", g)
        }
        mod_form <- as.formula(mod_form)
        if (which_wts=="one_set") {
            wts <- ipw_tvt_wts_t4$gtatt %>% filter(group==g, time==t, ctrl_group==ctrl_g) %>% pull(wts) %>% unlist()
        } else if (which_wts=="separate") {
            wts <- ipw_tvt_sepweights_list[[t]]$gtatt %>% filter(group==g, time==t, ctrl_group==ctrl_g) %>% pull(wts) %>% unlist()
        }
        ipw_msm_fit(msm_mod_form = mod_form, dat = dat, wts = wts) %>% extract_info_msm(gtatt = TRUE)
    }
    ### Parameters for GTATT estimation
    df_gtatt <- tidyr::crossing(
        group = 1:4,
        time = 1:4,
        ctrl_group = c("never", "notyet"),
        which_wts = c("one_set", "separate"),
        change_score = c(FALSE, TRUE)
    ) %>%
        filter(time >= group) %>%
        filter(!(change_score & group==1))

    estimates_gtatt <- lapply(seq_len(nrow(df_gtatt)), function(r) {
        estimate_gtatt(
            dat = sim_data,
            g = df_gtatt$group[r],
            t = df_gtatt$time[r],
            ctrl_g = df_gtatt$ctrl_group[r],
            which_wts = df_gtatt$which_wts[r],
            change_score = df_gtatt$change_score[r]
        )
    })
    df_gtatt$estimate <- sapply(estimates_gtatt, "[[", "estimate")
    df_gtatt$se <- sapply(estimates_gtatt, "[[", "se")
    
    ## Effect estimation using Callaway and Sant'Anna 2021 (did package)
    ## Reformat data for compatibility with the att_gt() function
    sim_data_long_did <- sim_data %>%
        pivot_longer(-c(id, group), names_to = c(".value", "time"), names_pattern = "(.)(.)") %>%
        mutate(time = as.numeric(time))

    x_formulas <- c(~1, ~X)
    gt_att_res <- lapply(x_formulas, function(x_form) {
        gt_atts_never_dr <- att_gt(yname = "Y", tname = "time", idname = "id", gname = "group", xformla = x_form, panel = TRUE, control_group = "nevertreated", est_method = "dr", data = sim_data_long_did) %>% extract_info_att_gt()
        gt_atts_never_ipw <- att_gt(yname = "Y", tname = "time", idname = "id", gname = "group", xformla = x_form, panel = TRUE, control_group = "nevertreated", est_method = "ipw", data = sim_data_long_did) %>% extract_info_att_gt()
        gt_atts_notyet_dr <- att_gt(yname = "Y", tname = "time", idname = "id", gname = "group", xformla = x_form, panel = TRUE, control_group = "notyettreated", est_method = "dr", data = sim_data_long_did) %>% extract_info_att_gt()
        gt_atts_notyet_ipw <- att_gt(yname = "Y", tname = "time", idname = "id", gname = "group", xformla = x_form, panel = TRUE, control_group = "notyettreated", est_method = "ipw", data = sim_data_long_did) %>% extract_info_att_gt()
        list(never_dr = gt_atts_never_dr, never_ipw = gt_atts_never_ipw, notyet_dr = gt_atts_notyet_dr, notyet_ipw = gt_atts_notyet_ipw)
    })
    names(gt_att_res) <- c("unadj", "adj")

    list(
        ipw_tvt_results = list(time1 = ipw_tvt_fit1, time2 = ipw_tvt_fit2, time3 = ipw_tvt_fit3, time4 = ipw_tvt_fit4, time1_alt = ipw_tvt_fit1_sepweights, time2_alt = ipw_tvt_fit2_sepweights, time3_alt = ipw_tvt_fit3_sepweights, gt_att = df_gtatt),
        gt_att_results = gt_att_res
    )
}

define_struct_equations <- function(multiplier) {
    m <- multiplier
    treat_conf_strength <- 0.5
    struct_eqns1 <- list(
        "U" = tibble(var = "U", when_rel = -1, coeff = 0.5),
        "X" = tibble(var = c("U", "D"), when_rel = c(0, -1), coeff = c(0.5, treat_conf_strength)),
        "D" = tibble(var = c("1", "X", "X"), when_rel = c(0, 0, -1), coeff = c(-1, -0.5, -0.5)*c(1, m, m)),
        "Y" = tibble(var = c("U", "D", "Y"), when_rel = c(0, 0, -1), coeff = c(0.5, 0.5, 0.5)*c(m, 1, 1))
    )

    # struct_eqns2 <- struct_eqns1
    # struct_eqns2$U <- tibble(var = c("U", "D"), when_rel = c(-1, -1), coeff = c(0.5, treat_conf_strength))
    # struct_eqns2$X <- tibble(var = "U", when_rel = 0, coeff = 0.5)

    # struct_eqns3 <- struct_eqns1
    # struct_eqns3$U <- NULL
    # struct_eqns3$X <- tibble(var = "D", when_rel = -1, coeff = treat_conf_strength)
    # struct_eqns3$Y$var <- c("X", "D", "Y")

    struct_eqns4 <- struct_eqns1
    struct_eqns4$W <- struct_eqns4$U
    struct_eqns4$W$var <- "W"
    struct_eqns4$D <- bind_rows(struct_eqns4$D, tibble(var = "W", when_rel = 0, coeff = -0.5*m))
    struct_eqns4$Y <- bind_rows(struct_eqns4$Y, tibble(var = "W", when_rel = 0, coeff = 0.5*m))
    struct_eqns4 <- struct_eqns4[c("U", "W", "X", "D", "Y")]

    # list(setup1 = struct_eqns1, setup2 = struct_eqns2, setup3 = struct_eqns3, setup4 = struct_eqns4)
    list(setup1 = struct_eqns1, setup4 = struct_eqns4)
}

get_n_all_strats <- function(dat, num_time_periods) {
    which_vars <- paste0("D", seq_len(num_time_periods))
    dat %>%
        count(across({{ which_vars }}))
}

run_simulations <- function(struct_eqns_list, n, num_time_periods) {
    ## Simulate data and estimate causal effects
    sim_data_list <- lapply(struct_eqns_list, simulate_data_from_equations, n = n, num_tp = num_time_periods, absorbing = TRUE, under_intervention = FALSE, treat_strategy = NULL)
    estimates <- lapply(sim_data_list, estimate_effects)
    samp_sizes_by_strat <- lapply(sim_data_list, get_n_all_strats, num_time_periods = num_time_periods)
    list(estimates = estimates, samp_sizes = samp_sizes_by_strat)
}

set.seed(141)
sim_results <- lapply(multipliers, function(multiplier) {
    ## Define structural equations for different data settings
    struct_eqns_list <- define_struct_equations(multiplier)

    ## Repeat simulations
    estimates_list <- replicate(n = num_reps, {
        list_over_sampsizes <- lapply(samp_sizes, function(this_n) {
            run_simulations(struct_eqns_list = struct_eqns_list, n = this_n, num_time_periods = num_time_periods)
        })
        names(list_over_sampsizes) <- samp_sizes
        list_over_sampsizes
    }, simplify = FALSE)

    ## "Simulate" expected values of variables under interventions
    truth_list <- lapply(struct_eqns_list, get_truth_all_strats, num_tp = num_time_periods, absorbing = TRUE)

    list(estimates = estimates_list, truth = truth_list)
})

write_rds(sim_results, file = "../data/sim_results.rds", compress = "gz")

## ============================================================================
## Session Info
sessionInfo()

