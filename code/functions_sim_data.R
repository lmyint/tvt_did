#' Generate a logical vector from a log-odds vector
#' @param x Vector of log odds
quant_to_bin <- function(x) {
    n <- length(x)
    odds <- exp(x)
    p <- odds/(1+odds)
    as.logical(rbinom(n, size = 1, prob = p))
}

#' Simulate data from a pre-specified set of structural equations for a 
#' user-defined number of time periods. Treatment variables are binary.
#' @param n Number of cases
#' @param num_tp Number of time periods
#' @param treat_conf_strength Strength of treatment-confounder feedback
#'                            (Coefficient on D_t-1 in the structural 
#'                             equation for X_t)
#' @param absorbing If TRUE, once treated, a unit stays treated

simulate_data_from_equations <- function(struct_eqns, n = NULL, num_tp = 4, absorbing = TRUE, under_intervention = FALSE, treat_strategy = NULL) {
    stopifnot(num_tp >= 2)

    ## Create tibble with correct number of rows using a placeholder variable
    if (under_intervention)
        sim_data <- tibble(temp_var = 0)
    else
        sim_data <- tibble(temp_var = rep(0, n))

    ## Loop over time points and variables
    for (tp in seq_len(num_tp)) {
        for (i in seq_along(struct_eqns)) {
            this_var <- names(struct_eqns)[i] # e.g., U
            this_var_this_time <- paste0(this_var, tp) # e.g., U2

            design_beta <- create_design_mat_beta(struct_eqns[[i]], tp, sim_data)
            design <- design_beta$design
            beta <- design_beta$beta

            noise <- create_noise(n, under_intervention)
            lin_comb <- as.numeric(design %*% beta)
            sim_data[[this_var_this_time]] <- lin_comb + noise

            if (this_var=="D") { ## Treatment variable
                if (under_intervention) {
                    sim_data[[this_var_this_time]] <- treat_strategy[tp]
                } else {
                    sim_data[[this_var_this_time]] <- quant_to_bin(lin_comb)
                    if (absorbing & tp >= 2) {
                        node_D_prev <- paste0("D", tp-1)
                        sim_data[[this_var_this_time]] <- ifelse(sim_data[[node_D_prev]], TRUE, sim_data[[this_var_this_time]])
                    }
                }
            }
        }
    }
    sim_data %>% select(-temp_var)
}

create_noise <- function(n, under_intervention) {
    if (under_intervention) {
        noise <- 0
    } else {
        noise <- rnorm(n, 0, 1)
    }
    noise
}

create_design_mat_beta <- function(df_eqns, time_point, dat) {
    ## Determine if intercept is zero or not
    if ("1" %in% df_eqns$var)
        intercept <- df_eqns %>% filter(var=="1") %>% pull(coeff)
    else
        intercept <- 0

    df_eqns_new <- df_eqns %>%
        filter(var != "1") %>%
        mutate(
            when_abs = when_rel + time_point,
            var = paste0(var, when_abs)
        ) %>%
        filter(when_abs >= 1)
    if (nrow(df_eqns_new)==0) {
        design <- model.matrix(~1, data = dat)
        beta <- matrix(intercept, nrow = 1, ncol = 1)
    } else {
        design_mod_form <- paste("~", paste(df_eqns_new$var, collapse = "+"))
        design <- model.matrix(as.formula(design_mod_form), data = dat)
        beta <- c(intercept, df_eqns_new$coeff)
    }
    list(design = design, beta = beta)
}

get_truth_all_strats <- function(struct_eqns, num_tp, absorbing = TRUE) {
    args_list <- rep(list(c(0,1)), num_tp)
    names(args_list) <- paste0("D", seq_len(num_tp))
    all_strategies <- do.call(tidyr::crossing, args_list)
    all_strategies <- as.matrix(all_strategies)
    all_strategies <- all_strategies[,names(args_list)]
    if (absorbing) {
        is_absorbing <- sapply(seq_len(nrow(all_strategies)), function(r) {
            !is.unsorted(all_strategies[r,])
        })
        all_strategies <- all_strategies[is_absorbing,,drop=FALSE]
    }

    truth_list <- lapply(seq_len(nrow(all_strategies)), function(r) {
        strat <- all_strategies[r,]
        simulate_data_from_equations(struct_eqns = struct_eqns, n = NULL, num_tp = num_tp, absorbing = absorbing, under_intervention = TRUE, treat_strategy = strat)
    })
    bind_rows(truth_list)
}
