get_static_strategies <- function(num_tp, absorbing = TRUE) {
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
    all_strategies
}

#' Generate a logical vector from a log-odds vector
#' @param x Vector of log odds
quant_to_bin <- function(x) {
    n <- length(x)
    odds <- exp(x)
    p <- odds/(1+odds)
    as.logical(rbinom(n, size = 1, prob = p))
}

simulate_all_variables <- function(fn_list, num_tp, N = NULL, sim_data_pop = NULL, absorbing = TRUE, multiplier, treat_strategy = NULL) {
    if (is.null(sim_data_pop)) {
        sim_data_pop <- tibble(temp_var = rep(0, N))
    } else {
        N <- nrow(sim_data_pop)
    }
    zeroes <- rep(0, N)

    ## Create uber-list of variables: variables at current and previous time, noise
    ## For use in supplied functions in fn_list
    var_list <- vector("list", 16)
    names(var_list) <- c(
        paste0(c("U", "W", "X", "D", "Y"), "_prev"),
        paste0(c("U", "W", "X", "D", "Y"), "_curr"),
        paste0("noise_", c("U", "W", "X", "Y")),
        "noise_X_invariant", "noise_U_invariant"
    )
    var_list$noise_X_invariant <- rnorm(N, 0, 1)
    var_list$noise_U_invariant <- rnorm(N, 0, 1)
    for (tp in seq_len(num_tp)) {
        if (tp==1) {
            var_list$U_prev <- zeroes
            var_list$W_prev <- zeroes
            var_list$X_prev <- zeroes
            var_list$D_prev <- zeroes
            var_list$Y_prev <- zeroes
        } else {
            var_list$U_prev <- sim_data_pop[[paste0("U", tp-1)]]
            var_list$W_prev <- sim_data_pop[[paste0("W", tp-1)]]
            var_list$X_prev <- sim_data_pop[[paste0("X", tp-1)]]
            var_list$D_prev <- sim_data_pop[[paste0("D", tp-1)]]
            var_list$Y_prev <- sim_data_pop[[paste0("Y", tp-1)]]
        }
        if (is.null(treat_strategy)) {
            var_list$noise_U <- rnorm(N, 0, 1)
            var_list$noise_W <- rnorm(N, 0, 1)
            var_list$noise_X <- rnorm(N, 0, 1)
            var_list$noise_Y <- rnorm(N, 0, 1)
        } else {
            var_list$noise_U <- sim_data_pop[[paste0("e_U", tp)]]
            var_list$noise_W <- sim_data_pop[[paste0("e_W", tp)]]
            var_list$noise_X <- sim_data_pop[[paste0("e_X", tp)]]
            var_list$noise_Y <- sim_data_pop[[paste0("e_Y", tp)]]
        }

        for (i in seq_along(fn_list)) {
            var_sim_fn <- fn_list[[i]] # e.g., simulate_U()
            the_var <- names(fn_list)[i] # e.g., "U"
            the_var_curr <- paste0(the_var, "_curr") # e.g., "U_curr"
            the_var_this_time <- paste0(the_var, tp) # e.g., "U1"
            the_var_err_this_time <- paste0("e_", the_var, tp) # e.g., "e_U1"
            the_var_noise <- paste0("noise_", the_var)

            if (the_var_this_time %in% colnames(sim_data_pop)) {
                ## Exogenous variable and simulating under intervention
                var_list[[the_var_curr]] <- sim_data_pop[[the_var_this_time]]
                next
            }

            if (the_var=="D") {
                if (is.null(treat_strategy)) { # Natural world
                    var_list[[the_var_curr]] <- var_sim_fn(var_list, absorbing = absorbing, mult = multiplier)
                    sim_data_pop[[the_var_this_time]] <- var_list[[the_var_curr]]
                } else { # Under intervention
                    var_list[[the_var_curr]] <- treat_strategy[tp]
                    sim_data_pop[[the_var_this_time]] <- treat_strategy[tp]
                }
            } else if (the_var=="Y") {
                var_list[[the_var_curr]] <- var_sim_fn(var_list, mult = multiplier)
                sim_data_pop[[the_var_this_time]] <- var_list[[the_var_curr]]
                sim_data_pop[[the_var_err_this_time]] <- var_list[[the_var_noise]]
            } else {
                var_list[[the_var_curr]] <- var_sim_fn(var_list)
                sim_data_pop[[the_var_this_time]] <- var_list[[the_var_curr]]
                sim_data_pop[[the_var_err_this_time]] <- var_list[[the_var_noise]]
            }
        }
    }

    sim_data_pop
}

###########################################################################
## Data setup 1: Time-invariant measured confounder
###########################################################################
simulate_X_setup1 <- function(var_list) {
    var_list$noise_X_invariant
}

simulate_D_setup1 <- function(var_list, absorbing, mult) {
    X_curr <- var_list$X_curr
    D_prev <- var_list$D_prev
    log_odds_D <- -1 - 0.5*mult*X_curr
    D_curr <- quant_to_bin(log_odds_D)
    if (absorbing) {
        D_curr <- ifelse(D_prev, TRUE, D_curr)
    }
    D_curr
}

simulate_Y_setup1 <- function(var_list, mult) {
    X_curr <- var_list$X_curr
    D_curr <- var_list$D_curr
    noise_Y <- var_list$noise_Y
    0.5*mult*X_curr + 0.5*D_curr + 0.1*X_curr*D_curr + noise_Y
}

#' @param N Population size
#' @param num_tp Number of time periods
#' @param multiplier Confounding strength multiplier
#' @param absorbing Is being treated an absorbing state?
simulate_data_setup1 <- function(N = 1e6, num_tp, multiplier, absorbing = TRUE) {
    fn_list <- list("X" = simulate_X_setup1, "D" = simulate_D_setup1, "Y" = simulate_Y_setup1)

    ## Simulate variables under NATURAL circumstances
    sim_data_pop <- simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = NULL, absorbing = absorbing, multiplier = multiplier, treat_strategy = NULL)

    ## Identify all exogenous variables
    exo_vars <- c(paste0("X", 1:4), colnames(sim_data_pop)[str_detect(colnames(sim_data_pop), "^e_")])

    ## Simulate variables under INTERVENTION
    treat_strategies <- get_static_strategies(num_tp = num_tp, absorbing = absorbing)
    sim_data_pop_intervention <- lapply(seq_len(nrow(treat_strategies)), function(i) {
        simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = sim_data_pop[,exo_vars], absorbing = absorbing, multiplier = multiplier, treat_strategy = treat_strategies[i,])
    })
    treat_strats_char <- sapply(seq_len(nrow(treat_strategies)), function(i) {
        paste(treat_strategies[i,], collapse = "")
    })
    names(sim_data_pop_intervention) <- treat_strats_char

    list(
        pop_natural = sim_data_pop,
        pop_intervention = sim_data_pop_intervention
    )
}

###########################################################################
## Data setup 2: Time-invariant unmeasured confounder
###########################################################################
simulate_U_setup2 <- function(var_list) {
    var_list$noise_U_invariant
}

simulate_D_setup2 <- function(var_list, absorbing, mult) {
    X_curr <- var_list$X_curr
    U_curr <- var_list$U_curr
    D_prev <- var_list$D_prev
    log_odds_D <- -1 - 0.5*mult*X_curr - 0.5*mult*U_curr
    D_curr <- quant_to_bin(log_odds_D)
    if (absorbing) {
        D_curr <- ifelse(D_prev, TRUE, D_curr)
    }
    D_curr
}

simulate_Y_setup2 <- function(var_list, mult) {
    X_curr <- var_list$X_curr
    U_curr <- var_list$U_curr
    D_curr <- var_list$D_curr
    noise_Y <- var_list$noise_Y
    0.5*mult*X_curr + 0.5*mult*U_curr + 0.5*D_curr + 0.1*X_curr*D_curr + noise_Y
}

#' @param N Population size
#' @param num_tp Number of time periods
#' @param multiplier Confounding strength multiplier
#' @param absorbing Is being treated an absorbing state?
simulate_data_setup2 <- function(N = 1e6, num_tp, multiplier, absorbing = TRUE) {
    fn_list <- list("X" = simulate_X_setup1, "U" = simulate_U_setup2, "D" = simulate_D_setup2, "Y" = simulate_Y_setup2)

    ## Simulate variables under NATURAL circumstances
    sim_data_pop <- simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = NULL, absorbing = absorbing, multiplier = multiplier, treat_strategy = NULL)

    ## Identify all exogenous variables
    exo_vars <- c(paste0("X", 1:4), paste0("U", 1:4), colnames(sim_data_pop)[str_detect(colnames(sim_data_pop), "^e_")])

    ## Simulate variables under INTERVENTION
    treat_strategies <- get_static_strategies(num_tp = num_tp, absorbing = absorbing)
    sim_data_pop_intervention <- lapply(seq_len(nrow(treat_strategies)), function(i) {
        simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = sim_data_pop[,exo_vars], absorbing = absorbing, multiplier = multiplier, treat_strategy = treat_strategies[i,])
    })
    treat_strats_char <- sapply(seq_len(nrow(treat_strategies)), function(i) {
        paste(treat_strategies[i,], collapse = "")
    })
    names(sim_data_pop_intervention) <- treat_strats_char

    list(
        pop_natural = sim_data_pop,
        pop_intervention = sim_data_pop_intervention
    )
}

###########################################################################
## Data setup 3: Time-varying confounder, measured proxy (X)
###########################################################################
simulate_U_setup3 <- function(var_list) {
    U_prev <- var_list$U_prev
    noise_U <- var_list$noise_U
    0.5*U_prev + noise_U
}

simulate_X_setup3 <- function(var_list) {
    U_curr <- var_list$U_curr
    D_prev <- var_list$D_prev
    noise_X <- var_list$noise_X
    0.5*U_curr + 0.5*D_prev + noise_X
}

simulate_D_setup3 <- function(var_list, absorbing, mult) {
    X_curr <- var_list$X_curr
    X_prev <- var_list$X_prev
    D_prev <- var_list$D_prev
    log_odds_D <- -1 - 0.5*mult*X_curr - 0.5*mult*X_prev + 0.5*D_prev
    D_curr <- quant_to_bin(log_odds_D)
    if (absorbing) {
        D_curr <- ifelse(D_prev, TRUE, D_curr)
    }
    D_curr
}

simulate_Y_setup3 <- function(var_list, mult) {
    U_curr <- var_list$U_curr
    D_curr <- var_list$D_curr
    Y_prev <- var_list$Y_prev
    noise_Y <- var_list$noise_Y
    0.5*mult*U_curr + 0.5*D_curr + 0.5*Y_prev + 0.1*U_curr*D_curr + noise_Y
}

simulate_data_setup3 <- function(N = 1e6, num_tp, multiplier, absorbing = TRUE) {
    fn_list <- list("U" = simulate_U_setup3, "X" = simulate_X_setup3, "D" = simulate_D_setup3, "Y" = simulate_Y_setup3)

    ## Simulate variables under NATURAL circumstances
    sim_data_pop <- simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = NULL, absorbing = absorbing, multiplier = multiplier, treat_strategy = NULL)

    ## Identify all exogenous variables
    exo_vars <- c("U1", colnames(sim_data_pop)[str_detect(colnames(sim_data_pop), "^e_")])

    ## Simulate variables under INTERVENTION
    treat_strategies <- get_static_strategies(num_tp = num_tp, absorbing = absorbing)
    sim_data_pop_intervention <- lapply(seq_len(nrow(treat_strategies)), function(i) {
        simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = sim_data_pop[,exo_vars], absorbing = absorbing, multiplier = multiplier, treat_strategy = treat_strategies[i,])
    })
    treat_strats_char <- sapply(seq_len(nrow(treat_strategies)), function(i) {
        paste(treat_strategies[i,], collapse = "")
    })
    names(sim_data_pop_intervention) <- treat_strats_char

    list(
        pop_natural = sim_data_pop,
        pop_intervention = sim_data_pop_intervention
    )
}


###########################################################################
## Data setup 4: Time-varying confounder, unmeasured
###########################################################################
simulate_W_setup4 <- function(var_list) {
    W_prev <- var_list$W_prev
    noise_W <- var_list$noise_W
    0.5*W_prev + noise_W
}

simulate_D_setup4 <- function(var_list, absorbing, mult) {
    X_curr <- var_list$X_curr
    X_prev <- var_list$X_prev
    W_curr <- var_list$W_curr
    D_prev <- var_list$D_prev
    log_odds_D <- -1 - 0.5*mult*X_curr - 0.5*mult*X_prev - 0.5*mult*W_curr + 0.5*D_prev
    D_curr <- quant_to_bin(log_odds_D)
    if (absorbing) {
        D_curr <- ifelse(D_prev, TRUE, D_curr)
    }
    D_curr
}

simulate_Y_setup4 <- function(var_list, mult) {
    U_curr <- var_list$U_curr
    W_curr <- var_list$W_curr
    D_curr <- var_list$D_curr
    Y_prev <- var_list$Y_prev
    noise_Y <- var_list$noise_Y
    0.5*mult*U_curr + 0.5*mult*W_curr + 0.5*D_curr + 0.5*Y_prev + 0.1*U_curr*D_curr + noise_Y
}

simulate_data_setup4 <- function(N = 1e6, num_tp, multiplier, absorbing = TRUE) {
    ## U and X functions are the same from setup 3
    fn_list <- list("U" = simulate_U_setup3, "W" = simulate_W_setup4, "X" = simulate_X_setup3, "D" = simulate_D_setup4, "Y" = simulate_Y_setup4)

    ## Simulate variables under NATURAL circumstances
    sim_data_pop <- simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = NULL, absorbing = absorbing, multiplier = multiplier, treat_strategy = NULL)

    ## Identify all exogenous variables
    exo_vars <- c("U1", "W1", colnames(sim_data_pop)[str_detect(colnames(sim_data_pop), "^e_")])

    ## Simulate variables under INTERVENTION
    treat_strategies <- get_static_strategies(num_tp = num_tp, absorbing = absorbing)
    sim_data_pop_intervention <- lapply(seq_len(nrow(treat_strategies)), function(i) {
        simulate_all_variables(fn_list = fn_list, num_tp = num_tp, N = N, sim_data_pop = sim_data_pop[,exo_vars], absorbing = absorbing, multiplier = multiplier, treat_strategy = treat_strategies[i,])
    })
    treat_strats_char <- sapply(seq_len(nrow(treat_strategies)), function(i) {
        paste(treat_strategies[i,], collapse = "")
    })
    names(sim_data_pop_intervention) <- treat_strats_char

    list(
        pop_natural = sim_data_pop,
        pop_intervention = sim_data_pop_intervention
    )
}
