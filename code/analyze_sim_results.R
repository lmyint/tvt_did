library(tidyverse)
library(did)
library(gridExtra)
library(ggthemes)

all_sim_results <- read_rds("../data/sim_results.rds")

compute_truth_for_estimands <- function(df_means_by_strat_all_vars) {
    df <- df_means_by_strat_all_vars
    df <- df %>%
        select(starts_with("D"), starts_with("Y")) %>%
        unique() %>%
        mutate(group = 5-(D1+D2+D3+D4), group = ifelse(group==5, 0, group))
    
    get_true_gtatt <- function(g, t) {
        y_var <- paste0("Y", t)
        untreated_po <- df %>% filter(group==0) %>% pull({{ y_var }})
        treated_po <- df %>% filter(group==g) %>% pull({{ y_var }})
        treated_po - untreated_po
    }
    df_truth_gt <- tibble(
        group = c(rep(1, 4), rep(2, 3), rep(3, 2), 4),
        time = c(1:4, 2:4, 3:4, 4)
    )
    df_truth_gt$truth <- sapply(seq_len(nrow(df_truth_gt)), function(r) {
        get_true_gtatt(df_truth_gt$group[r], df_truth_gt$time[r])
    })

    get_true_msm_coeff <- function(coeff, t) {
        # D1 = group1 - group2
        # D2 = group2 - group3
        # D3 = group3 - group4
        # D4 = group4 - group0
        # intcpt = group0
        y_var <- paste0("Y", t)
        if (coeff=="Intercept") {
            term1 <- df %>% filter(group==0) %>% pull({{ y_var }})
            term2 <- 0
        } else if (coeff=="D1") {
            term1 <- df %>% filter(group==1) %>% pull({{ y_var }})
            term2 <- df %>% filter(group==2) %>% pull({{ y_var }})
        } else if (coeff=="D2") {
            term1 <- df %>% filter(group==2) %>% pull({{ y_var }})
            term2 <- df %>% filter(group==3) %>% pull({{ y_var }})
        } else if (coeff=="D3") {
            term1 <- df %>% filter(group==3) %>% pull({{ y_var }})
            term2 <- df %>% filter(group==4) %>% pull({{ y_var }})
        } else if (coeff=="D4") {
            term1 <- df %>% filter(group==4) %>% pull({{ y_var }})
            term2 <- df %>% filter(group==0) %>% pull({{ y_var }})
        }
        term1 - term2
    }

    df_truth_msm <- tibble(
        coeff = c(rep("Intercept", 4), rep("D1", 4), rep("D2", 3), rep("D3", 2), "D4"),
        time = c(1:4, 1:4, 2:4, 3:4, 4)
    )
    df_truth_msm$truth <- sapply(seq_len(nrow(df_truth_msm)), function(r) {
        get_true_msm_coeff(df_truth_msm$coeff[r], df_truth_msm$time[r])
    })

    list(gtatt = df_truth_gt, msm = df_truth_msm)
}

get_results_for_one_mult <- function(sim_results, what = c("bias", "n")) {
    data_setup_names <- names(sim_results$truth)
    list_over_setups <- lapply(data_setup_names, function(data_setup_name) {
        get_results_for_one_data_setup(sim_results, data_setup_name, what)
    })
    names(list_over_setups) <- data_setup_names
    bind_rows(list_over_setups, .id = "data_setup")
}

get_results_for_one_data_setup <- function(sim_results, data_setup_name, what) {
    truth <- compute_truth_for_estimands(sim_results$truth[[data_setup_name]])

    list_over_reps <- lapply(sim_results$estimates, get_results_for_one_sim_rep, data_setup_name = data_setup_name, df_truth_gt = truth$gtatt, df_truth_msm = truth$msm, what = what)
    names(list_over_reps) <- seq_along(list_over_reps)
    bind_rows(list_over_reps, .id = "sim_rep") %>% mutate(sim_rep = as.integer(sim_rep))
}

get_results_for_one_sim_rep <- function(sim_results, data_setup_name, df_truth_gt, df_truth_msm, what) {
    list_over_sampsizes <- lapply(sim_results, get_results_for_one_sampsize, data_setup_name = data_setup_name, df_truth_gt = df_truth_gt, df_truth_msm = df_truth_msm, what = what)
    bind_rows(list_over_sampsizes, .id = "samp_size") %>% mutate(samp_size = as.integer(samp_size))
}

get_results_for_one_sampsize <- function(sim_results, data_setup_name, df_truth_gt, df_truth_msm, what) {
    if (what=="n") {
        df_n <- sim_results$samp_sizes[[data_setup_name]] %>%
            mutate(group = D1+D2+D3+D4, n_frac = n/sum(n)) %>%
            select(group, n_frac) %>%
            pivot_wider(names_from = group, names_prefix = "n_frac_g", values_from = n_frac)
        return(df_n)
    }
    sim_results <- sim_results$estimates[[data_setup_name]]

    ## Results for TVT framework
    res_ipw_tvt <- sim_results$ipw_tvt_results
    bool_msm <- names(res_ipw_tvt) != "gt_att"
    res_ipw_tvt_msm <- res_ipw_tvt[bool_msm]
    res_ipw_tvt_msm <- bind_rows(res_ipw_tvt_msm, .id = "method") %>%
        separate(col = "method", into = c("time", "method"), sep = "_", fill = "right") %>%
        mutate(
            method = ifelse(is.na(method), "one_wts", "mult_wts"),
            coeff = str_remove_all(coeff, "TRUE|\\(|\\)"),
            time = as.integer(str_remove(time, "time"))
        )
    res_ipw_tvt_msm <- res_ipw_tvt_msm %>%
        left_join(df_truth_msm, by = c("coeff", "time")) %>%
        mutate(
            bias = estimate - truth,
            abs_perc_bias = 100*abs(bias)/abs(truth),
            estimand = paste0(coeff, "_t", time),
            method = paste0("ipw_tvt_", method)
        ) %>%
        select(method, estimand, bias, abs_perc_bias, se)
    
    res_ipw_tvt_gtatt <- res_ipw_tvt$gt_att %>%
        left_join(df_truth_gt, by = c("group", "time")) %>%
        mutate(
            bias = estimate - truth,
            abs_perc_bias = 100*abs(bias)/abs(truth),
            estimand = paste0("g", group, "_t", time),
            method = paste0(ctrl_group, "_", which_wts, "_", ifelse(change_score, "cs", "nocs"))
        ) %>% 
        select(method, estimand, bias, abs_perc_bias, se)

    ## Results for Callaway & Sant'Anna framework
    res_gt_att <- sim_results$gt_att_results
    res_gt_att <- lapply(res_gt_att, function(list_by_adj) {
        bind_rows(list_by_adj, .id = "ctrl_grp_est_meth")
    }) %>% bind_rows(.id = "adj_type")
    res_gt_att <- res_gt_att %>%
        filter(group != "dynamic") %>%
        mutate(group = as.integer(group)) %>%
        left_join(df_truth_gt, by = c("group", "time")) %>%
        filter(!is.na(truth)) %>%
        mutate(
            bias = estimate - truth,
            abs_perc_bias = 100*abs(bias)/abs(truth),
            method = paste0(adj_type, "_", ctrl_grp_est_meth),
            estimand = paste0("g", group, "_t", time)
        ) %>% 
        select(method, estimand, bias, abs_perc_bias, se)
    
    bind_rows(res_ipw_tvt_msm, res_ipw_tvt_gtatt, res_gt_att)
}

all_bias_results_list <- lapply(all_sim_results, get_results_for_one_mult, what = "bias")
names(all_bias_results_list) <- 1:3
all_bias_results <- bind_rows(all_bias_results_list, .id = "multiplier") %>%
    mutate(multiplier = as.numeric(multiplier))

samp_size_list <- lapply(all_sim_results, get_results_for_one_mult, what = "n")
names(samp_size_list) <- 1:3
samp_sizes <- bind_rows(samp_size_list, .id = "multiplier") %>%
    mutate(multiplier = as.numeric(multiplier))
samp_sizes %>%
    group_by(data_setup, multiplier, samp_size) %>%
    summarize(across(starts_with("n_frac"), mean)) %>%
    as.data.frame()
samp_sizes %>%
    group_by(data_setup, multiplier, samp_size) %>%
    summarize(across(starts_with("n_frac"), mean)) %>%
    filter(samp_size==3000) %>%
    select(-samp_size) %>%
    as.data.frame()
samp_sizes %>%
    pivot_longer(starts_with("n_frac"), names_to = "treat_strat", values_to = "n_frac") %>%
    mutate(treat_strat = str_remove(treat_strat, "n_frac_")) %>%
    ggplot(aes(x = data_setup, y = n_frac)) +
        geom_boxplot() +
        facet_grid(multiplier ~ treat_strat) +
        theme_classic()

## Exploring results
pdf("../results/explore_bias_results.pdf", width = 2*12, height = 2*4)
for (m in 1:3) {
    all_bias_results_msm <- all_bias_results %>%
        filter(multiplier==m, !str_detect(estimand, "^g"))
    all_bias_results_gtatt <- all_bias_results %>%
        filter(multiplier==m, str_detect(estimand, "^g")) %>%
        mutate(cs_or_notyet = str_detect(method, "_cs$") | str_detect(method, "_notyet_"))
    method_colors <- c(rep(c("lightgreen", "seagreen", "lightblue1", "dodgerblue"), 2), rep(c("lightpink", "deeppink"), each = 4))
    names(method_colors) <- all_bias_results_gtatt$method %>% unique()
    p1 <- ggplot(all_bias_results_msm, aes(x = estimand, y = bias, color = method)) +
        geom_tufteboxplot(median.type = "point") +
        facet_grid(samp_size ~ data_setup, labeller = label_both) +
        geom_hline(yintercept = 0, color = "red") +
        theme_classic() +
        scale_color_manual(values = c("ipw_tvt_one_wts" = "dodgerblue", "ipw_tvt_mult_wts" = "deepskyblue")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(title = paste("Multiplier:", m))
    p2 <- ggplot(all_bias_results_gtatt, aes(x = estimand, y = bias, color = method, linetype = cs_or_notyet)) +
        geom_tufteboxplot(median.type = "point") +
        facet_grid(samp_size ~ data_setup, labeller = label_both) +
        geom_hline(yintercept = 0, color = "red") +
        theme_classic() +
        scale_color_manual(values = method_colors) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(title = paste("Multiplier:", m))
    print(p1)
    print(p2)
}
dev.off()

## ======================================================================
## Exploring relationship between bias and coefficient multiplier (confounding strength)
pdf("../results/explore_bias_results_by_conf_strength.pdf", width = 4*2, 24*2)
all_bias_results_subs <- all_bias_results %>%
    filter(method %in% c("ipw_tvt_one_wts", "notyet_one_set_nocs", "notyet_one_set_cs", "unadj_notyet_ipw", "adj_notyet_ipw"), samp_size==3000)
for (setup in c(1,4)) {
    p <- all_bias_results_subs %>%
        filter(data_setup==paste0("setup", setup)) %>%
        ggplot(aes(x = factor(multiplier), y = bias)) +
            geom_tufteboxplot() +
            facet_grid(estimand ~ method, scales = "free_y") +
            geom_hline(yintercept = 0, color = "red") +
            theme_classic() +
            labs(title = paste("Data setup:", setup))
    print(p)
}
dev.off()

## ======================================================================
## Figures for paper
## Simplify presentation by only showing n = 3000; one set of weights (versus
## separate sets for each time point); IPW only for `did` package (versus also
## showing DR estimator)
methods_ordered <- c("ipw_tvt_one_wts", "never_one_set_nocs", "notyet_one_set_nocs", "never_one_set_cs", "notyet_one_set_cs", "unadj_never_ipw", "unadj_notyet_ipw", "adj_never_ipw", "adj_notyet_ipw")
estimands_ordered <- c(paste0("Intercept_t", 1:4), paste0("D1_t", 1:4), paste0("D2_t", 2:4), paste0("D3_t", 3:4), "D4_t4", paste0("g1_t", 1:4), paste0("g2_t", 2:4), paste0("g3_t", 3:4), "g4_t4")
df_estimands <- tibble(estimand = estimands_ordered) %>%
    mutate(
        estimand_neat = str_replace_all(estimand, "_t", "\nt="),
        estimand_neat = str_replace_all(estimand_neat, "^g", "g="),
        estimand_neat = str_replace_all(estimand_neat, "Intercept", "Int")
    )
all_bias_results_subs <- all_bias_results %>%
    filter(samp_size==3000, method != "ipw_tvt_mult_wts", !str_detect(method, "_separate_|_dr$")) %>%
    mutate(
        estimand_type = case_when(
            str_detect(estimand, "^D[0-9]_t[0-9]") ~ "TVT",
            str_detect(estimand, "^Intercept") ~ "TVT",
            str_detect(estimand, "^g[0-9]_t[0-9]") ~ "TVT_GTATT"
        ),
        data_setup = case_when(
            data_setup=="setup1" ~ "Seq. exch. holds",
            data_setup=="setup4" ~ "Seq. exch. does not hold"
        ),
        data_setup = factor(data_setup, levels = c("Seq. exch. holds", "Seq. exch. does not hold")),
        method = factor(method, levels = methods_ordered)
    ) %>%
    left_join(df_estimands) %>%
    mutate(estimand_neat = factor(estimand_neat, levels = df_estimands$estimand_neat))
method_colors <- c("dodgerblue", "lightblue1", "dodgerblue", "lightgreen", "seagreen", "#f28500", "#ffbb00", "firebrick4", "firebrick1")
names(method_colors) <- unique(all_bias_results_subs$method)

plot_list <- list()
index <- 1
for (what in c("bias", "se")) {
    if (what=="bias") {
        ylab <- "Bias"
    } else if (what=="se") {
        ylab <- "SE"
    }
    for (est_type in c("TVT", "TVT_GTATT")) {
        if (est_type=="TVT") {
            title <- "Estimands specific to TVT framework (MSM coefficients)"
            xlab <- "Coefficient"
        } else if (est_type=="TVT_GTATT") {
            title <- "Group-time average treatment effects"
            xlab <- "ATT(g,t)"
        }
        plot_data <- all_bias_results_subs %>%
            filter(multiplier==1, estimand_type==est_type)
        if (est_type=="TVT_GTATT") {
            plot_data <- plot_data %>%
                mutate(
                    is_tvt_method = !str_detect(method, "adj"),
                    estimand_neat = paste0(estimand_neat, ifelse(is_tvt_method, " ", "  "))
                )
        }
        p <- ggplot(plot_data, aes(x = estimand_neat, y = .data[[what]], color = method)) +
            geom_tufteboxplot() +
            theme_classic() +
            scale_color_manual(
                values = method_colors, 
                breaks = unique(all_bias_results_subs$method),
                labels = c(
                    "TVT: IPW MSM",
                    "TVT: not change score\nControls: never treated",
                    "TVT: not change score\nControls: not yet treated",
                    "TVT: change score\nControls: never treated",
                    "TVT: change score\nControls: not yet treated",
                    "CS2021: unadjusted\nControls: never treated",
                    "CS2021: unadjusted\nControls: not yet treated",
                    "CS2021: adjusted\nControls: never treated",
                    "CS2021: adjusted\nControls: not yet treated"
                )
            ) +
            facet_grid(. ~ data_setup) +
            labs(x = xlab, y = ylab, title = title) +
            theme(text = element_text(size = 14))
        if (est_type=="TVT") {
            p <- p +
                geom_vline(xintercept = c(4,8,11,13)+0.5, color = "black", lty = "dashed")
        } else if (est_type=="TVT_GTATT") {
            p <- p + geom_vline(xintercept = c(4,10,14)+0.5, color = "black", lty = "dashed")
        }
        if (what=="bias") {
            p <- p + geom_hline(yintercept = 0, color = "red")
        }
        plot_list[[index]] <- p
        index <- index + 1
    }
}

pdf("../results/results_main.pdf", width = 24, height = 8)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], nrow = 2, ncol = 2)
dev.off()

plot_list <- list()
index <- 1
plot_data <- all_bias_results_subs %>%
    filter(estimand=="g2_t2")
for (what in c("bias", "se")) {
    if (what=="bias") {
        ylab <- "Bias"
    } else if (what=="se") {
        ylab <- "SE"
    }
    p <- ggplot(plot_data, aes(x = factor(multiplier), y = .data[[what]], color = method)) +
        geom_tufteboxplot() +
        theme_classic() +
        scale_color_manual(
            values = method_colors, 
            breaks = unique(all_bias_results_subs$method),
            labels = c(
                "TVT: IPW MSM",
                "TVT: not change score\nControls: never treated",
                "TVT: not change score\nControls: not yet treated",
                "TVT: change score\nControls: never treated",
                "TVT: change score\nControls: not yet treated",
                "CS2021: unadjusted\nControls: never treated",
                "CS2021: unadjusted\nControls: not yet treated",
                "CS2021: adjusted\nControls: never treated",
                "CS2021: adjusted\nControls: not yet treated"
            )
        ) +
        facet_grid(. ~ data_setup) +
        labs(x = "Confounding strength multiplier", y = ylab) +
        theme(text = element_text(size = 14))
    if (what=="bias") {
        p <- p + geom_hline(yintercept = 0, color = "red")
    }
    plot_list[[index]] <- p
    index <- index + 1
}

pdf("../results/results_conf_strength_mult.pdf", width = 14, height = 6)
grid.arrange(plot_list[[1]], plot_list[[2]], ncol = 2)
dev.off()
