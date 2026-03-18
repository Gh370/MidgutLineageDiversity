
# Robustness analysis of the age-associated clonal diversity decline in ISCs of the adult Drosophila posterior midgut

rm(list = ls())

library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(readr)
library(phangorn)
library(patchwork)

# Load data

load('~/.../data/Phylotime/Phylotime_PosteriorMidgut.RData')

tree_names <- ls(pattern = "^S\\d+_P_\\d+d$")

extract_day <- function(x) {
  as.integer(str_match(x, "_(\\d+)d$")[, 2])
}

extract_rep <- function(x) {
  paste0("s", str_match(x, "^S(\\d+)_P_")[, 2])
}
extract_rep_num <- function(x) {
  as.integer(str_match(x, "^S(\\d+)_P_")[, 2])
}

dl_days <- c(33, 43, 53, 63)
dl_tree_names <- tree_names[extract_day(tree_names) %in% dl_days]
dl_tree_names <- dl_tree_names[order(extract_day(dl_tree_names), extract_rep_num(dl_tree_names))]

cutoff_day <- 2
calibrate_tip_to_sample_day <- TRUE
n_rarefy <- 500
n_boot <- 2000
thresholds <- c(2, 3)   

set.seed(2026)

compute_metrics_from_counts <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[counts > 0]
  total <- sum(counts)
  
  if (length(counts) == 0 || total == 0) {
    return(c(total = 0, richness = 0, shannon = 0, inv_max = 0))
  }
  
  p <- counts / total
  shannon <- exp(-sum(p * log(p)))
  inv_max <- 1 / max(p)
  
  c(total = total, richness = length(counts), shannon = shannon, inv_max = inv_max)
}

get_node_age <- function(tr, sample_day, calibrate = TRUE) {
  depth <- node.depth.edgelength(tr)
  if (!calibrate) return(depth)
  
  n_tip <- Ntip(tr)
  tip_mean <- mean(depth[1:n_tip], na.rm = TRUE)
  
  if (is.finite(tip_mean) && tip_mean > 0) {
    depth * (sample_day / tip_mean)
  } else {
    depth
  }
}

assign_lineage_cutoff <- function(tr, node_age, cutoff = 2) {
  n_tip <- Ntip(tr)
  n_total <- n_tip + tr$Nnode
  
  parent <- rep(NA_integer_, n_total)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  
  get_internal_ancestors_root2tip <- function(tip_id) {
    anc <- integer(0)
    cur <- tip_id
    while (!is.na(parent[cur])) {
      cur <- parent[cur]
      anc <- c(cur, anc)   
    }
    anc[anc > n_tip]
  }
  
  lineage_id <- vapply(seq_len(n_tip), function(tip_id) {
    anc <- get_internal_ancestors_root2tip(tip_id)
    if (length(anc) == 0) return(NA_character_)
    
    anc_age <- node_age[anc]
    eligible <- anc[anc_age <= cutoff]
    
    if (length(eligible) > 0) {
      chosen <- eligible[which.max(node_age[eligible])]
    } else {
      chosen <- anc[which.min(abs(anc_age - cutoff))]
    }
    as.character(chosen)
  }, character(1))
  
  names(lineage_id) <- tr$tip.label
  lineage_id
}

extract_clone_counts_from_tree <- function(tree_name, cutoff = 2, calibrate = TRUE) {
  tr <- get(tree_name)
  
  if (!inherits(tr, "phylo")) {
    stop(paste0(tree_name, " NO phylo "))
  }
  
  day <- extract_day(tree_name)
  rep_id <- extract_rep(tree_name)
  
  node_age <- get_node_age(tr, sample_day = day, calibrate = calibrate)
  lineage <- assign_lineage_cutoff(tr, node_age = node_age, cutoff = cutoff)
  
  tab <- table(lineage)
  
  tibble(
    sample = tree_name,
    day = day,
    rep = rep_id,
    lineage = names(tab),
    n = as.integer(tab)
  )
}

calc_metrics_table <- function(clone_tbl, min_clone_size = 1) {
  all_samples <- clone_tbl %>% distinct(sample, day, rep)
  
  m <- clone_tbl %>%
    filter(n >= min_clone_size) %>%
    group_by(sample, day, rep) %>%
    summarise(
      total = sum(n),
      richness = n(),
      shannon = {
        p <- n / sum(n)
        exp(-sum(p * log(p)))
      },
      inv_max = {
        p <- n / sum(n)
        1 / max(p)
      },
      .groups = "drop"
    )
  
  all_samples %>%
    left_join(m, by = c("sample", "day", "rep")) %>%
    mutate(
      total = ifelse(is.na(total), 0, total),
      richness = ifelse(is.na(richness), 0, richness),
      shannon = ifelse(is.na(shannon), 0, shannon),
      inv_max = ifelse(is.na(inv_max), 0, inv_max)
    )
}

to_long_metric <- function(df) {
  df %>%
    pivot_longer(
      cols = c(richness, shannon, inv_max),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        richness = "Lineage richness",
        shannon  = "Shannon diversity",
        inv_max  = "1 / max lineage frequency"
      )
    )
}

clone_tbl <- bind_rows(lapply(
  dl_tree_names,
  extract_clone_counts_from_tree,
  cutoff = cutoff_day,
  calibrate = calibrate_tip_to_sample_day
))

raw_metrics <- calc_metrics_table(clone_tbl, min_clone_size = 1)
raw_long <- to_long_metric(raw_metrics)

# depth-matched downsampling + rarefaction

rarefy_sample_means <- function(counts, target_depth, n_iter = 500) {
  counts <- as.integer(counts)
  total <- sum(counts)
  
  if (total < target_depth) {
    stop(" NO target_depth ")
  }
  
  expanded <- rep(seq_along(counts), counts)
  
  res <- replicate(n_iter, {
    draw <- sample(expanded, size = target_depth, replace = FALSE)
    tab <- tabulate(draw, nbins = length(counts))
    met <- compute_metrics_from_counts(tab)
    unname(c(met["richness"], met["shannon"], met["inv_max"]))
  }, simplify = "matrix")
  
  if (is.null(dim(res))) {
    res <- matrix(res, nrow = 3)
  }
  if (nrow(res) != 3) {
    if (ncol(res) == 3) {
      res <- t(res)
    } else {
      stop(" NO rarefaction ")
    }
  }
  
  c(
    richness = mean(res[1, ]),
    shannon  = mean(res[2, ]),
    inv_max  = mean(res[3, ])
  )
}

sample_total <- clone_tbl %>%
  group_by(sample, day, rep) %>%
  summarise(total = sum(n), .groups = "drop")

target_depth <- min(sample_total$total)
cat("\nDepth-matched target depth =", target_depth, "\n")

rarefied_metrics <- clone_tbl %>%
  group_by(sample, day, rep) %>%
  summarise(counts = list(as.integer(n)), .groups = "drop") %>%
  mutate(
    met = purrr::map(
      counts,
      ~ rarefy_sample_means(.x, target_depth = target_depth, n_iter = n_rarefy)
    ),
    richness = purrr::map_dbl(met, ~ .x[["richness"]]),
    shannon  = purrr::map_dbl(met, ~ .x[["shannon"]]),
    inv_max  = purrr::map_dbl(met, ~ .x[["inv_max"]])
  ) %>%
  select(sample, day, rep, richness, shannon, inv_max)

depth_panel_sum <- rarefied_metrics %>%
  to_long_metric() %>%
  group_by(day, metric) %>%
  summarise(
    mean = mean(value),
    lwr  = quantile(value, 0.025),
    upr  = quantile(value, 0.975),
    .groups = "drop"
  )

# leave-one-replicate-out + bootstrap

loo_lines <- bind_rows(lapply(sort(unique(raw_metrics$rep)), function(rp) {
  raw_long %>%
    filter(rep != rp) %>%
    group_by(day, metric) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    mutate(leave_out = rp)
}))

bootstrap_metric <- function(df_metric, B = 2000) {
  ages <- sort(unique(df_metric$day))
  boot_age <- vector("list", B)
  boot_slope <- numeric(B)
  
  for (b in seq_len(B)) {
    db <- bind_rows(lapply(ages, function(a) {
      sub <- df_metric %>% filter(day == a)
      sub[sample(seq_len(nrow(sub)), size = nrow(sub), replace = TRUE), , drop = FALSE]
    }))
    
    boot_age[[b]] <- db %>%
      group_by(day) %>%
      summarise(mean = mean(value), .groups = "drop") %>%
      mutate(iter = b)
    
    fit <- lm(value ~ day, data = db)
    boot_slope[b] <- coef(fit)[["day"]]
  }
  
  boot_age_df <- bind_rows(boot_age)
  
  age_ci <- boot_age_df %>%
    group_by(day) %>%
    summarise(
      lwr = quantile(mean, 0.025),
      upr = quantile(mean, 0.975),
      .groups = "drop"
    )
  
  list(age_ci = age_ci, slope = boot_slope)
}

boot_results <- lapply(unique(raw_long$metric), function(met_name) {
  d <- raw_long %>% filter(metric == met_name)
  bt <- bootstrap_metric(d, B = n_boot)
  
  tibble(
    metric = met_name,
    boot_slope_mean = mean(bt$slope),
    boot_slope_lwr  = quantile(bt$slope, 0.025),
    boot_slope_upr  = quantile(bt$slope, 0.975)
  ) %>%
    mutate(age_ci = list(bt$age_ci))
})

boot_slope_tbl <- bind_rows(lapply(boot_results, function(x) x %>% select(-age_ci)))
boot_ci_df <- bind_rows(lapply(boot_results, function(x) {
  x$age_ci[[1]] %>% mutate(metric = x$metric[1])
}))

loo_slope_tbl <- bind_rows(lapply(unique(raw_long$metric), function(met_name) {
  d <- raw_long %>% filter(metric == met_name)
  
  full_slope <- coef(lm(value ~ day, data = d))[["day"]]
  loo_slopes <- sapply(sort(unique(d$rep)), function(rp) {
    coef(lm(value ~ day, data = d %>% filter(rep != rp)))[["day"]]
  })
  
  tibble(
    metric = met_name,
    full_slope = full_slope,
    loo_slope_min = min(loo_slopes),
    loo_slope_max = max(loo_slopes)
  )
}))

effect_tbl <- loo_slope_tbl %>%
  left_join(boot_slope_tbl, by = "metric")

# min clone-size sensitivity (>=2, >=3 only)

threshold_df <- bind_rows(lapply(thresholds, function(k) {
  calc_metrics_table(clone_tbl, min_clone_size = k) %>%
    select(sample, day, rep, richness, shannon, inv_max) %>%
    mutate(threshold = paste0("Min clone size >= ", k))
}))

threshold_sum <- threshold_df %>%
  to_long_metric() %>%
  group_by(day, metric, threshold) %>%
  summarise(
    mean = mean(value),
    lwr  = quantile(value, 0.025),
    upr  = quantile(value, 0.975),
    .groups = "drop"
  )

# Align factor levels
depth_panel_sum <- depth_panel_sum %>% mutate(metric = factor(metric, levels = metric_levels))
boot_ci_df      <- boot_ci_df %>% mutate(metric = factor(metric, levels = metric_levels))
loo_lines       <- loo_lines %>% mutate(metric = factor(metric, levels = metric_levels))
threshold_sum   <- threshold_sum %>% mutate(metric = factor(metric, levels = metric_levels))
