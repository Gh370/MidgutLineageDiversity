
# Validation of temporal mapping schemes for Dl-GAL4 dataset

rm(list=ls())

library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(readr)
library(phangorn)

load('~/.../data/Phylotime/Phylotime_PosteriorMidgut.RData')

# Collect trees from the workspace

tree_names <- ls(pattern = "^S\\d+_P_\\d+d$")

if (length(tree_names) == 0) {
  stop("No objects matching ^S\\d+_P_\\d+d$ found in the loaded RData.")
}

tub_manifest <- tibble(tree_name = tree_names) %>%
  mutate(
    replicate = as.integer(str_match(tree_name, "^S(\\d+)_P_")[,2]),
    timepoint_day = as.integer(str_match(tree_name, "_(\\d+)d$")[,2]),
    sample_id = tree_name,
    tree = map(tree_name, ~ get(.x, envir = .GlobalEnv))
  )

tub_manifest <- tub_manifest %>%
  filter(timepoint_day %in% c(33, 43, 53, 63))

if (nrow(tub_manifest) == 0) {
  stop("After filtering timepoints, tub_manifest is empty. Check timepoint_day values.")
}

bad <- tub_manifest %>% filter(!map_lgl(tree, ~ inherits(.x, "phylo")))
if (nrow(bad) > 0) {
  stop("Some objects are not of class 'phylo': ", paste(bad$tree_name, collapse = ", "))
}

rf_distance <- function(tr1, tr2) {
  tr1u <- unroot(tr1)
  tr2u <- unroot(tr2)
  as.numeric(dist.topo(tr1u, tr2u, method = "PH85"))
}

run_upgma_invariance_sim <- function(
    n_rep = 1000,
    n_tips = 60,
    offset_range = c(0.1, 10),
    scale_range = c(0.5, 3),
    seed = 1
) {
  set.seed(seed)
  out <- vector("list", n_rep)
  
  for (i in seq_len(n_rep)) {
    tr_u <- rcoal(n_tips)  
    D <- cophenetic(tr_u)
    
    
    tr0 <- phangorn::upgma(as.dist(D))
    
    
    c0 <- runif(1, offset_range[1], offset_range[2])
    D_off <- D
    idx <- row(D_off) != col(D_off)
    D_off[idx] <- D_off[idx] + c0
    diag(D_off) <- 0
    
    D_off <- (D_off + t(D_off)) / 2
    tr_off <- phangorn::upgma(as.dist(D_off))
    
    
    s0 <- runif(1, scale_range[1], scale_range[2])
    D_scl <- D * s0
    tr_scl <- phangorn::upgma(as.dist(D_scl))
    
    out[[i]] <- tibble(
      rep = i,
      offset = c0,
      scale = s0,
      rf_offset = rf_distance(tr0, tr_off),
      rf_scale  = rf_distance(tr0, tr_scl)
    )
  }
  bind_rows(out)
}

sim_res <- run_upgma_invariance_sim(n_rep = 1000, n_tips = 60, seed = 1)

get_root_node <- function(tr) {
  parents <- tr$edge[, 1]
  childs  <- tr$edge[, 2]
  setdiff(parents, childs)[1]
}

compute_node_depth <- function(tr) {
  node.depth.edgelength(tr)
}

apply_time_mapping_tree <- function(tr, tip_age_day, scheme) {
  depths <- compute_node_depth(tr)
  n_tip  <- length(tr$tip.label)
  H <- max(depths[seq_len(n_tip)])  
  
  if (!is.finite(H) || H <= 0) stop("Invalid tree height H; check branch lengths.")
  
  if (scheme == "eclosion_shift") {
    node_time <- depths + (tip_age_day - H)
    cut_map <- function(day_cut) day_cut
  } else if (scheme == "maturation_day2_shift") {
    node_time <- (depths + (tip_age_day - H)) - 2
    cut_map <- function(day_cut) day_cut - 2
  } else if (scheme == "raw_scale") {
    node_time <- depths
    cut_map <- function(day_cut) H * (day_cut / tip_age_day)
  } else if (scheme == "tip_age_scale") {
    node_time <- depths / H * tip_age_day
    cut_map <- function(day_cut) day_cut
  } else {
    stop("Unknown scheme: ", scheme)
  }
  
  list(node_time = node_time, cut_map = cut_map)
}

assign_lineages_at_cutoff <- function(tr, node_time, cutoff_time) {
  n_tip  <- length(tr$tip.label)
  n_node <- tr$Nnode
  n_all  <- n_tip + n_node
  
  root <- get_root_node(tr)
  children <- split(tr$edge[, 2], tr$edge[, 1])
  
  lineage_root <- rep(NA_integer_, n_all)
  lineage_root[root] <- NA_integer_
  
  stack <- list(root)
  while (length(stack) > 0) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    
    ch <- children[[as.character(node)]]
    if (is.null(ch)) next
    
    for (kid in ch) {
      start_new <- (node_time[kid] >= cutoff_time) && (node_time[node] < cutoff_time)
      lineage_root[kid] <- if (start_new) kid else lineage_root[node]
      stack[[length(stack) + 1]] <- kid
    }
  }
  
  tip_lin <- lineage_root[seq_len(n_tip)]
  if (all(is.na(tip_lin))) {
    tip_lin <- rep(root, n_tip)
  } else {
    tip_lin[is.na(tip_lin)] <- seq_len(n_tip)[is.na(tip_lin)]
  }
  
  tibble(tip = seq_len(n_tip), lineage_root = tip_lin)
}

diversity_metrics <- function(lineage_assign_df) {
  tab <- lineage_assign_df %>%
    count(lineage_root, name = "n_cells") %>%
    mutate(p = n_cells / sum(n_cells)) %>%
    arrange(desc(p))
  
  
  p <- tab$p
  p <- p[p > 0]
  
  shannon_entropy <- -sum(p * log(p)) 
  q1 <- exp(shannon_entropy)           
  
  tibble(
    richness = nrow(tab),
    q1 = q1,
    shannon_entropy = shannon_entropy, 
    inv_maxfreq = 1 / max(tab$p),
    dominant_frac = max(tab$p)
  )
}

run_empirical_sensitivity_df <- function(mf, schemes, cutoff_days = 1:5) {
  grid <- expand_grid(scheme = schemes, cutoff_day = cutoff_days)
  
  mf %>%
    crossing(grid) %>%
    mutate(
      mapping = pmap(list(tree, timepoint_day, scheme), apply_time_mapping_tree),
      node_time = map(mapping, "node_time"),
      cut_time = map2_dbl(mapping, cutoff_day, ~ .x$cut_map(.y)),
      lineage_assign = pmap(list(tree, node_time, cut_time), assign_lineages_at_cutoff),
      met = map(lineage_assign, diversity_metrics)
    ) %>%
    unnest(met) %>%
    select(sample_id, timepoint_day, replicate, scheme, cutoff_day,
           richness, q1, shannon_entropy, inv_maxfreq, dominant_frac)
}

tub_schemes <- c("eclosion_shift", "maturation_day2_shift", "raw_scale", "tip_age_scale")
tub_res <- run_empirical_sensitivity_df(tub_manifest, schemes = tub_schemes, cutoff_days = 2)

df_S2B <- tub_res %>%
  filter(cutoff_day == 2) %>%
  mutate(day_f = factor(timepoint_day, levels = c(33, 43, 53, 63)))

df_S2B_mean <- df_S2B %>%
  group_by(scheme, day_f) %>%
  summarise(
    q1_mean = mean(q1, na.rm = TRUE),
    .groups = "drop"
  )

dodge_w <- 0.75

scheme_labels <- c(
  eclosion_shift = "Eclosion anchoring (shift)",
  maturation_day2_shift = "Day-2 anchoring (maturation)",
  raw_scale = "Raw Phylotime scale",
  tip_age_scale = "Tip-age calibration (scaling)"
)

df_S2B <- df_S2B %>%
  mutate(scheme_label = recode(scheme, !!!scheme_labels))

df_S2B_mean <- df_S2B_mean %>%
  mutate(scheme_label = recode(scheme, !!!scheme_labels))

my_cols <- c(
  "Eclosion anchoring (shift)" = "#FBB7C0",
  "Day-2 anchoring (maturation)" = "#4C956C",
  "Raw Phylotime scale" = "#B27092",
  "Tip-age calibration (scaling)" = "#B3AD70"
)

p_S2B <- ggplot() +
  geom_col(
    data = df_S2B_mean ,
    aes(x = day_f, y = q1_mean, fill = scheme_label),
    position = position_dodge(width = dodge_w),
    width = 0.65,
    alpha = 0.65
  ) +
  geom_point(
    data = df_S2B ,
    aes(x = day_f, y = q1, color = scheme_label),
    position = position_jitterdodge(
      jitter.width = 0.18,
      jitter.height = 0,
      dodge.width = dodge_w
    ),
    size = 1.5
  )  +
  scale_fill_manual(values = my_cols, breaks = names(my_cols), drop = FALSE) +
  scale_color_manual(values = my_cols, breaks = names(my_cols), drop = FALSE) +
  scale_y_continuous(
    breaks = c(0, 5, 10),
    expand=c(0,0),
    limits=c(0,11)
  ) +
  labs(
    x = "Days post-eclosion",
    y = "Shannon diversity",
    fill = "Temporal mapping",
    color = "Temporal mapping"
  ) +
  theme_classic()+
  theme(
    axis.title = element_text(face = "bold", size = 11),
    axis.text  = element_text(face = "bold", size = 9),
    
    legend.title = element_text(face = "bold", size = 10),
    legend.text  = element_text(face = "bold", size = 9),
    legend.position = c(0.98, 0.98),  # (x,y) 0-1
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.box = "vertical"
    
  )

p_S2B
