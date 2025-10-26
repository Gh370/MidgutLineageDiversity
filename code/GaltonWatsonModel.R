
# Galton–Watson model reveals mechanisms underlying age‑dependent clonal extinction

rm(list=ls())

library(data.table)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(randomForest)
library(tidyverse)
library(ggtext)
library(patchwork)
library(scales)

set.seed(52637)

# Ne
load('~/.../code/RandomForest_WholeMidgut_data.RData')

exp_ne   <- data[,c(1,5,2)]
colnames(exp_ne) <- c("gen", "age", "Ne")
ages <- sort(unique(exp_ne$age))

sample_1 <- exp_ne[exp_ne$age=='3',]
sample_2 <- exp_ne[exp_ne$age=='13',]
sample_3 <- exp_ne[exp_ne$age=='23',] 
sample_4 <- exp_ne[exp_ne$age=='33',] 

target_division1 <- sort(unique(round(sample_1$gen,0)))
target_division2 <- sort(unique(round(sample_2$gen,0)))
target_division3 <- sort(unique(round(sample_3$gen,0)))
target_division4 <- sort(unique(round(sample_4$gen,0)))

rf_model_1 <- randomForest(Ne ~ gen, data = sample_1, ntree = 500)
rf_predictions_1 <- predict(rf_model_1, newdata = data.frame(gen = target_division1))

result_df1 <- data.frame(
  gen = target_division1,
  age = rep(3,length(target_division1)),
  Ne = rf_predictions_1
)

rf_model_2 <- randomForest(Ne ~ gen, data = sample_2, ntree = 500)
rf_predictions_2 <- predict(rf_model_2, newdata = data.frame(gen = target_division2))

result_df2 <- data.frame(
  gen = target_division2,
  age = rep(13,length(target_division2)),
  Ne = rf_predictions_2
)

rf_model_3 <- randomForest(Ne ~ gen, data = sample_3, ntree = 500)
rf_predictions_3 <- predict(rf_model_3, newdata = data.frame(gen = target_division3))

result_df3 <- data.frame(
  gen = target_division3,
  age = rep(23,length(target_division3)),
  Ne = rf_predictions_3
)

rf_model_4 <- randomForest(Ne ~ gen, data = sample_4, ntree = 500)
rf_predictions_4 <- predict(rf_model_4, newdata = data.frame(gen = target_division4))

result_df4 <- data.frame(
  gen = target_division4,
  age = rep(33,length(target_division4)),
  Ne = rf_predictions_4
)

exp_ne <- do.call(rbind, list(result_df1, result_df2, result_df3,result_df4))

gen_max_list <- c("3" = 23, "13" = 24, "23" = 30, "33" = 35)

# Galton–Watson model
gw_sim2 <- function(gen_max, n0, p2_0, lam1, lam2, g_star, 
                    p0_y = 0.05, p0_inc = 0.25, 
                    simulate = TRUE, return_lin = TRUE, age = NULL) {
  
  age <- as.numeric(age)
  print(paste("Age in gw_sim2:", age))
  
  if(age == 23) {
    stable_period <- 15  
    lam1 <- lam1 * 0.5   
    lam2 <- lam2 * 0.7   
  }
  
  if(!simulate) {
    Ne <- numeric(gen_max + 1); Ne[1] <- n0
    for(g in seq_len(gen_max)) {
      
      if (age == 23 && g > stable_period) {
        decay <- lam1 * stable_period  
      } else {
        decay <- if(g <= g_star) lam1 * g else lam1 * g_star + lam2 * (g - g_star)
      }
      p2 <- p2_0 * exp(-decay)
      p0 <- p0_y + p0_inc * (1 - exp(-lam2 * max(g - g_star, 0)))
      if (p2 + p0 >= 1) { 
        fac <- 0.999 / (p2 + p0)
        p2 <- p2 * fac
        p0 <- p0 * fac
      }
      Ne[g + 1] <- Ne[g] * (1 - p0 + p2)
    }
    if (return_lin) stop("no lin_size，simulate=TRUE")
    return(Ne)
  }
  
  lin <- rep(1, n0)
  Ne <- numeric(gen_max + 1); Ne[1] <- n0
  lineage_sizes <- list()  
  for(g in seq_len(gen_max)) {
    
    if (age == 23 && g > stable_period) {
      decay <- lam1 * stable_period  
    } else {
      decay <- if(g <= g_star) lam1 * g else lam1 * g_star + lam2 * (g - g_star)
    }
    p2 <- p2_0 * exp(-decay)
    p0 <- p0_y + p0_inc * (1 - exp(-lam2 * max(g - g_star, 0)))
    if (p2 + p0 >= 1) { 
      fac <- 0.999 / (p2 + p0)
      p2 <- p2 * fac
      p0 <- p0 * fac
    }
    p1 <- 1 - p0 - p2
    prob <- c(p0, p1, p2)
    lin <- sapply(lin, \(c) if(c == 0) 0 else sum(sample(c(0, 1, 2), c, TRUE, prob)))
    
    
    lineage_sizes[[g]] <- lin
    
    Ne[g + 1] <- sum(lin)
  }
  if (return_lin) return(list(Ne = Ne, lin_size = lin, lineage_sizes = lineage_sizes))  # 返回每个谱系的大小
  return(Ne)
}

# Loss
param_age <- exp_ne %>% 
  group_by(age) %>% 
  group_split() %>% 
  map_dfr(function(df_age) {
    age <- unique(df_age$age)
    gen_max <- gen_max_list[as.character(age)]
    n0_fix <- df_age$Ne[df_age$gen == 0][1]
    
    peak_obs <- df_age$gen[which.max(df_age$Ne)]  
    g_slope <- peak_obs + 5                       
    slope_obs <- (log(df_age$Ne[df_age$gen == g_slope] + 1) - log(df_age$Ne[df_age$gen == peak_obs] + 1)) / 5
    
    w_vec <- {
      center <- peak_obs
      sigma <- 0.25 * gen_max
      w <- dnorm(df_age$gen, mean = center, sd = sigma)
      w[df_age$gen >= 0.8 * gen_max] <- w[df_age$gen >= 0.8 * gen_max] * 1.5
      w / max(w)
    }
    
    loss <- function(par) {
      p2_0 <- par[1]
      lam1 <- par[2]
      lam2 <- par[3]
      gs <- par[4]
      
      
      if (p2_0 < .35 || p2_0 > .95 || lam1 < .005 || lam1 > .30 || lam2 < 1e-4 || lam2 > .30 || gs < 3 || gs > 25) return(1e12)
      
      
      sim <- gw_sim2(gen_max, n0_fix, p2_0, lam1, lam2, gs, simulate = FALSE, return_lin = FALSE, age = age)
      if (any(!is.finite(sim))) return(1e12)
      
     
      rss <- sum(w_vec * (log10(sim[df_age$gen + 1] + 1) - log10(df_age$Ne + 1))^2)
      
      
      peak_sim <- which.max(sim) - 1
      peak_pen <- 30 * (peak_sim - peak_obs)^2
      
      
      if (peak_obs + 5 > gen_max || peak_sim + 5 > gen_max) {
        slope_pen <- 0  
      } else {
        slope_sim <- (log(sim[peak_sim + 6] + 1) - log(sim[peak_sim + 1] + 1)) / 5
        slope_pen <- 50 * (slope_sim - slope_obs)^2
      }
      
      
      idx_tail <- which(df_age$gen >= peak_obs + 6)
      if (length(idx_tail) == 0 || any(idx_tail + 1 > length(sim))) {
        tail_pen2 <- 0
      } else {
        tail_pen2 <- 100 * mean((log10(sim[idx_tail + 1] + 1) - log10(df_age$Ne[idx_tail] + 1))^2)
      }
      
      
      total <- rss + peak_pen + slope_pen + tail_pen2
      if (!is.finite(total)) total <- 1e12
      total
    }
    
    
    set.seed(age)
    starts <- matrix(runif(80 * 4), 80, 4)
    starts[,1] <- 0.35 + starts[,1] * 0.60   
    starts[,2] <- 0.005 + starts[,2] * 0.295  
    starts[,3] <- 0.0001 + starts[,3] * 0.2999 
    starts[,4] <- 3 + starts[,4] * 22  
    
    best <- starts[1,]; best_rss <- Inf
    for(i in seq_len(nrow(starts))) {
      opt <- optim(starts[i,], loss, method = "Nelder-Mead", control = list(maxit = 1500))
      if(opt$value < best_rss) {
        best_rss <- opt$value
        best <- opt$par
      }
    }
    
    tibble(
      age = age,
      n0 = n0_fix,
      p2_0 = best[1],
      lam1 = best[2],
      lam2 = best[3],
      g_star = best[4],
      rss = best_rss
    )
  })

fwrite(param_age, "parameters_twoSlope_fixN0.csv")
print(param_age)

# Fit the model and plot the fitted curve

sim_curves <- param_age %>% 
  pmap_dfr(function(age, n0, p2_0, lam1, lam2, g_star, rss) {
    gen_max <- gen_max_list[as.character(age)]
    
    Ne_sim <- gw_sim2(gen_max, n0, p2_0, lam1, lam2, g_star, simulate = FALSE, return_lin = FALSE, age = age)
    tibble(
      age = age,
      gen = 0:gen_max,
      Ne_sim = Ne_sim  
    )
  })

# Figure 2d-2e
load('~/.../code/RandomForest_WholeMidgut_data.RData')
exp_ne1   <- data[,c(1,5,2)]
colnames(exp_ne1) <- c("gen", "age", "Ne")
ages <- sort(unique(exp_ne1$age))

exp_ne1_33 <- exp_ne1[exp_ne1$age==33,]
sim_curves_33 <- sim_curves[sim_curves$age==33,]

exp_ne1_3 <- exp_ne1[exp_ne1$age!=33,]
sim_curves_3 <- sim_curves[sim_curves$age!=33,]

# Figure 2d

plot_2d <- ggplot() +
  ## Observed
  geom_line(
    data = exp_ne1_33,
    aes(gen, Ne, linetype = "Observed"),
    size   = 0.6,
    colour = "#457b9d",
    alpha  = 0.8
  ) +
  ## Predicted
  geom_line(
    data = sim_curves_33,
    aes(gen, Ne_sim, linetype = "Predicted (Galton–Watson model)"),
    size   = 0.6,
    colour = "#99582a"
  ) + 
  
  ## Facet
  facet_wrap(
    ~ age,
    nrow   = 2,
    scales = "free_y",
    labeller = labeller(age = \(x) paste0("Day ", x))
  ) +
  
  
  scale_y_log10(
    limits = c(1, 2000),     
    breaks = 10^(0:3),        
    labels = parse(text = paste0("10^", 0:3))
  )  + scale_x_continuous(breaks = seq(0, 35, by = 5))+
  
  
  scale_linetype_manual(
    values = c("Observed" = "solid",
               "Predicted (Galton–Watson model)" = "dashed"),
    name   = NULL
  ) +
  
  
  theme_minimal(base_size = 12) +
  theme(
    panel.border   = element_rect(fill = NA, colour = "black"),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    
    
    axis.text      = element_text(size = 9,  colour = "black"),
    axis.title     = element_text(size = 11, colour = "black"),
    
    strip.background = element_rect(fill = NA, colour = "black"),
    strip.text       = element_text(size = 10, colour = "black"),
    legend.text  = element_text(size = 9,  colour = "black"),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.margin       = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    plot.margin         = ggplot2::margin(-4, 10, 10, 10, unit = "pt"), 
    plot.title = element_blank()
  ) +
  
  labs(
    x = "Lineage distance",
    y = "Effective population size"
  )

plot_2d

# Simulation
simulations <- 1000
gen_max_33 <- gen_max_list["33"]
param_33 <- param_age[param_age$age == 33,]

n0_33 <- 100
p2_0_33 <- param_33$p2_0
lam1_33 <- param_33$lam1
lam2_33 <- param_33$lam2
g_star_33 <- param_33$g_star

lineage_diversities <- list()

for (i in 1:simulations) {
  sim_result <- gw_sim2(gen_max_33, n0_33, p2_0_33, lam1_33, lam2_33, g_star_33, simulate = TRUE, return_lin = TRUE, age = 33)
  
  
  lineage_sizes <- sim_result$lineage_sizes
  
  
  diversity_per_generation <- sapply(lineage_sizes, function(lin_size) {
    
    p <- lin_size / sum(lin_size)  
    p <- p[p!=0]
    diversity <- exp(-sum(p * log(p)))  
    return(diversity)
  })
  
  
  lineage_diversities[[i]] <- diversity_per_generation
}


lineage_diversities_df <- as.data.frame(t(simplify2array(lineage_diversities)))

colnames(lineage_diversities_df) <- 1:gen_max_33 


lineage_diversities_long <- data.frame(simulation = integer(), gen = integer(), diversity = numeric())


for (sim in 1:simulations) {
  
  for (gen in 1:gen_max_33) {
    lineage_diversities_long <- rbind(lineage_diversities_long, 
                                      data.frame(simulation = sim, 
                                                 gen = gen, 
                                                 diversity = lineage_diversities_df[sim, gen]))
  }
}


ggplot(lineage_diversities_long, aes(x = gen, y = diversity,group=gen)) +
  geom_boxplot() +
  labs(x = "Generation", y = "Lineage Diversity (Shannon Index)", title = "Lineage Diversity Over Generations (1000 Simulations)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

# Figure 2e
max(lineage_diversities_long$diversity)
min(lineage_diversities_long$diversity)

p <- ggplot(lineage_diversities_long, aes(x = gen, y = diversity,group=gen)) +
  geom_boxplot(size=0.3,outlier.size = 0.1,fill = "white", color = "black", alpha = 0.2) +
  labs(x = "Lineage distance", y = "Shannon diversity\nin simulated cell populations") + scale_y_continuous(limits =c(50,100),breaks=c(50,60,70,80,90,100))+ scale_x_continuous(breaks = seq(0, 35, by = 5))

plot_2e <- p+theme_bw()+theme( strip.text = element_text(size = 9),       
                               panel.grid.major = element_blank(),         
                               panel.grid.minor = element_blank(),
                               panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
                               
                               axis.line.x = element_blank(),  
                               axis.line.y = element_blank(),
                               axis.text.y = element_text(size = 9,color='black'),axis.text.x = element_text(size = 9,color='black'),axis.title.y = element_text(size = 11,vjust=2,color='black'),axis.title.x = element_text(size = 11,color='black'),legend.position="none",strip.text.y = element_text(angle = 0))
plot_2e












