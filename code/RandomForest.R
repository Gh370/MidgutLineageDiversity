# Framework for modeling Ne dynamics using random forest 

rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(scales)
library(randomForest)

load('~/.../code/Ne_WholeMidgut_data.RData')

# Take the samples obtained on day 3 (S1_3d, S2_3d and S3_3d) as an example
data_1 <- data[data$age=='3',]

colnames(data_1) <- c('division','mean','ci_upper','ci_lower','div','age','sample')

sample_1 <- data_1[data_1$sample=='S1_3d',1:4]
sample_2 <- data_1[data_1$sample=='S2_3d',1:4]
sample_3 <- data_1[data_1$sample=='S3_3d',1:4] 

# Determine the range
target_division <- sort(unique(c(sample_1$division, sample_2$division, sample_3$division)))
length(target_division)
nrow(sample_1)+nrow(sample_2)+nrow(sample_3)

# Random forest fitting
all_data <- do.call(rbind, list(sample_1, sample_2, sample_3))
rf_model <- randomForest(mean ~ division, data = all_data, ntree = 500)
rf_predictions <- predict(rf_model, newdata = data.frame(division = target_division))

# Generate confidence intervals using the distribution of prediction residuals
simulate_ci <- function(predictions, residuals, alpha = 0.05, n_sim = 1000) {
  
  sim_values <- replicate(n_sim, predictions + sample(residuals, size = length(predictions), replace = TRUE))
  ci_lower <- apply(sim_values, 1, quantile, probs = alpha / 2)
  ci_upper <- apply(sim_values, 1, quantile, probs = 1 - alpha / 2)
  return(list(lower = ci_lower, upper = ci_upper))
}

rf_residuals <- all_data$mean - predict(rf_model, newdata = all_data)
bayesian_ci <- simulate_ci(rf_predictions, rf_residuals)

smoothed_mean <- loess(rf_predictions ~ target_division, span = 0.3)$fitted
smoothed_bayesian_ci_lower <- loess(bayesian_ci$lower ~ target_division, span = 0.3)$fitted
smoothed_bayesian_ci_upper <- loess(bayesian_ci$upper ~ target_division, span = 0.3)$fitted

# Gamma distribution correction for the lower bound
gamma_ci_lower <- qgamma(p = 0.025, 
                         shape = (smoothed_mean^2) / (smoothed_mean * 0.2)^2, 
                         rate = smoothed_mean / (smoothed_mean * 0.2)^2)

smoothed_ci_lower <- pmax(smoothed_bayesian_ci_lower, gamma_ci_lower)

gamma_ci_upper <- qgamma(p = 0.975, 
                         shape = (smoothed_mean^2) / (smoothed_mean * 0.2)^2, 
                         rate = smoothed_mean / (smoothed_mean * 0.2)^2)

smoothed_ci_upper <- pmin(smoothed_bayesian_ci_upper, gamma_ci_upper)

smoothed_ci_lower <- pmax(smoothed_ci_lower, 1e-6)  
smoothed_ci_upper <- pmax(smoothed_ci_upper, smoothed_ci_lower) 

# Result integration
result_df1 <- data.frame(
  division = target_division,
  mean = smoothed_mean,
  ci_lower = smoothed_ci_lower,
  ci_upper = smoothed_ci_upper
)

# plot
p <- ggplot(result_df1,aes(x = division)) + geom_ribbon(aes(ymin=ci_lower,ymax=ci_upper),fill='#f5cac3',colour = NA, alpha = 0.6)+geom_line(aes(y=mean),linetype="solid",linewidth=1)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_x_continuous(limits=c(0,25.2),breaks=c(5,10,15,20,25),expand=c(0,0))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(1,5000))+ylab("Effective population size")+xlab("Lineage distance")
p
p <- p+theme_bw()+theme(panel.grid.major = element_line(color = "grey90", linewidth = 1),legend.title = element_blank(),legend.text = element_text(size = 9,face="bold"),
                            panel.grid.minor = element_blank(),
                            legend.background = element_blank(),axis.text.y = element_text(size = 9,face = "bold",color='black'),axis.text.x = element_text(size = 9,face = "bold",color='black'),axis.title.y = element_text(size = 11,vjust=2,face = "bold",color='black'),axis.title.x = element_text(size = 11,face = "bold",color='black'),legend.position="top",legend.direction = "horizontal", plot.title = element_text(size=12, face="bold",hjust=0.5))
p

