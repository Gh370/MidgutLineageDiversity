# Validation of random forest-based framework for estimating effective population size trends

rm(list=ls())

library(randomForest)
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)

set.seed(123)

target_division <- seq(0.5,35, length.out = 100)

# Initialize result storage
simulated_means <- matrix(NA, nrow = length(target_division), ncol = 1000)
simulated_ci_lower <- matrix(NA, nrow = length(target_division), ncol = 1000)
simulated_ci_upper <- matrix(NA, nrow = length(target_division), ncol = 1000)
originaldata_ci_width <- vector()
originaldata_ci_div <- vector()

# Generate simulated data and repeatedly execute the framework
for (i in 1:1000) {
  
  generate_random_sample <- function(n, division_range) {
    
    division <- sort(runif(n, min(division_range), max(division_range)))
    
    # Generate data in stages
    n_up <- floor(n * 0.3)     
    n_flat <- floor(n * 0.4)    
    n_down <- n - n_up - n_flat 
    
    mean_up <- seq(200, 1000, length.out = n_up) + rnorm(n_up, 0, 50)
    
    mean_flat <- rep(1000, n_flat) + rnorm(n_flat, 0, 50)
    
    mean_down <- seq(1000, 500, length.out = n_down) + rnorm(n_down, 0, 50)
    
    mean <- c(mean_up, mean_flat, mean_down)
    mean <- pmax(mean, 100) 
    
    ci_width <- runif(n, 0.1, 0.2) * mean
    
    ci_lower <- pmax(mean - ci_width, 1e-6)  # 保证下界为正
    ci_upper <- mean + ci_width
    
    return(data.frame(division = division, mean = mean, ci_lower = ci_lower, ci_upper = ci_upper))
  }
  
  # Generate simulated sample data
  simulated_sample <- generate_random_sample(100, c(0.5, 35))
  originaldata_ci_width <- append(originaldata_ci_width,simulated_sample$mean)
  originaldata_ci_div <- append(originaldata_ci_div,simulated_sample$division)
  
  repeat {
    
    cut_points <- sort(c(0, runif(2, min = 0, max = 1), 1))

    random_values <- diff(cut_points)
    
    group_labels <- sample(1:3, size = 100, replace = TRUE, prob = random_values)
    grouped_data <- split(1:100, group_labels)
    
    if (length(grouped_data) == 3) break
  }
    
  simulated_sample_1 <- simulated_sample[grouped_data[[1]],]
  simulated_sample_2 <- simulated_sample[grouped_data[[2]],]
  simulated_sample_3 <- simulated_sample[grouped_data[[3]],]
  
  simulated_data <- rbind(simulated_sample_1, simulated_sample_2, simulated_sample_3)
  
  # Random forest fitting
  rf_model <- randomForest(mean ~ division, data = simulated_data, ntree = 500)
  rf_predictions <- predict(rf_model, newdata = data.frame(division = target_division))
  
  # Generate confidence intervals using the distribution of prediction residuals
  simulate_ci <- function(predictions, residuals, alpha = 0.05, n_sim = 1000) {
   
    sim_values <- replicate(n_sim, predictions + sample(residuals, size = length(predictions), replace = TRUE))
    ci_lower <- apply(sim_values, 1, quantile, probs = alpha / 2)
    ci_upper <- apply(sim_values, 1, quantile, probs = 1 - alpha / 2)
    return(list(lower = ci_lower, upper = ci_upper))
  }
  
  rf_residuals <- simulated_data$mean - predict(rf_model, newdata = simulated_data)
  bayesian_ci <- simulate_ci(rf_predictions, rf_residuals)
  
  smoothed_mean <- loess(rf_predictions ~ target_division, span = 0.3)$fitted
  smoothed_ci_lower <- loess(bayesian_ci$lower ~ target_division, span = 0.3)$fitted
  smoothed_ci_upper <- loess(bayesian_ci$upper ~ target_division, span = 0.3)$fitted
  
  # Gamma distribution correction for the lower bound
  gamma_ci_lower <- qgamma(p = 0.025, 
                           shape = (smoothed_mean^2) / (smoothed_mean * 0.2)^2, 
                           rate = smoothed_mean / (smoothed_mean * 0.2)^2)
  
  smoothed_ci_lower <- pmax(smoothed_ci_lower, gamma_ci_lower)
  
  gamma_ci_upper <- qgamma(p = 0.975, 
                           shape = (smoothed_mean^2) / (smoothed_mean * 0.2)^2, 
                           rate = smoothed_mean / (smoothed_mean * 0.2)^2)
  
  smoothed_ci_upper <- pmin(smoothed_ci_upper, gamma_ci_upper)
  
  smoothed_ci_lower <- pmax(smoothed_ci_lower, 1e-6) 
  smoothed_ci_upper <- pmax(smoothed_ci_upper, smoothed_ci_lower)
  
  simulated_means[, i] <- smoothed_mean
  simulated_ci_lower[, i] <- smoothed_ci_lower
  simulated_ci_upper[, i] <- smoothed_ci_upper
}

# plot
melted_means <- melt(simulated_means)
colnames(melted_means) <- c("Division_Index", "Simulation", "Mean")

melted_lower <- melt(simulated_ci_lower)
colnames(melted_lower) <- c("Division_Index", "Simulation", "lower")
melted_upper <- melt(simulated_ci_upper)
colnames(melted_upper) <- c("Division_Index", "Simulation", "upper")
melted_lower$upper <- melted_upper$upper

result_summary <- data.frame(
  division = originaldata_ci_div,
  mean = originaldata_ci_width)

z12 <- ggplot(melted_means)+ geom_ribbon(data=melted_lower,aes(x = target_division[Division_Index], ymin = lower, ymax = upper,group = Simulation),fill='#f5cac3',colour = NA, alpha = 0.6)+geom_line(aes(x = target_division[Division_Index], y = Mean, group = Simulation),alpha = 0.1, color = "#403d39")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_x_continuous(limits=c(0,35.2),breaks=c(5,10,15,20,25,30),expand=c(0,0))+ylab("Predicted population size")+xlab("Lineage distance")
z12
z12 <- z12+theme_bw()+theme(panel.grid.major = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 9,face="bold"),
                            panel.grid.minor = element_blank(), 
                            legend.background = element_blank(),axis.text.y = element_text(size = 25,face = "bold",color='black'),axis.text.x = element_text(size = 25,face = "bold",color='black'),axis.title.y = element_text(size = 28,vjust=2,face = "bold",color='black'),axis.title.x = element_text(size = 28,face = "bold",color='black'),legend.position="top",legend.direction = "horizontal", plot.title = element_text(size=12, face="bold",hjust=0.5))+scale_y_continuous(expand=c(0,0),limits=c(0,1500),breaks=c(250,500,750,1000,1250))
z12

avg_mean <- rowMeans(simulated_means)
result_summary1 <- data.frame(
  division = target_division,
  mean = avg_mean
)

z1 <- ggplot(result_summary1)+ geom_point(data = result_summary, aes(x = division, y = mean), color='grey',shape=21,fill='white',size=0.7)+geom_line(aes(x = division, y = mean),color = "#403d39",lty=2,size=2)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_x_continuous(limits=c(0,35.2),breaks=c(5,10,15,20,25,30),expand=c(0,0))+ylab("True population size")+xlab("Lineage distance")+scale_y_continuous(limits=c(0,1500),breaks=c(250,500,750,1000,1250),expand=c(0,0))
z1
z1 <- z1+theme_bw()+theme(panel.grid.major = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 9,face="bold"),
                            panel.grid.minor = element_blank(), 
                            legend.background = element_blank(),axis.text.y = element_text(size = 25,face = "bold",color='black'),axis.text.x = element_text(size = 25,face = "bold",color='black'),axis.title.y = element_text(size = 28,vjust=2,face = "bold",color='black'),axis.title.x = element_text(size = 28,face = "bold",color='black'),legend.position="top",legend.direction = "horizontal", plot.title = element_text(size=12, face="bold",hjust=0.5))
z1








