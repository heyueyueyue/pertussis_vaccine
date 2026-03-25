
library(splines)
library(tidyverse)

# read data
data <- openxlsx::read.xlsx("../data/mortality_2010.xlsx")
data$Age_mid <- (data$Age_start + data$Age_end) / 2

age_values <- seq(min(data$Age_start), max(data$Age_end), by = 0.01)
spline_model <- lm(Mortality ~ ns(Age_mid, df = 10), data = data)
predicted_mortality <- predict(spline_model, newdata = data.frame(Age_mid = age_values))
predicted_mortality <- pmax(predicted_mortality, 0)

plot(age_values, predicted_mortality, type = 'l', xlab = "Age", ylab = "Mortality")

# create age groups
Age_start = c(0, 1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 1.5, 6, 12, 21, 31, 41, 51, 61, 71)

age_groups <- tibble(
     Age_start = Age_start,
     Age_end = c(Age_start[-1], 100)
)

# estimate mortality for each age group
estimated_mortality <- sapply(1:nrow(age_groups), function(i) {
     group_start <- age_groups$Age_start[i]
     group_end <- age_groups$Age_end[i]
     mask <- age_values >= group_start & age_values <= group_end
     
     group_mortality <- predicted_mortality[mask]
     age_subset <- age_values[mask]
     
     # Calculate the weight for each age
     weights <- c(diff(age_subset), last = 1)  # Assuming equal weight for the last interval
     weighted_mortality <- sum(group_mortality * weights) / sum(weights)
     
     return(weighted_mortality)
})

result_df <- data.frame(
     Age_start = age_groups$Age_start,
     Age_end = age_groups$Age_end,
     Estimated_Mortality = estimated_mortality
) |> 
     mutate(age_durations = Age_end - Age_start)

# visualize the results
ggplot() +
     geom_step(data = data,
               aes(x = Age_end, y = Mortality, color = 'Original'),
               direction = 'vh') +
     geom_line(aes(x = age_values, y = predicted_mortality, color = 'Predicted')) +
     geom_step(data = result_df,
               aes(x = Age_end, y = Estimated_Mortality, color = 'Estimated'),
               direction = 'vh') +
     scale_x_continuous(breaks = seq(0, 100, by = 10))+
     scale_y_continuous(limits = c(0, NA),
                        expand = expansion(mult = c(0, 0.1))) +
     labs(x = "Age", y = "Mortality", color = NULL) +
     theme(legend.position = "bottom")


ggsave("../outcome/fig_s1_age_mortality.png",
       width = 8, height = 6,
       dpi = 300)

write.csv(result_df,
          "../data/mortality_group.csv",
          row.names = FALSE)
