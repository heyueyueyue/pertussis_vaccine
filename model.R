library(rstan)
library(ggplot2)

# options ---------------------------------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

# load data ---------------------------------------------------------------

## read age mortality data
data_mortality <- read.csv("../data/mortality_group.csv")
age_mortality <- data_mortality$Estimated_Mortality/1000
mu_age <- 1 / (data_mortality$age_durations * 365)
N <- nrow(data_mortality)
age_breaks <- sort(unique(c(data_mortality$Age_start, data_mortality$Age_end)))

## read contact matrix
contact_matrix <- read.csv("../data/contact_changsha_baseline.csv")
rownames(contact_matrix) <- contact_matrix$X
contact_matrix <- as.matrix(contact_matrix[, -1])
colnames(contact_matrix) <- rownames(contact_matrix)

## read incidence data
data_incidence <- read.csv("../data/incidence_2020_2024.csv")
data_incidence <- data_incidence |> 
     mutate(Date = as.Date(Date, format = "%Y/%m/%d"),
            AgeStart = case_when(str_detect(Age, '月') ~ as.numeric(str_extract(Age, '\\d+'))/12 - 1/12,
                                 Age == '1岁' ~ 11/12,
                                 str_detect(Age, '岁') ~ as.numeric(str_extract(Age, '\\d+')) - 1,
                                 TRUE ~ as.numeric(str_extract(Age, '\\d+'))),
            AgeEnd = case_when(str_detect(Age, '月') ~ as.numeric(str_extract(Age, '\\d+'))/12,
                               Age == '1岁' ~ 1,
                               str_detect(Age, '岁') ~ as.numeric(str_extract(Age, '\\d+')),
                               TRUE ~ as.numeric(str_extract(Age, '\\d+')) + 1),
            # sample one age between AgeStart and AgeEnd
            AgeModel = AgeStart + runif(n(), 0, 1) * (AgeEnd - AgeStart),
            AgeModel = round(AgeModel, 4),
            # group age
            AgeGroup = cut(AgeModel, breaks = age_breaks, labels = 1:N),
            Year = year(Date),
            Month = month(Date))

## overview of the data
data_incidence |> 
     group_by(Year, Month) |>
     summarise(Incidence = n(),
               .groups = 'drop') |>
     arrange(Year, Month) |> 
     mutate(Date = as.Date(paste(Year, Month, 1, sep = "-"))) |>
     ggplot(aes(x = Date, y = Incidence)) +
     geom_line() +
     theme_minimal() +
     labs(title = "Incidence of pertussis",
          x = "Date",
          y = "Incidence")

## subset the data: 2023 - 2024
data_incidence_model <- data_incidence |> 
     filter(Year >= 2023) |> 
     group_by(AgeGroup, Date) |>
     summarise(Incidence = n(),
               .groups = 'drop') |> 
     mutate(Index = as.integer(Date - min(Date) + 1)) |> 
     complete(AgeGroup, Index = 1:730, fill = list(Incidence = 0))

incidence_matrix <- matrix(data_incidence_model$Incidence, nrow = N)

## Stan model --------------------------------------------------

# estimate the beta
R0 <- 5.5
sigma <- 0.2/2
gamma <- 0.1/2
alpha <- 1/log(8 * 365)
delta <- mean(age_mortality)
contact <- mean(contact_matrix)

beta <- R0 * (sigma + delta) * (gamma + delta) / (sigma * contact)

data_list <- list(
     N = N,
     contact_matrix = contact_matrix,
     delta = age_mortality,
     v = 0.95,  # 疫苗接种率
     B = 41.20*10^4/365,  # 新生人口数量（假设每年出生100人）
     N_pop = 6640*10^4,  # 人口总数
     beta = beta,
     lambda = 0.002,  # 假设疫苗效果每天衰减0.2%
     mu = mu_age,
     'T' = 730,  # 总天数
     observed_cases = incidence_matrix
)

# compile the model
sirs_model <- stan(file = "model.stan", data = data_list, iter = 2000, chains = 4)

# extract the results
sirs_results <- rstan::extract(sirs_model)

# extract the posterior samples
posterior_beta <- sirs_results$beta
posterior_sigma <- sirs_results$sigma
posterior_gamma <- sirs_results$gamma
posterior_alpha <- sirs_results$alpha

# extract the model output
S_output <- sirs_results$y[, 1:16] 
I_output <- sirs_results$y[, 49:64]
E_output <- sirs_results$y[, 17:32]
R_output <- sirs_results$y[, 65:80]
time_steps <- 1:nrow(S_output)
age_groups <- rep(1:16, times = length(time_steps))

# compute the mean and 95% CI for each state
compute_stats <- function(data) {
     mean_vals <- apply(data, 1, mean)
     lower_CI <- apply(data, 1, function(x) quantile(x, 0.025))
     upper_CI <- apply(data, 1, function(x) quantile(x, 0.975))
     
     data.frame(
          Mean = mean_vals,
          LowerCI = lower_CI,
          UpperCI = upper_CI
     )
}

S_stats <- compute_stats(S_output)
I_stats <- compute_stats(I_output)
E_stats <- compute_stats(E_output)
R_stats <- compute_stats(R_output)

plot_data <- data.frame(
     Time = rep(time_steps, each = 16),
     AgeGroup = rep(1:16, times = length(time_steps)),
     Susceptible_Mean = as.vector(S_stats$Mean),
     Susceptible_LowerCI = as.vector(S_stats$LowerCI),
     Susceptible_UpperCI = as.vector(S_stats$UpperCI),
     Infected_Mean = as.vector(I_stats$Mean),
     Infected_LowerCI = as.vector(I_stats$LowerCI),
     Infected_UpperCI = as.vector(I_stats$UpperCI),
     Exposed_Mean = as.vector(E_stats$Mean),
     Exposed_LowerCI = as.vector(E_stats$LowerCI),
     Exposed_UpperCI = as.vector(E_stats$UpperCI),
     Recovered_Mean = as.vector(R_stats$Mean),
     Recovered_LowerCI = as.vector(R_stats$LowerCI),
     Recovered_UpperCI = as.vector(R_stats$UpperCI)
)

plot_data_1 <- plot_data[plot_data$Time <= 365*2, ]

# visualization -----------------------------------------------------------

fig <- ggplot(plot_data_1,
       mapping = aes(x = Time)) +
     geom_line(aes(y = Susceptible_Mean, color = "Susceptible")) +
     geom_ribbon(aes(ymin = Susceptible_LowerCI, ymax = Susceptible_UpperCI, fill = "Susceptible"), alpha = 0.2) +
     geom_line(aes(y = Infected_Mean, color = "Infected")) +
     geom_ribbon(aes(ymin = Infected_LowerCI, ymax = Infected_UpperCI, fill = "Infected"), alpha = 0.2) +
     geom_line(aes(y = Exposed_Mean, color = "Exposed")) +
     geom_ribbon(aes(ymin = Exposed_LowerCI, ymax = Exposed_UpperCI, fill = "Exposed"), alpha = 0.2) +
     geom_line(aes(y = Recovered_Mean, color = "Recovered")) +
     geom_ribbon(aes(ymin = Recovered_LowerCI, ymax = Recovered_UpperCI, fill = "Recovered"), alpha = 0.2) +
     facet_wrap(~AgeGroup, scales = "free_y") +  # 按年龄组分面
     labs(title = "Epidemic Model - SEIRS-VSEIR Dynamics Over Time",
          x = "Time (in steps)",
          y = "Population Count",
          color = "State",
          fill = "State") +
     theme_minimal() +
     theme(legend.position = "bottom")

ggsave(fig,
       filename = "../outcome/epidemic_model.png",
       width = 12,
       height = 8, 
       dpi = 300)
