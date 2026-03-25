library(tidyverse)
library(openxlsx)
library(socialmixr)
library(data.table)
library(patchwork)

# reading data
participants <- read.xlsx("./data/contact_changsha_baseline.xlsx", sheet = "Participants")
contacts <- read.xlsx("./data/contact_changsha_baseline.xlsx", sheet = "Contacts")
population <- read.xlsx("./data/contact_changsha_baseline.xlsx", sheet = "Population")

# Convert age intervals into a usable format for `cut` function
population$age_limits <- sub("\\[", "", population$age)
population$age_limits <- sub("\\)", "", population$age_limits)
age_breaks <- as.numeric(unlist(strsplit(paste(population$age_limits, collapse=","), ",")))
age_breaks <- unique(age_breaks)
age_breaks <- sort(c(age_breaks, max(age_breaks) + 5))  # Ensure the last interval is included

# Prepare participant data
participants <- participants |> 
     mutate(country = "China",
            age_group = cut(part_age, breaks = age_breaks, right = FALSE, labels = FALSE))

# Prepare contacts data
set.seed(20241227)

contacts <- contacts |> 
     filter(!is.na(cnt_age_est_min) & !is.na(cnt_age_est_max)) |>
     mutate(cnt_age_exact = map2_dbl(cnt_age_est_min, cnt_age_est_max, ~ round(runif(1, min = .x, max = .y))))

# Function to calculate bootstrapped contact matrices
calculate_matrices <- function(data_p, data_c, population, breaks, n_bootstrap = 100) {
     pop_freq <- population$freq / sum(population$freq)
     sample_sizes <- round(nrow(data_p) * pop_freq)
     
     bootstrap_matrices <- vector("list", n_bootstrap)
     
     for (i in seq_len(n_bootstrap)) {
          sampled_participants <- vector("list", length(sample_sizes))
          sampled_contacts <- vector("list", length(sample_sizes))
          
          for (age_group in seq_along(sample_sizes)) {
               ids <- sample(data_p$part_id[data_p$age_group == age_group], size = sample_sizes[age_group], replace = TRUE)
               sampled_participants[[age_group]] <- data_p[data_p$part_id %in% ids,]
               sampled_contacts[[age_group]] <- data_c[data_c$part_id %in% ids,]
          }
          
          sampled_participants <- bind_rows(sampled_participants)
          sampled_contacts <- bind_rows(sampled_contacts)
          
          survey_results <- socialmixr::survey(sampled_participants, sampled_contacts)
          m <- socialmixr::contact_matrix(survey_results,
                                          age.limits = breaks,
                                          symmetric = FALSE,
                                          missing.contact.age = "sample")
          bootstrap_matrices[[i]] <- m$matrix
     }
     
     # Calculate the mean matrix
     mean_matrix <- Reduce(`+`, bootstrap_matrices) / length(bootstrap_matrices)
     return(mean_matrix)
}

# Execute the function
mean_contact_matrix <- calculate_matrices(participants, contacts, population, age_breaks)

# reshape contacts --------------------------------------------------------

data_mortality <- read.csv("./data/mortality_group.csv")
new_age_breaks <- unique(c(data_mortality$Age_start, data_mortality$Age_end))

reshape_contact_matrix <- function(original_matrix, old_age_breaks, new_age_breaks) {
     n <- length(new_age_breaks) - 1
     new_matrix <- matrix(0, nrow = n, ncol = n)
     
     for (i in 1:n) {
          for (j in 1:n) {
               new_start_i <- new_age_breaks[i]
               new_end_i <- new_age_breaks[i + 1]
               new_start_j <- new_age_breaks[j]
               new_end_j <- new_age_breaks[j + 1]
               
               for (k in 1:(length(old_age_breaks) - 1)) {
                    for (l in 1:(length(old_age_breaks) - 1)) {
                         old_start_k <- old_age_breaks[k]
                         old_end_k <- old_age_breaks[k + 1]
                         old_start_l <- old_age_breaks[l]
                         old_end_l <- old_age_breaks[l + 1]
                         
                         overlap_i <- max(0, min(new_end_i, old_end_k) - max(new_start_i, old_start_k))
                         overlap_j <- max(0, min(new_end_j, old_end_l) - max(new_start_j, old_start_l))
                         
                         if (overlap_i > 0 && overlap_j > 0) {
                              weight <- (overlap_i / (new_end_i - new_start_i)) * (overlap_j / (new_end_j - new_start_j))
                              new_matrix[i, j] <- new_matrix[i, j] + weight * original_matrix[k, l]
                         }
                    }
               }
          }
     }
     
     return(new_matrix)
}

adjusted_contact_matrix <- reshape_contact_matrix(
     original_matrix = mean_contact_matrix,
     old_age_breaks = age_breaks,
     new_age_breaks = new_age_breaks
)

# names matrix rows and columns
# 0-1m, 2-2m, 3-3m, 4-4m, 5-5m, 6-6m, 7-17m, 1.5-5y, 6-11y
# 12-20y, 21-30y, 31-40y, 41-50y, 51-60y, 61-70y, 71-100y
new_age_names <- c("[0m, 2m)", "[2m, 3m)", "[3m, 4m)", "[4m, 5m)", "[5m, 6m)", "[6m, 7m)", "[7m, 17m)", "[1.5y, 5y)",
                   "[6y, 11y)", "[12y, 20y)", "[21y, 30y)", "[31y, 40y)", "[41y, 50y)", "[51y, 60y)", "[61y, 70y)", "[71y, 100y)")

rownames(adjusted_contact_matrix) <- new_age_names
colnames(adjusted_contact_matrix) <- new_age_names

# visualize contacts -----------------------------------------------------

plot_contact_matrix <- function(contact_matrix, title = "Contact Matrix") {
     df <- as.data.frame(as.table(contact_matrix))
     names(df) <- c("age_of_participant", "age_of_contact", "frequency")
     
     # 将年龄组数据转为字符型以避免在绘图时被解释为数值
     df$age_of_participant <- as.character(df$age_of_participant)
     df$age_of_contact <- as.character(df$age_of_contact)
     
     p <- ggplot(df, aes(x = age_of_participant, y = age_of_contact, fill = frequency)) +
          geom_tile() +
          scale_fill_gradient(low = "white", high = "red", na.value = "grey50") +
          scale_x_discrete(limits = rownames(contact_matrix)) +
          scale_y_discrete(limits = rownames(contact_matrix)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                axis.title = element_blank()) +
          labs(fill = "Contacts per Day", title = title) +
          coord_fixed()
     
     return(p)
}

fig1 <- plot_contact_matrix(mean_contact_matrix, "Original Contact Matrix")
fig2 <- plot_contact_matrix(adjusted_contact_matrix, "Adjusted Contact Matrix")

# save-----------------------------------------------------

write.csv(adjusted_contact_matrix,
          "./data/contact_changsha_baseline.csv")

ggsave("./outcome/fig_s2_contact_matrix.png",
       fig1 + fig2,
       width = 12,
       height = 6)