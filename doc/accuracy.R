## ----lambda_comparison, fig.cap = "Mean Absolute Error vs Lambda", fig.align = "center", message=FALSE, warning=FALSE, results="hide"----
library(ggplot2)
library(fastei)
library(dplyr)
library(showtext)
font_add_google("Lato", "lato")
showtext_auto()

lambda_values <- seq(0, 1, by = 0.1)
methods <- c("exact", "mult", "mvn_cdf", "mvn_pdf", "mcmc")

# Prepare an empty list to hold all results
all_results <- list()

for (method in methods) {
    mae_values <- numeric(length(lambda_values))

    for (i in seq_along(lambda_values)) {
        lambda <- lambda_values[i]

        real <- simulate_election(num_ballots = 50, num_candidates = 2, num_groups = 3, seed = 42, lambda = lambda)

        estimate <- run_em(X = real$X, W = real$W, method = method)
        mae <- mean(abs(real$real_prob - estimate$prob))
        mae_values[i] <- mae
    }

    df_method <- data.frame(
        lambda = lambda_values,
        mae = mae_values,
        method = method
    )

    all_results[[method]] <- df_method
}

# Combine all method to one dataframe
df_all <- do.call(rbind, all_results)

# Rename the labels for the plot
label_lookup <- c(
    "exact" = "Exact",
    "mult" = "Multinomial",
    "mvn_cdf" = "MVN-CDF",
    "mvn_pdf" = "MVN-PDF",
    "mcmc" = "MCMC"
)

df_all <- df_all %>%
    mutate(method_label = factor(label_lookup[method], levels = unname(label_lookup)))

ggplot(df_all, aes(x = lambda, y = mae, color = method_label)) +
    geom_line(linewidth = 1.2, alpha = 0.4) +
    geom_point(aes(shape = method_label), size = 2.5) +
    labs(
        title = "Mean Absolute Error vs Lambda for Different Methods",
        subtitle = "Among 50 ballot boxes, 2 candidates and 3 groups",
        x = "Lambda",
        y = "Mean Absolute Error",
        color = "Method",
        shape = "Method"
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_minimal(base_size = 12, base_family = "lato")

