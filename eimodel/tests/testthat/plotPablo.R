library(ggplot2)

# Store q values from each iteration
q_list <- list()

precisionTests <- 50
precisionsResults <- vector("list", precisionTests)

# Generate random test cases
# for (i in 1:precisionTests) {
#    precisionsResults[[i]] <- fastei:::.random_samples(
#        c(200, 400), # Ballots range
#        c(2, 10), # Candidates range
#        c(2, 10), # Demographic group range
#        c(200, 400), # Votes per ballot range
#        seed = i
#    )
# }
for (i in 1:precisionTests) {
    precisionsResults[[i]] <- fastei::simulate_election(
        num_ballots = 50,
        num_candidates = 2,
        num_groups = 2,
        seed = i
    )
}
# Loop through each precision test
precisionTests <- 50
for (i in 1:precisionTests) {
    model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
    result <- run_em(model, stop_threshold = 0.000000000000000005, maxiter = 2500, method = "mult")

    # Store q array
    q_list[[i]] <- result$q

    rm(model, result)
    gc()
}

# Convert list of q arrays into a matrix (each row is an iteration)
q_matrix <- do.call(rbind, q_list) # Create matrix with rows as iterations

# Compute average across all iterations (column-wise mean)
q_mean <- colMeans(q_matrix, na.rm = TRUE) # Mean of each q element
q_mean <- c(NA, abs(diff(q_mean))) # NA para alinear índices
# Create dataframe for plotting
q_df <- data.frame(Index = seq_along(q_mean), AvgQ = q_mean)

# Plot the averaged q values
ggplot(q_df, aes(x = Index, y = AvgQ)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 1, alpha = 0.7) +
    scale_x_log10() +
    labs(
        title = "Delta Log-likelihood promedio por iteración",
        subtitle = "Multinomial",
        x = "Iteración",
        y = "Delta Log-likelihood"
    ) +
    coord_cartesian(ylim = c(0.0005, 0.11)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0.005, linetype = "dashed", color = "gray") +
    annotate("text", x = 5, y = 0.001, label = "0.005", vjust = -0.5, hjust = 0, color = "black") +
    annotate("text", x = 5, y = 0.01, label = "0.01", vjust = -0.5, hjust = 0, color = "black") +
    annotate("text", x = 5, y = 0.095, label = "0.1", vjust = -0.5, hjust = 0, color = "black") +
    theme_minimal()
# Save the plot
ggsave("q_values_average.png", width = 7, height = 5, dpi = 300, bg = "white")

cat("Plot saved as 'q_values_average.png'!\n")
