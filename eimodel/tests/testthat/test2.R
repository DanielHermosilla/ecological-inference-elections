library(ggplot2)

# Define K values (logarithmic scale for better coverage)
K_values <- seq(0.1, 0.00001, length.out = 50)

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
# Store errors and times for each K
results <- data.frame(K = numeric(), Error = numeric(), Time = numeric(), Iterations = numeric())

# Loop through each K value
difference <- 0
total_time <- 0
total_iterations <- 0
for (i in 1:precisionTests) {
    model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
    result <- run_em(model, stop_threshold = 0.05, maxiter = 10000, method = "mvn_cdf")

    # Compute error
    difference <- difference + mean(abs(result$prob - precisionsResults[[i]]$prob))

    # Collect time taken
    total_time <- total_time + result$time
    total_iterations <- total_iterations + result$iterations

    rm(model)
    rm(result)
    gc()
}

error <- difference / precisionTests
avg_time <- total_time / precisionTests
avg_iterations <- total_iterations / precisionTests
cat("Error:\t", error, "\nTiempo:\t", avg_time, "\nIteraciones:\t", avg_iterations, "\n")

# Store results
results <- rbind(results, data.frame(K = k, Error = error, Time = avg_time, Iterations = avg_iterations))
# ---- 1. Plot Error vs K ----
p1 <- ggplot(results, aes(x = K, y = Error)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    # scale_y_log10() + # Log scale for Error
    labs(
        title = "MVN CDF",
        subtitle = "Error vs. Log-likelihood threshold",
        x = "Convergence Threshold",
        y = "Error (MAE)"
    ) +
    theme_minimal()

# Save the updated error plot
ggsave("error.png", plot = p1, width = 7, height = 5, dpi = 300, bg = "white")

# ---- 2. Plot Time vs K ----
p2 <- ggplot(results, aes(x = K, y = Time)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
        title = "MVN CDF",
        subtitle = "Tiempo vs. Log-likelihood threshold",
        x = "Convergence Threshold",
        y = "Tiempo (segundos)"
    ) +
    theme_minimal()

# Save the plot
ggsave("time.png", plot = p2, width = 7, height = 5, dpi = 300, bg = "white")

p3 <- ggplot(results, aes(x = K, y = Iterations)) +
    geom_line(color = "orange", size = 1) +
    geom_point(color = "orange", size = 2) +
    labs(
        title = "MVN CDF",
        subtitle = "Iteraciones vs. Log-likelihood threshold",
        x = "Convergence Threshold",
        y = "Iteraciones"
    ) +
    theme_minimal()

ggsave("iteraciones.png", plot = p3, width = 7, height = 5, dpi = 300, bg = "white")

# Print message confirming the files were saved
cat("Plots saved as 'error.png' and 'time.png' in your working directory.\n")
