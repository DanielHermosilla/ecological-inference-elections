library(jsonlite)
library(ggplot2)

# Set the directory containing JSON files
json_dir <- "/Users/daniel/ecological-inference-elections/instances"

# Get a list of all JSON files containing "G2" and "I2"
json_files <- list.files(json_dir, pattern = "G5.*I2.*\\.json$", full.names = TRUE)

# Initialize a list to store all q arrays
q_list <- list()

# Loop through each filtered JSON file
for (file in json_files) {
    filename <- basename(file)
    match <- regmatches(filename, regexec("G([0-9]+)_.*I([0-9]+)_.*seed([0-9]+)", filename))

    if (length(match[[1]]) == 4) {
        G <- as.integer(match[[1]][2])
        I <- as.integer(match[[1]][3])
        seed <- as.integer(match[[1]][4])

        # Load JSON data
        data <- fromJSON(file)

        # Ensure required fields exist
        if (!("X" %in% names(data)) || !("W" %in% names(data)) || !("p" %in% names(data))) {
            next # Skip file if missing data
        }

        # Run EM algorithm
        result <- run_em(X = data$X, W = data$W, method = "mult", stop_threshold = 0.001)
        result <- get_agg_proxy(object = result, verbose = TRUE, sd_threshold = 0.001)

        # Store the q array
        q_list[[seed]] <- result$q
    }
}

# Convert the list to a dataframe for plotting
q_df <- do.call(cbind, q_list)
q_df <- as.data.frame(q_df)

# Add iteration index
q_df$Iteration <- 1:nrow(q_df)

# Convert to long format for ggplot
q_long <- tidyr::pivot_longer(q_df, cols = -Iteration, names_to = "Seed", values_to = "Q")

# Plot all 20 curves
ggplot(q_long, aes(x = Iteration, y = Q, color = Seed, group = Seed)) +
    geom_line(alpha = 0.5) +
    labs(
        title = "Log-likelihood con MVN CDF (C2, G2)",
        subtitle = "Instancias del paper",
        x = "Iteración",
        y = "Log-likelihood",
        color = "Seed"
    ) +
    scale_x_log10() + # Log scale for better visualization
    theme_minimal()

# Save the plot
ggsave("q_values_G2_I2.png", width = 7, height = 5, dpi = 300, bg = "white")
