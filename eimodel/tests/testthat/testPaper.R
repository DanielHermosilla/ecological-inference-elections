library(jsonlite)
library(dplyr)
library(ggplot2)
library(tidyr)

# Set the directory containing JSON files
json_dir <- "/Users/daniel/ecological-inference-elections/instances"

# Get a list of all JSON files in the directory
json_files <- list.files(json_dir, pattern = "\\.json$", full.names = TRUE)

# Function to compute results for a given stop_threshold

compute_results <- function(threshold) {
    results_df <- data.frame(G = integer(), I = integer(), AvgDifference = numeric(), AvgTime = numeric())

    differences_list <- list()
    time_list <- list()

    # Loop through each JSON file
    for (file in json_files) {
        filename <- basename(file)
        match <- regmatches(filename, regexec("G([0-9]+)_.*I([0-9]+)_.*seed([0-9]+)", filename))

        if (length(match[[1]]) == 4) {
            G <- as.integer(match[[1]][2])
            I <- as.integer(match[[1]][3])
            seed <- as.integer(match[[1]][4])

            # Load JSON data
            data <- fromJSON(file)
            cat("\n--- Processing File:", file, "---\n")
            print(str(data)) # Print JSON structure
            cat("\nThreshold:", threshold, "\n")

            # Check if matrices exist
            if (!("X" %in% names(data)) || !("W" %in% names(data)) || !("p" %in% names(data))) {
                stop("Error: JSON file does not contain required fields (X, W, p)")
            }

            # Check matrix structures
            print(dim(data$X))
            print(dim(data$W))
            print(dim(data$p))

            # Run EM algorithm
            result <- tryCatch(
                {
                    run_em(X = data$X, W = data$W, method = "mvn_cdf", stop_threshold = threshold)
                },
                error = function(e) {
                    cat("\nERROR in run_em():", conditionMessage(e), "\n")
                    return(NULL)
                }
            )

            if (!is.null(result)) {
                difference <- mean(abs(result$prob - data$p))
                time_taken <- result$time
                rm(result)

                key <- paste(G, I, sep = "_")
                if (!key %in% names(differences_list)) {
                    differences_list[[key]] <- c()
                    time_list[[key]] <- c()
                }
                differences_list[[key]] <- c(differences_list[[key]], difference)
                time_list[[key]] <- c(time_list[[key]], time_taken)
            }
        }
    }

    # Compute the average difference and time for each (G, I)
    for (key in names(differences_list)) {
        G_I_values <- as.integer(strsplit(key, "_")[[1]])
        G <- G_I_values[1]
        I <- G_I_values[2]

        avg_diff <- mean(differences_list[[key]], na.rm = TRUE)
        avg_time <- mean(time_list[[key]], na.rm = TRUE)

        results_df <- rbind(results_df, data.frame(G = G, I = I, AvgDifference = avg_diff, AvgTime = avg_time))
    }

    colnames(results_df)[colnames(results_df) == "I"] <- "C"
    results_df <- results_df[, c("C", "G", "AvgDifference", "AvgTime")]
    results_df <- results_df %>% arrange(C, G)

    return(results_df)
}
# Compute results for stop_threshold = 0.05 and 0.01
results_05 <- compute_results(0.05)
results_01 <- compute_results(0.01)

# Merge results on C and G
comparison_df <- merge(results_05, results_01, by = c("C", "G"), suffixes = c("_05", "_01"))

# Transform dataframe for plotting
comparison_long <- comparison_df %>%
    pivot_longer(cols = c(AvgDifference_05, AvgDifference_01), names_to = "Threshold", values_to = "AvgDifference") %>%
    mutate(Threshold = ifelse(Threshold == "AvgDifference_05", "0.05", "0.01"))

# Plot AvgDifference for stop_threshold = 0.05 vs. 0.01
p1 <- ggplot(comparison_long, aes(x = factor(C), y = AvgDifference, fill = Threshold)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        title = "Error (AvgDifference) Comparison",
        x = "C",
        y = "AvgDifference",
        fill = "Stop Threshold"
    ) +
    theme_minimal()

# Save the plot
ggsave("error_comparison.png", plot = p1, width = 7, height = 5, dpi = 300, bg = "white")

# Transform dataframe for time plotting
time_long <- comparison_df %>%
    pivot_longer(cols = c(AvgTime_05, AvgTime_01), names_to = "Threshold", values_to = "AvgTime") %>%
    mutate(Threshold = ifelse(Threshold == "AvgTime_05", "0.05", "0.01"))

# Plot AvgTime for stop_threshold = 0.05 vs. 0.01
p2 <- ggplot(time_long, aes(x = factor(C), y = AvgTime, fill = Threshold)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        title = "Computation Time Comparison",
        x = "C",
        y = "AvgTime (seconds)",
        fill = "Stop Threshold"
    ) +
    theme_minimal()

# Save the plot
ggsave("time_comparison.png", plot = p2, width = 7, height = 5, dpi = 300, bg = "white")
