library(jsonlite)
library(ggplot2)

# Set the directory containing JSON files
json_dir <- "/Users/daniel/ecological-inference-elections/instances"

# Get a list of all JSON files containing "G2" and "I2"
json_files <- list.files(json_dir, pattern = "G2.*I3.*\\.json$", full.names = TRUE)

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
        seed <- 1

        # Load JSON data
        data <- fromJSON(file)

        # Ensure required fields exist
        if (!("X" %in% names(data)) || !("W" %in% names(data)) || !("p" %in% names(data))) {
            next # Skip file if missing data
        }

        # Run EM algorithm
        result <- run_em(X = data$X, W = data$W, verbose = FALSE, method = "hnr", maxiter = 1000, stop_threshold = 0.0001, p_method = "group_proportional")
        # result <- get_agg_proxy(object = result, verbose = TRUE, sd_threshold = 0.001)
        print("The real probability is:")
        print(data$p)
        # Store the q array
        q_list[[seed]] <- result$q
    }
}
