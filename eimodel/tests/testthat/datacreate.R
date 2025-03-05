electores <- read.csv("/Users/daniel/ecological-inference-elections/2021_11_Presidencial/output/2021_11_Presidencial_ELECTORES.csv") # nolint
votos <- read.csv("/Users/daniel/ecological-inference-elections/2021_11_Presidencial/output/2021_11_Presidencial_VOTOS.csv") # nolint

library(dplyr)

chile_election_2021 <- electores %>%
    inner_join(votos, by = c("CIRCUNSCRIPCION.ELECTORAL", "MESA"))

chile_election_2021 <- chile_election_2021 %>% select(-REGION.x, -LOCAL.x, -REGION.y, -LOCAL.y, -NULO.BLANCO)

chile_election_2021 <- chile_election_2021 %>% select(CIRCUNSCRIPCION.ELECTORAL, MESA, GABRIEL.BORIC, JOSE.ANTONIO.KAST, YASNA.PROVOSTE, SEBASTIAN.SICHEL, EDUARDO.ARTES, MARCO.ENRIQUEZ.OMINAMI, FRANCO.ALDO.PARISI, VOTOS.EN.BLANCO, VOTOS.NULOS, X18.19, X20.29, X30.39, X40.49, X50.59, X60.69, X70.79, X80.) # nolint

chile_election_2021 <- chile_election_2021 %>%
    rename(
        ELECTORAL.DISTRICT = CIRCUNSCRIPCION.ELECTORAL,
        BALLOT.BOX = MESA,
        C1 = GABRIEL.BORIC,
        C2 = JOSE.ANTONIO.KAST,
        C3 = YASNA.PROVOSTE,
        C4 = SEBASTIAN.SICHEL,
        C5 = EDUARDO.ARTES,
        C6 = MARCO.ENRIQUEZ.OMINAMI,
        C7 = FRANCO.ALDO.PARISI,
        BLANK.VOTES = VOTOS.EN.BLANCO,
        NULL.VOTES = VOTOS.NULOS
    )

chile_election_2021 <- chile_election_2021 %>%
    mutate(MISMATCH = rowSums(across(c("X18.19", "X20.29", "X30.39", "X40.49", "X50.59", "X60.69", "X70.79", "X80."))) !=
        rowSums(across(c("C1", "C2", "C3", "C4", "C5", "C6", "C7"))))

library(fastei)
# Load necessary library
library(dplyr)

# Filter rows where ELECTORAL.DISTRICT == "ANTOFAGASTA NORTE"
filtered_df <- chile_election_2021 %>%
    filter(ELECTORAL.DISTRICT == "NIEBLA")

# Create the X matrix with selected columns
X <- as.matrix(filtered_df[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7")])

# Create the W matrix with selected columns
W <- as.matrix(filtered_df[, c("X18.19", "X20.29", "X30.39", "X40.49", "X50.59", "X60.69", "X70.79", "X80.")])

# Print the matrices to verify
# print(X)
# print(W)

solution <- run_em(X = X, W = W, allow_mismatch = TRUE)
solution2 <- get_agg_proxy(X = X, W = W, allow_mismatch = TRUE, verbose = TRUE, sd_threshold = 0.0001, sd_statistic = "mean", seed = 1, nboot = 400)
