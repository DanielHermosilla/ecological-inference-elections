## -----------------------------------------------------------------------------
library(fastei)

el_golf <- get_XW_chile(elect_district = "EL GOLF")

el_golf

## -----------------------------------------------------------------------------
results <- run_em(el_golf)
results$prob

## -----------------------------------------------------------------------------
results <- bootstrap(results, seed = 42, nboot = 30)
results$sd

## -----------------------------------------------------------------------------
navidad <- get_XW_chile(elect_district = "NAVIDAD")
navidad <- bootstrap(navidad, seed = 42, nboot = 30)
navidad$sd

## -----------------------------------------------------------------------------
navidad_proxy <- get_agg_proxy(navidad, seed = 42)
navidad_proxy$group_agg

## -----------------------------------------------------------------------------
mean(navidad$sd) - mean(navidad_proxy$sd)

## -----------------------------------------------------------------------------
navidad_opt <- get_agg_opt(navidad, seed = 42)
navidad_opt$group_agg

## -----------------------------------------------------------------------------
apoquindo <- get_XW_chile(elect_district = "APOQUINDO")

comparison <- welchtest(
    X = el_golf$X,
    W = el_golf$W,
    X2 = apoquindo$X,
    W2 = apoquindo$W,
    method = "mult",
    nboot = 30,
    seed = 42
)

nonsignificant <- comparison$pvals >= 0.05
nonsignificant

## -----------------------------------------------------------------------------
samples <- simulate_election(num_groups = 3, num_candidates = 3, num_ballots = 50, seed = 42)
object <- eim(samples$X, samples$W)

## -----------------------------------------------------------------------------
multinomial <- run_em(object, method = "mult")
multinomial_mae <- mean(abs(multinomial$prob - samples$real_prob))
multinomial_mae
multinomial$time

## -----------------------------------------------------------------------------
multivariate_pdf <- run_em(object, method = "mvn_pdf")
multivariate_pdf_mae <- mean(abs(multivariate_pdf$prob - samples$real_prob))
multivariate_pdf_mae
multivariate_pdf$time

## -----------------------------------------------------------------------------
multivariate_cdf <- run_em(object, method = "mvn_cdf", mc_error = 1e-8)
multivariate_cdf_mae <- mean(abs(multivariate_cdf$prob - samples$real_prob))
multivariate_cdf_mae
multivariate_cdf$time

