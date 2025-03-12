library(fastei)

tests <- 50
results <- vector("list", tests)
for (i in 1:tests) {
    results[[i]] <- fastei:::.random_samples(
        c(30, 50), # Ballots range
        c(2, 4), # Candidates range
        c(2, 4), # Demographic group range
        c(30, 60), # Votes per ballot range
        seed = i
    )
}

available_methods <- c("exact", "mult", "hnr", "mvn_cdf", "mvn_pdf")
method_division <- tests %/% 5
print(paste0("The method division is", method_division))
for (k in 1:5) {
    for (i in 1:method_division) {
        print(paste0("Going in with i=", i, "and k = ", k))
        print("Trying with X:")
        print(results[[i]]$X)
        print("Trying with W:")
        print(results[[i]]$W)
        print("The sum is:")
        print(sum(results[[i]]$X))
        result <- compute(X = results[[i]]$X, W = results[[i]]$W, method = available_methods[k])
        rm(result)
    }
}


object <- eim(samp$X, samp$W)
print("The object to pass is")
print(object)
compute(object = object, method = "exact")
print("\nrunnning the other one (directly from matrices)")
compute(X = samp$X, W = samp$W, method = "exact")
