test_that("get_eim_chile loads district data", {
    eim_obj <- get_eim_chile(elect_district = "APOQUINDO", remove_mismatch = FALSE)

    expect_s3_class(eim_obj, "eim")
    expect_true(is.matrix(eim_obj$X))
    expect_true(is.matrix(eim_obj$W))
    expect_equal(nrow(eim_obj$X), nrow(eim_obj$W))
    expect_true(ncol(eim_obj$X) >= 2)
    expect_true(ncol(eim_obj$W) >= 2)
})
