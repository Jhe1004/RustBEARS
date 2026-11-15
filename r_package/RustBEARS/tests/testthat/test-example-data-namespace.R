test_that("the retained example belongs to the installed package namespace", {
    example <- new.env(parent=emptyenv())
    data("Psychotria_ML_DEC", package="RustBEARS", envir=example)
    fit <- example$Psychotria_ML_DEC
    expect_s3_class(fit, "calc_loglike_sp_results")
    expect_identical(attr(class(fit$inputs$BioGeoBEARS_model_object), "package"), "RustBEARS")
    expect_true(methods::validObject(fit$inputs$BioGeoBEARS_model_object))
    expect_true(is.finite(fit$total_loglikelihood))
})
