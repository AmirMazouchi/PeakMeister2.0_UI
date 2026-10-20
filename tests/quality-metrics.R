# Rscript --vanilla tests/quality-metrics.R
# Exercises the actual CV/PCA helpers and Calculate CVs observer. No raw mz5,
# FactoMineR, irlba, or Matrix package is required for this calculation.
suppressPackageStartupMessages({ library(shiny); library(plotly); library(jsonlite) })

is_assignment <- function(x, name) {
  is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
}
server <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(server[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
helper_names <- c("action_label", "report_action_error", "with_action_error", "observe_action",
                  "observe_updates", "extract_mass", "calculate_cv", "pareto_scale", "remove_zero_variance",
                  "perform_pca", "compute_quality_metrics", "build_quality_cv_plot", "build_quality_pca_plot")
helpers <- new.env(parent = globalenv())
for (name in helper_names) eval(assignment(name), helpers)
notifications <- character()
helpers$showNotification <- function(ui, ...) { notifications <<- c(notifications, ui); invisible(NULL) }
helpers$PCA <- function(...) stop("The failing FactoMineR/irlba backend must not be called")
same <- function(x, y, tolerance = 1e-9) stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance)))

# Independent reference: eigenvalues of the correlation matrix. This is the
# previous PCA(..., scale.unit=TRUE) normalization after Pareto scaling.
set.seed(42)
numeric_data <- as.data.frame(matrix(rnorm(12 * 20), nrow = 12))
numeric_data <- as.data.frame(sweep(numeric_data, 2, seq_len(20), "*") + 100)
colnames(numeric_data) <- paste0(seq(300, 110, by = -10), "_Metabolite")
pca <- helpers$perform_pca(numeric_data)
reference <- eigen(cor(numeric_data), symmetric = TRUE)
same(unname(pca$eig[, 1]), reference$values[1:5])
same(unname(pca$eig[, 2]), reference$values[1:5] / sum(reference$values) * 100)
centered <- scale(numeric_data, scale = FALSE)
standardized <- sweep(centered, 2, sqrt(colMeans(centered^2)), "/")
expected_coordinates <- standardized %*% reference$vectors[, 1:5]
same(unname(abs(pca$ind$coord)), unname(abs(expected_coordinates)))
stopifnot(identical(dim(pca$ind$coord), c(12L, 5L)))

quality_data <- data.frame(file.name = c("Run01", rep("", 11)), peak.number = seq_len(12),
  `PBM Sample ID` = c(paste0("QC_NIST_", 1:3), paste0("Patient", 1:9)),
  numeric_data, check.names = FALSE)
metrics <- helpers$compute_quality_metrics(quality_data, "QC=NIST", "Blank")
stopifnot(nrow(metrics$qc_samples) == 3L, nrow(metrics$non_qc_samples) == 9L,
          !is.null(metrics$pca_result))
expected_cv <- vapply(numeric_data[1:3, ], function(x) sd(x) / mean(x) * 100, numeric(1))
same(metrics$qc_cv_values, expected_cv)

# Sorting and column exclusion must not attach a QC value to another metabolite.
pruned <- quality_data; pruned[[colnames(numeric_data)[2]]][1:3] <- NA_real_
pruned_metrics <- helpers$compute_quality_metrics(pruned, "QC=NIST", "")
stopifnot(!colnames(numeric_data)[2] %in% names(pruned_metrics$qc_cv_values))
plot <- plotly_build(helpers$build_quality_cv_plot(pruned_metrics))
qc_trace <- Filter(function(x) identical(x$name, "QC CV"), plot$x$data)[[1]]
expected_order <- order(as.numeric(sub("_.*", "", names(pruned_metrics$qc_cv_values))))
same(as.numeric(qc_trace$y), unname(pruned_metrics$qc_cv_values[expected_order]))
stopifnot(length(qc_trace$x) == 19L)
warnings <- character()
withCallingHandlers({
  invisible(ggplot2::ggplot_build(helpers$build_quality_pca_plot(metrics)))
}, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
stopifnot(!any(grepl("Too few|singular|positive definite", warnings, ignore.case = TRUE)))

# Invalid inputs return a readable notification, never an uncaught exception.
reject <- function(data, qc = "QC=NIST", excluded = "") {
  before <- length(notifications)
  stopifnot(is.null(helpers$compute_quality_metrics(data, qc, excluded)), length(notifications) > before)
}
reject(quality_data, "QC="); reject(quality_data, ""); reject(quality_data, NULL)
reject(quality_data, "QC=Unknown"); reject(quality_data, "QC=."); reject(quality_data, "QC=[")
reject(quality_data, excluded = paste0("QC_NIST_", 1:3, collapse = ","))
reject(quality_data[, names(quality_data) != "PBM Sample ID"])
missing_ids <- quality_data; missing_ids$`PBM Sample ID` <- NA_character_; reject(missing_ids)
missing_value <- quality_data; missing_value[1, 4] <- NA_real_; reject(missing_value)
bad_value <- quality_data; bad_value[1, 4] <- "<LOD"; reject(bad_value)
infinite <- quality_data; infinite[1, 4] <- Inf; reject(infinite)
empty <- quality_data[0, ]; reject(empty); reject(NULL)
unusable <- quality_data; unusable[, 4:ncol(unusable)] <- NA_real_; reject(unusable)

# Too few features/rank prevents a 2D PCA, while valid CVs remain usable.
single <- quality_data[, 1:4]
single_metrics <- helpers$compute_quality_metrics(single, "QC=NIST", "")
stopifnot(is.null(single_metrics$pca_result), length(single_metrics$qc_cv_values) == 1L,
          is.finite(single_metrics$qc_cv))
constant <- quality_data; constant[, 4:ncol(constant)] <- 2
constant_metrics <- helpers$compute_quality_metrics(constant, "QC=NIST", "")
stopifnot(is.null(constant_metrics$pca_result), constant_metrics$qc_cv == 0)
collinear <- quality_data[, 1:5]; collinear[, 4] <- seq_len(12); collinear[, 5] <- 2 * seq_len(12)
collinear_metrics <- helpers$compute_quality_metrics(collinear, "QC=NIST", "")
stopifnot(is.null(collinear_metrics$pca_result), is.finite(collinear_metrics$qc_cv))
stopifnot(is.na(helpers$calculate_cv(data.frame(x = c(-1, 1)))[1]))
stopifnot(is.null(helpers$perform_pca(numeric_data[1:2, ])),
          is.null(helpers$perform_pca(data.frame(x = c("a", "b", "c")))))

# Exercise actual Shiny observer clicks, including failure and recovery. The
# artificial failure uses the exact NULL-object error reported by the user.
observer <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) &&
                     identical(x[[2]], quote(input$compute_cv)), statements)[[1]]
result_folder <- tempfile("quality-metrics-"); dir.create(file.path(result_folder, "Data"), recursive = TRUE)
testServer(function(input, output, session) {
  # App dependencies can mask shiny::validate (jsonlite is loaded after shiny
  # in ui.R). Reproduce the reported S3 validate(NULL) failure as well.
  validate <- function(object, ...) UseMethod("validate")
  for (name in helper_names) eval(assignment(name), envir = environment())
  showNotification <- helpers$showNotification
  peak_areas_data <- reactiveVal(quality_data)
  main_folder2 <- function() result_folder
  force_failure <- reactiveVal(FALSE)
  real_scaling <- pareto_scale
  pareto_scale <- function(data) {
    if (force_failure()) stop("Error in oknum: LENGTH or similar applied to NULL object")
    real_scaling(data)
  }
  eval(assignment("average_cv"), envir = environment())
  eval(observer, envir = environment())
}, {
  session$setInputs(category_input = "QC=NIST", exclude_ids = "Blank", compute_cv = 1)
  stopifnot(grepl("Technical CV", output$average_cv_table), !session$isClosed())
  invisible(output$cv_plot); invisible(output$pca_plot)
  stopifnot(file.exists(file.path(result_folder, "Data", "QC_PCA_plot.png")))
  force_failure(TRUE)
  session$setInputs(compute_cv = 2)
  stopifnot(grepl("Technical CV", output$average_cv_table), is.null(average_cv()$pca_result), !session$isClosed())
  invisible(output$cv_plot)
  stopifnot(any(grepl("PCA unavailable:.*oknum", notifications)))
  pca_output <- try(output$pca_plot, silent = TRUE)
  stopifnot(inherits(pca_output, "try-error"), grepl("CV results are shown separately", pca_output))
  session$setInputs(category_input = "QC=", compute_cv = 3)
  stopifnot(!session$isClosed(), is.null(average_cv()))
  force_failure(FALSE)
  session$setInputs(category_input = "QC=NIST", compute_cv = 4)
  stopifnot(!is.null(average_cv()$pca_result), grepl("Technical CV", output$average_cv_table), !session$isClosed())
  invisible(output$pca_plot)
})
unlink(result_folder, recursive = TRUE)
cat("PASS: native PCA matches correlation eigenvalues/scores; QC/input validation; single/constant/rank-deficient data keep CVs; metabolite-aligned plots; actual Calculate CVs observer survives the reported NULL error and recovers.\n")
