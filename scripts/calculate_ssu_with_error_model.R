#!/usr/bin/env Rscript

# ============================================================
# Libraries
# ============================================================
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "glue", "tidyverse", "data.table", "vroom", "gtools", "parallel")
invisible(lapply(packages, quiet_library))

# ============================================================
# Model description
# ============================================================
model_help_description <- glue(r"(
SSU Error Model
================

Goal:
    Estimate a reliable SSU for each variant and base position by
    accounting for sequencing/sampling noise and replicate-specific
    technical variation.

Input:
    Base coverage and total coverage for each variant, base position,
    and replicate.

1. Bayesian SSU estimation
    └─ Add eps = 0.5 to avoid 0/1 estimates:

       SSU_adj = (base_cov + 0.5) / (total_cov + 1)

       Equivalent to the posterior mean under a
       Beta(0.5, 0.5) Jeffreys prior.
    │
    ▼
2. Logit transformation
    └─ Transform SSU to an unbounded scale:

       theta = logit(SSU_adj)
             = log(SSU_adj / (1 - SSU_adj))

       This provides an approximately Gaussian scale for the
       error model.
    │
    ▼
3. Sampling / multiplicative variance
    └─ Estimate variance caused by finite sequencing depth:

       var_mult = 1 /
                  (total_cov * SSU_adj * (1 - SSU_adj))

       Higher coverage → lower sampling variance.
       Lower coverage  → higher sampling variance.
    │
    ▼
4. Collapse identical likelihood patterns
    └─ For computational efficiency, base positions are collapsed
       only when their complete likelihood information is identical:

       var_id
       + logit(SSU) for every replicate
       + sampling variance for every replicate
       + replicate missing/valid pattern

       n_obs records how many original base positions are represented.

       This preserves the original likelihood exactly while reducing
       the number of rows used during optimisation.
    │
    ▼
5. Replicate-subset likelihood
    └─ Construct all replicate combinations containing ≥2 replicates.

       For 3 replicates:
           (1,2)
           (1,3)
           (2,3)
           (1,2,3)

       Subsets allow replicate-specific technical variances to
       be identified.
    │
    ▼
6. Estimate replicate-specific additive variance
    └─ Each replicate has its own technical variance:

       var_addi[r] = sigma²_rep,r

       Total variance:

       var_total(v,r) =
           var_mult(v,r) + var_addi[r]

       The additive variances are estimated jointly by minimising
       the profiled Gaussian negative log-likelihood across all
       replicate subsets.

       Optimisation is performed on log(sigma²_rep) to ensure
       positive variances.
    │
    ▼
7. Error-corrected SSU
    └─ Return to the original base-level observations and combine
       replicates using inverse-variance weighting:

       theta_hat(v) =
           Σ [theta(v,r) / var_total(v,r)]
           --------------------------------
           Σ [1 / var_total(v,r)]

       var_theta(v) =
           1 / Σ [1 / var_total(v,r)]

       Replicates with lower total variance receive greater weight.

       At least 2 valid replicates are required.
    │
    ▼
8. Empirical Bayes shrinkage
    └─ Shrink uncertain estimates toward the global mean:

       mu_global = mean(theta_hat)
       tau²       = var(theta_hat)

       lambda(v) =
           tau² / (tau² + var_theta(v))

       theta_shrunk(v) =
           lambda(v) * theta_hat(v)
           + (1 - lambda(v)) * mu_global

       Low-precision estimates → stronger shrinkage.
       High-precision estimates → weaker shrinkage.
    │
    ▼
9. Final corrected SSU
    └─ Transform back to the [0,1] scale:

       SSU_corrected = plogis(theta_shrunk)
    │
    ▼
10. 95% confidence interval
    └─ Calculate the interval on the logit scale and transform
       back to SSU:

       theta_shrunk ± 1.96 * sqrt(var_theta_shrunk)

       SSU_lwr = plogis(lower limit)
       SSU_upr = plogis(upper limit)

Output:
    One row per variant + base position containing the corrected
    SSU, uncertainty, number of replicates used, replicate-level
    SSU/coverage, and 95% confidence interval.

Computational optimisation:
    Identical likelihood patterns are collapsed before optimisation,
    and the replicate-subset likelihood is evaluated in parallel.
    The final corrected SSU is always calculated at the original
    variant + base-position level.
)")

# ============================================================
# Options
# ============================================================
option_list <- list(make_option(c("-r", "--rscript_dir"),     type = "character",    help = "directory path of R scripts", default = NULL),
                    make_option(c("-s", "--sample_id"),       type = "character",    help = "list of sample IDs",          default = NULL),
                    make_option(c("-d", "--ssu_counts"),      type = "character",    help = "list of SSU counts",          default = NULL),
                    make_option(c("-o", "--output_dir"),      type = "character",    help = "output directory",            default = getwd()),
                    make_option(c("-p", "--prefix"),          type = "character",    help = "output prefix",               default = "sample"),
                    make_option(c("-m", "--model_help"),      action = "store_true", help = "print model description",     default=FALSE))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
    print_help(opt_parser)
    quit(status = 1)
}

if(opt$model_help)
{
    print(model_help_description)
    quit(status = 0)
}

# ============================================================
# Check options
# ============================================================
if(is.null(opt$rscript_dir)) stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$sample_id))   stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$ssu_counts))  stop("-d, list of splicing counts is required!", call. = FALSE)

# ============================================================
# Modules
# ============================================================
source(file.path(opt$rscript_dir, "report_utils.R"))

# ============================================================
# Inputs
# ============================================================
sample_reps      <- unlist(strsplit(opt$sample_id, ","))
files_ssu_counts <- unlist(strsplit(opt$ssu_counts, ","))

sample_reps      <- mixedsort(sample_reps)
files_ssu_counts <- sort_paths_by_filename(files_ssu_counts)

# ============================================================
# Outputs
# ============================================================
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- paste0(opt$prefix, ".ssu_per_base")

# ============================================================
# 1. Read input files
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "1. reading input files ...")

ssu_counts <- list()
for(i in seq_along(sample_reps)) {    
    ssu_counts[[sample_reps[i]]] <- as.data.table(
        vroom(
            files_ssu_counts[i], 
            delim = "\t", 
            comment = "#", 
            col_names = TRUE, 
            show_col_types = FALSE
        )
    )
}

dt_ssu <- rbindlist(ssu_counts, idcol = "reps")

rm(ssu_counts)
invisible(gc())

# ============================================================
# 2. Calculate SSU with Bayesian prior
# Note: For low-count variants, SSU is pulled toward 0.5. 
#       This is intentional Bayesian shrinkage.
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "2. calculate SSU with eps ...")
eps <- 0.5
dt_ssu[, ssu_eps := NA_real_] # avoid value issue when max_cov = 0
dt_ssu[max_cov > 0, ssu_eps := (base_cov + eps) / (max_cov + 2 * eps)]

# ============================================================
# 3. Logit SSU and multiplicative variance
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "3. calculate logit SSU and multiplicative variance ...")
dt_ssu[, logit_ssu := NA_real_]
dt_ssu[max_cov > 0, logit_ssu := qlogis(ssu_eps)]
 
dt_ssu[, var_mult := NA_real_]
dt_ssu[max_cov > 0, var_mult := 1 / (max_cov * ssu_eps * (1 - ssu_eps))]
 
# ============================================================
# 4. Build wide likelihood table
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "4. build likelihood table ...")
dt_ssu_wide <- dcast(dt_ssu, var_id + base_pos ~ reps, value.var = c("logit_ssu", "var_mult"))
setorder(dt_ssu_wide, var_id, base_pos)

n_rows_original <- nrow(dt_ssu_wide)
n_reps <- length(sample_reps)

cols_logit_ssu <- paste0("logit_ssu_", sample_reps)
cols_var_mult <- paste0("var_mult_",  sample_reps)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> original likelihood rows: ", format(n_rows_original, big.mark = ","))
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> replicates: ", n_reps)

# ============================================================
# 5. Collapse identical likelihood patterns to speed up
# ============================================================
# Collapse rows only when ALL likelihood-relevant quantities are
# identical across the complete replicate vector.
#
# Input table:
#   var_id   base_pos   logit_ssu_1   var_mult_1   logit_ssu_2   var_mult_2   ...
#   VAR1     1          2.31          0.12         2.28          0.10
#   VAR1     2          2.31          0.12         2.28          0.10
#   VAR1     3          2.31          0.12         2.28          0.10
#   ...
#
# Rows are collapsed only when the following columns are identical:
#
#   var_id
#   logit_ssu_1, var_mult_1
#   logit_ssu_2, var_mult_2
#   ...
#   logit_ssu_R, var_mult_R
#
# NA values are also part of the grouping, so the replicate
# validity/missingness pattern is preserved automatically.
#
# Example:
#
# Before collapsing:
#   var_id   base_pos   logit_ssu_1   var_mult_1   logit_ssu_2   var_mult_2
#   VAR1     1          0             0.01         0             0.02
#   VAR1     2          0             0.01         0             0.02
#   VAR1     3          0             0.01         0             0.02
#   VAR1     4          0             0.01         0             0.02
#   VAR1     5          0             0.01         0             0.02
#   VAR1     6          0             0.01         NA            NA
#   VAR1     7          0             0.01         NA            NA
#
# After collapsing:
#   var_id   logit_ssu_1   var_mult_1   logit_ssu_2   var_mult_2   n_obs
#   VAR1     0             0.01         0             0.02         5
#   VAR1     0             0.01         NA            NA            2
#
# The first row represents base positions 1-5, where both replicates
# are valid. The second row represents positions 6-7, where replicate 2
# is missing.
#
# Collapsing across the complete replicate vector therefore preserves
# the original likelihood exactly, while reducing the number of rows
# that need to be evaluated.
# ============================================================

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "5. collapse identical likelihood patterns ...")

likelihood_cols <- c(cols_logit_ssu, cols_var_mult)

dt_ssu_collapsed <- dt_ssu_wide[, .(n_obs = .N), by = c("var_id", likelihood_cols)]
n_rows_collapsed <- nrow(dt_ssu_collapsed)

collapsed_ratio <- 1 - n_rows_collapsed / n_rows_original

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> rows before being callapsed: ", format(n_rows_original, big.mark = ","))
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> rows after being callapsed: ", format(n_rows_collapsed, big.mark = ","))
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> collapsed ratio: ", signif(collapsed_ratio, 4))
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> total n_obs: ", format(sum(dt_ssu_collapsed$n_obs), big.mark = ","))

# ============================================================
# 6. Convert collapsed table to matrices
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "6. build collapsed likelihood matrices ...")

mat_logit_ssu_collapsed <- as.matrix(dt_ssu_collapsed[, ..cols_logit_ssu])
mat_var_mult_collapsed <- as.matrix(dt_ssu_collapsed[, ..cols_var_mult])

check_val_collapsed <- is.finite(mat_logit_ssu_collapsed) & is.finite(mat_var_mult_collapsed)

# Invalid values are not used in the likelihood.
mat_logit_ssu_collapsed[!check_val_collapsed] <- 0
mat_var_mult_collapsed[!check_val_collapsed] <- 0

n_rows_collapsed <- nrow(mat_logit_ssu_collapsed)
n_obs <- dt_ssu_collapsed$n_obs

# ============================================================
# 7. Build replicate subsets
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "7. build replicate subsets ...")
rep_subsets <- unlist(
    lapply(
        2:n_reps,
        function(k) {
            combn(
                seq_len(n_reps),
                k,
                simplify = FALSE
            )
        }
    ),
    recursive = FALSE
)

n_subsets <- length(rep_subsets)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
    "    |----> number of replicate subsets: ", n_subsets
)

# ============================================================
# 8. Precompute subset data
# ============================================================
# Everything here is independent of the optimisation parameters.
# Do this ONCE rather than repeatedly inside optim().
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "8. precompute subset matrices to speed up calculation ...")

subset_data <- lapply(
    rep_subsets,
    function(idx) {
        mat_logit_ssu_sub <- mat_logit_ssu_collapsed[, idx, drop = FALSE]
        mat_var_mult_sub <- mat_var_mult_collapsed[, idx, drop = FALSE]
        check_val_sub <- check_val_collapsed[, idx, drop = FALSE]

        # A collapsed pattern contributes to the subset likelihood
        # only if it has >=2 valid replicate observations.
        # This corresponds to the original row-wise likelihood.
        keep <- rowSums(check_val_sub) >= 2

        list(idx   = idx,
             mat_L = mat_logit_ssu_sub,
             mat_V = mat_var_mult_sub,
             check = check_val_sub,
             n_obs = n_obs,
             keep  = keep)
    }
)

# ============================================================
# 9. Fast subset likelihood
# ============================================================
# For each collapsed row:
#     var_total_r = var_mult_r + var_addi_r
#
# Weighted theta:
#     theta =
#       sum(n_obs * logit_ssu / var_total) /
#       sum(n_obs / var_total)
#
# Profiled NLL:
#     0.5 * sum( n_obs * [ log(var_total) + logit_ssu² / var_total ] )
#
# after profiling theta:
#     NLL = 0.5 * sum( n_obs * [ log(var_total) + C - B²/A ] )
#
# where:
#     A = sum(1 / var_total)
#     B = sum(logit_ssu / var_total)
#     C = sum(logit_ssu² / var_total)
#
# Because n_obs is included in all sufficient statistics,
# this is equivalent to using every original base position.
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "9. build subset likelihood estimation ...")

subset_nll <- function(log_var_addi_sub, data) {
    var_addi_sub <- exp(log_var_addi_sub)

    mat_logit_ssu_sub <- data$mat_L
    mat_var_mult_sub  <- data$mat_V
    check_val_sub     <- data$check
    n_obs             <- data$n_obs
    keep              <- data$keep

    n_rows <- nrow(mat_var_mult_sub)

    # --------------------------------------------------------
    # Add replicate-specific additive variance.
    # --------------------------------------------------------
    mat_var_total_sub <- mat_var_mult_sub + rep(var_addi_sub, each = n_rows)
    mat_var_total_sub_inv <- 1 / mat_var_total_sub

    # --------------------------------------------------------
    # Weighted sufficient statistics
    # --------------------------------------------------------
    weight <- n_obs * check_val_sub

    A <- rowSums(weight * mat_var_total_sub_inv)
    B <- rowSums(weight * mat_logit_ssu_sub * mat_var_total_sub_inv)
    C <- rowSums(weight * mat_logit_ssu_sub * mat_logit_ssu_sub * mat_var_total_sub_inv)
    L <- rowSums(weight * log(mat_var_total_sub))

    # --------------------------------------------------------
    # Keep only patterns with >=2 valid replicates.
    # --------------------------------------------------------
    A_keep <- A[keep]
    B_keep <- B[keep]
    C_keep <- C[keep]
    L_keep <- L[keep]

    # --------------------------------------------------------
    # Profiled Gaussian NLL
    # --------------------------------------------------------
    sum( 0.5 * (L_keep + C_keep - B_keep * B_keep / A_keep) )
}

# ============================================================
# 10. Parallel workers and joint likelihood function
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "10. start parallel workers and build joint likelihood function ...")
n_cores <- 4

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),"    |----> using ", n_cores, " CPU cores (3 replicates)")
cl <- makeCluster(n_cores)

clusterExport(
    cl,
    varlist = c("subset_data", "subset_nll"),
    envir = environment()
)

# Assign subsets to workers with 3 replicates:
#   worker 1 -> subset 1
#   worker 2 -> subset 2
#   worker 3 -> subset 3
#   worker 4 -> subset 4

chunks <- split(
    seq_along(subset_data),
    rep(seq_len(n_cores), length.out = n_subsets)
)

joint_nll <- function(log_var_addi) {
    partials <- clusterApply(
        cl,
        chunks,
        function(idxs, log_var_addi) {
            var_total <- 0
            for (i in idxs) {
                data <- subset_data[[i]]
                log_var_addi_sub <- log_var_addi[data$idx]
                var_total <- var_total + subset_nll(log_var_addi_sub,data)
            }
            var_total
        },
        log_var_addi = log_var_addi
    )
    sum(unlist(partials, use.names = FALSE))
}

# ============================================================
# 11. Estimate additive replicate variances
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "11. estimate additive replicate variances ...")

init <- rep(log(0.01), n_reps)
names(init) <- sample_reps

fit <- optim(
    par = init,
    fn = joint_nll,
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e7)
)

stopCluster(cl)

var_addi_est <- exp(fit$par)
names(var_addi_est) <- sample_reps

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> optimisation convergence: ", fit$convergence)
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> function evaluations: ", fit$counts[["function"]])
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "    |----> estimated additive variances: ", fit$counts[["function"]])
for (srep in sample_reps) {
    message(
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
        "    |----> ", srep, ": ", format(var_addi_est[srep], scientific = TRUE)
    )
}

# ============================================================
# 12. Error-corrected SSU
# ============================================================
# Return to the ORIGINAL base-level dt_ssu.
# Final SSU remains one result per var_id + base_pos
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "12. calculate error-corrected SSU ...")
dt_ssu[, var_addi := var_addi_est[reps]]
dt_ssu[, var_total := var_mult + var_addi]

dt_ssu_corrected <- dt_ssu[
    ,
    {
        check_val <- is.finite(logit_ssu) & is.finite(var_total)
        n_check_val <- sum(check_val)

        if (n_check_val >= 2) {
            weight      <- 1 / var_total[check_val]
            denominator <- sum(weight)
            theta       <- sum(logit_ssu[check_val] * weight) / denominator
            var_theta   <- 1 / denominator
        } else {
            theta <- NA_real_
            var_theta <- NA_real_
        }

        .(
            theta       = theta,
            var_theta   = var_theta,
            n_reps_used = n_check_val
        )
    },
    by = .(var_id, base_pos)
]

dt_ssu_corrected[, ssu_est := plogis(theta)]

dt_ssu_wide <- dcast(dt_ssu, var_id + base_pos ~ reps, value.var = c("base_ssu", "max_cov"))

cols_base_ssu <- paste0("base_ssu_", sample_reps)
cols_max_cov <- paste0("max_cov_", sample_reps)

setnames(dt_ssu_wide, cols_base_ssu, paste0("ssu", seq_along(cols_base_ssu)))
setnames(dt_ssu_wide, cols_max_cov, paste0("mcov", seq_along(cols_max_cov)))

dt_ssu_corrected <- merge(dt_ssu_wide, dt_ssu_corrected, by = c("var_id", "base_pos"), all.x = TRUE)

rm(dt_ssu_wide)
invisible(gc())

# ============================================================
# 13. Empirical Bayes shrinkage
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "13. shrinkage (empirical Bayes) ...")

mu_global <- mean(dt_ssu_corrected$theta, na.rm = TRUE)
tau2 <- var(dt_ssu_corrected$theta, na.rm = TRUE)

# for each variant per pos: λ(v,p) = tau2 / (tau2 + var_theta(v,p))
dt_ssu_corrected[, shrinkage := tau2 / (tau2 + var_theta)]

# θ(shrunk)​ = λ(v,p)​θ(v,p)​ + (1−λ(v,p)​) * mu_global
dt_ssu_corrected[, theta_shrunk := shrinkage * theta + (1 - shrinkage) * mu_global]
dt_ssu_corrected[, var_theta_shrunk := shrinkage^2 * var_theta]
dt_ssu_corrected[, ssu_corrected := plogis(theta_shrunk)]

# ============================================================
# 14. 95% confidence intervals
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "14. calculate confidence intervals ...")

# θ(hat​) = logit(SSU) ∼ N(θ,Var(θ(est)​))
# a standard normal variable Z ~ N(0,1), so P(|Z| <= 1.96) = 0.95
# then 95% CI = mean ± 1.96 × SD
z <- 1.96
dt_ssu_corrected[, ssu_corrected_lwr := plogis(theta_shrunk - z * sqrt(var_theta_shrunk))]
dt_ssu_corrected[, ssu_corrected_upr := plogis(theta_shrunk + z * sqrt(var_theta_shrunk))]
dt_ssu_corrected[, ssu_ci_width := ssu_corrected_upr - ssu_corrected_lwr]

dt_ssu_corrected[, ssu_precision_class := fcase(
    ssu_ci_width <= 0.10, "high",
    ssu_ci_width <= 0.20, "moderate",
    default = "low"
)]

# ============================================================
# 15. Output
# ============================================================
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "15. output ...")
num_cols <- names(dt_ssu_corrected)[vapply(dt_ssu_corrected, is.numeric, logical(1))]
dt_ssu_corrected[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]

output_file <- file.path(opt$output_dir, paste0(sample_prefix, ".details.tsv"))
fwrite(dt_ssu_corrected, file = output_file, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE)

message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "All done ...")
