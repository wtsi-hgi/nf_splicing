#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "glue", "tidyverse", "data.table", "vroom", "gtools", "parallel")
invisible(lapply(packages, quiet_library))

model_help_description <- glue(r"(
1. Raw read count per base per variant per replicate
    │
    ▼
2. Calculate SSU with Bayesian prior (eps = 0.5)
    └─ SSU_adj = (base_cov + eps) / (total_cov + 2 * eps)
    │
    ▼
3. Logit transformation
    └─ logit(SSU) = log(SSU_adj / (1 - SSU_adj))
    │
    ▼
4. Compute multiplicative variance (sampling noise)
    └─ var_mult = 1 / (n_total * SSU_adj * (1 - SSU_adj))
    │
    ▼
5. Replicate subsets
    └─ Build all combinations of replicates
    │
    ▼
6. Estimate additive replicate variance (σ²_rep)
    └─ Optimize joint negative log-likelihood across subsets
    │
    ▼
7. Error-corrected SSU (MLE across replicates)
    └─ Inverse variance weighting:
        logit(SSU[v]) = Σ(logit(SSU[v|r]) / var(SSU[v|r])) / Σ(1 / var(SSU[v|r]))
        var_theta = 1 / Σ(1 / var(SSU[v|r]))
    │
    ▼
8. Empirical Bayes shrinkage
    └─ Shrink low-confidence variants toward global mean:
        theta_shrunk = λ * theta + (1 - λ) * mu_global
        λ = tau² / (tau² + var_theta)
    │
    ▼
9. Final SSU estimate:
    └─ ssu_corrected = plogis(theta_shrunk)
)")

# -- options -- #
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

# -- check options -- #
if(is.null(opt$rscript_dir)) stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$sample_id))   stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$ssu_counts))  stop("-d, list of splicing counts is required!", call. = FALSE)

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))

# -- inputs -- #
sample_reps      <- unlist(strsplit(opt$sample_id, ","))
files_ssu_counts <- unlist(strsplit(opt$ssu_counts, ","))

sample_reps      <- mixedsort(sample_reps)
files_ssu_counts <- sort_paths_by_filename(files_ssu_counts)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- paste0(opt$prefix, ".ssu_per_base")

# -- 1. reading files and formating -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "1. reading input files ...")

ssu_counts <- list()
for(i in seq_along(sample_reps))
{    
    ssu_counts[[sample_reps[i]]] <- as.data.table(vroom(files_ssu_counts[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
}
dt_ssu <- rbindlist(ssu_counts, idcol = "reps")

# -- 2. calculate SSU with eps -- #
# Note: For low-count variants, SSU is pulled toward 0.5. This is intentional Bayesian shrinkage.
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "2. calculate SSU with eps ...")
eps <- 0.5
dt_ssu[max_cov == 0, max_cov := 0.1]
dt_ssu[, ssu_eps := (base_cov + eps) / (max_cov + 2 * eps)]

# -- 3. calculate logit(ssu) and variance multiplicative -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "3. calculate logit(ssu) and variance multiplicative ...")
dt_ssu[, logit_ssu := qlogis(ssu_eps)]
dt_ssu[, var_mult := 1 /(max_cov * ssu_eps * (1 - ssu_eps))]

dt_ssu_wide <- dcast(dt_ssu, var_id + base_pos ~ reps, value.var = c("logit_ssu", "var_mult"))
setorder(dt_ssu_wide, var_id, base_pos)

# -- 4. build replicate subsets and likelihood matrices -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "4. replicate subsets and likelihood matrices ...")
rep_subsets <- unlist(lapply(2:length(sample_reps), function(k) combn(sample_reps, k, simplify = FALSE)), recursive = FALSE)

# build full matrices once: n_rows x n_reps
Y_full <- as.matrix(dt_ssu_wide[, paste0("logit_ssu_", sample_reps), with = FALSE])
V_full <- as.matrix(dt_ssu_wide[, paste0("var_mult_",  sample_reps), with = FALSE])
colnames(Y_full) <- sample_reps
colnames(V_full) <- sample_reps

OK_full <- is.finite(Y_full) & is.finite(V_full) 
Y_full[!OK_full] <- 0 
V_full[!OK_full] <- 0

n_rows <- nrow(Y_full)

rm(dt_ssu_wide)
gc()

# -- 5. set likelihood functions and parallel workers -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "5. set likelihood functions and parallel workers  ...")

subset_nll <- function(repset) {
    m   <- length(repset)
    Yc  <- Y_full[, repset,  drop = FALSE]
    Vc  <- V_full[, repset,  drop = FALSE]
    OKc <- OK_full[, repset, drop = FALSE] * 1  # Multiplying by 1 coerces it to numeric
    keep <- rowSums(OKc) >= 2

    function(a_log_sub) {
        va   <- rep(a_log_sub, each = n_rows)
        VV   <- Vc + exp(va)
        W    <- OKc / VV
        WY   <- W * Yc

        A    <- .rowSums(W,  n_rows, m)
        B    <- .rowSums(WY, n_rows, m)
        Csum <- .rowSums(WY * Yc, n_rows, m)
        Lsum <- .rowSums(OKc * log(VV), n_rows, m)

        Ak <- A[keep]
        Bk <- B[keep]

        sum(0.5 * (Lsum[keep] + Csum[keep] - Bk^2 / Ak))
    }
}

n_cores <- 4
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c("Y_full", "V_full", "OK_full", "n_rows", "rep_subsets", "subset_nll"), envir = environment())
clusterEvalQ(cl, { nll_funs <- lapply(rep_subsets, subset_nll); NULL })

chunk_id <- cut(seq_along(rep_subsets), n_cores, labels = FALSE)
chunks <- split(seq_along(rep_subsets), chunk_id)

joint_nll <- function(a_log) {
    partials <- clusterApply(
        cl, 
        chunks, 
        function(idxs, a_log) 
        { 
            sum(
                vapply(
                    idxs, 
                    function(i) nll_funs[[i]](a_log[rep_subsets[[i]]]), 
                    numeric(1)
                )
            ) 
        }, 
        a_log)

    sum(unlist(partials))
}


# nll_funs <- lapply(rep_subsets, subset_nll)

# joint_nll <- function(a_log) {
#     total <- 0
#     for (i in seq_along(rep_subsets)) {
#         total <- total + nll_funs[[i]](a_log[rep_subsets[[i]]])
#     }
#     total
# }

# -- 6. estimate error variances -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "6. estimate error variances (long running) ...")

init <- rep(log(0.01), length(sample_reps))
names(init) <- sample_reps

fit <- optim(par = init, fn = joint_nll, method = "L-BFGS-B")
stopCluster(cl)

var_addi_est <- exp(fit$par)

save.image("/lustre/scratch126/gengen/projects_v2/lehner_splicing/test/test_out/correct_ssu/test.img")

# -- 7. calculate error-corrected SSU -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "7. calculate error-corrected SSU ...")
dt_ssu[, var_addi := var_addi_est[reps]]
dt_ssu[, var_total := var_mult + var_addi]

dt_ssu_wide <- dcast(dt_ssu, var_id + base_pos ~ reps, value.var = c("base_ssu", "max_cov"))
base_ssu_cols <- grep("^base_ssu_", names(dt_ssu_wide), value = TRUE)
max_cov_cols <- grep("^max_cov_", names(dt_ssu_wide), value = TRUE)
setnames(dt_ssu_wide, base_ssu_cols, paste0("ssu", seq_along(base_ssu_cols)))
setnames(dt_ssu_wide, max_cov_cols, paste0("mcov", seq_along(max_cov_cols)))

# inverse variance weighting
dt_ssu_corrected <- dt_ssu[, .(theta = sum(logit_ssu / var_total) / sum(1 / var_total),
                               var_theta = 1 / sum(1 / var_total)),
                               by = .(var_id, base_pos)]
dt_ssu_corrected[, ssu_est := plogis(theta)]

dt_ssu_corrected <- merge(dt_ssu_wide, dt_ssu_corrected, by = c("var_id", "base_pos"), all.x = TRUE)

# -- 8. shrinkage (empirical Bayes) -- #
# We now shrink variant estimates toward a global mean, exactly as DiMSum does for fitness.
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "8. shrinkage (empirical Bayes) ...")

mu_global <- mean(dt_ssu_corrected$theta, na.rm = TRUE)
tau2 <- var(dt_ssu_corrected$theta, na.rm = TRUE)

# shrinkage factor
# for each variant per pos: λ(v,p) = tau2 / (tau2 + var_theta(v,p))
dt_ssu_corrected[, shrinkage := tau2 / (tau2 + var_theta)]

# shrunk estimates
# θ(shrunk)​ = λ(v,p)​θ(v,p)​ + (1−λ(v,p)​) * mu_global
dt_ssu_corrected[, theta_shrunk := shrinkage * theta + (1 - shrinkage) * mu_global]
dt_ssu_corrected[, var_theta_shrunk := shrinkage^2 * var_theta]
dt_ssu_corrected[, ssu_corrected := plogis(theta_shrunk)]

# -- 9. confidence intervals -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "9. calculate confidence intervals ...")

# θ(hat​) = logit(SSU) ∼ N(θ,Var(θ(est)​))
# a standard normal variable Z ~ N(0,1), so P(|Z| <= 1.96) = 0.95
# then 95% CI = mean ± 1.96 × SD
z <- 1.96
dt_ssu_corrected[, ssu_corrected_lwr := plogis(theta_shrunk - z * sqrt(var_theta_shrunk))]
dt_ssu_corrected[, ssu_corrected_upr := plogis(theta_shrunk + z * sqrt(var_theta_shrunk))]

# -- 10. output -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "10. output ...")
num_cols <- names(dt_ssu_corrected)[sapply(dt_ssu_corrected, is.numeric)]
dt_ssu_corrected[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
output_file <- file.path(opt$output_dir, paste0(sample_prefix, ".details.tsv"))
fwrite(dt_ssu_corrected, file = output_file, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE)
