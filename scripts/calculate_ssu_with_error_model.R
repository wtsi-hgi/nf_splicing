#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "glue", "tidyverse", "data.table", "vroom", "gtools")
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

# -- 4. build replicate subsets -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "4. build replicate subsets ...")
rep_subsets <- unlist(lapply(2:length(sample_reps), function(k) combn(sample_reps, k, simplify = FALSE)), recursive = FALSE)

subset_nll <- function(a_log, dt_sub)
{
    var_addi <- exp(a_log)
    dt_sub[, var_total := var_mult + var_addi[reps]]

    dt_theta <- dt_sub[, .(theta = sum(logit_ssu / var_total) / sum(1 / var_total)), by = var_id]
    dt_sub <- merge(dt_sub, dt_theta, by = "var_id")

    return(sum(0.5 * (log(dt_sub$var_total) + (dt_sub$logit_ssu - dt_sub$theta)^2 / dt_sub$var_total)))
}

joint_nll <- function(a_log, dt, rep_subsets)
{
    total_nll <- 0
    
    for(repset in rep_subsets) 
    {
        dt_sub <- dt[reps %in% repset]
        a_log_sub <- a_log[repset]
        total_nll <- total_nll + subset_nll(a_log_sub, dt_sub)
    }

    return(total_nll)
}

# -- 5. estimate error variances -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "5. estimate error variances (long running) ...")

# initial values for log(a{rep}^2)
init <- rep(log(0.01), length(sample_reps))
names(init) <- sample_reps

# estimates of additive variances for each replicate
fit <- optim(par = init, fn = joint_nll, dt = dt_ssu, rep_subsets = rep_subsets, method = "BFGS")
var_addi_est <- exp(fit$par)

# -- 6. calculate error-corrected SSU -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "6. calculate error-corrected SSU ...")
dt_ssu[, var_addi := var_addi_est[reps]]
dt_ssu[, var_total := var_mult + var_addi]

dt_wide <- dcast(dt_ssu, var_id ~ reps, value.var = c("base_ssu", "max_cov"))
base_ssu_cols <- grep("^base_ssu_", names(dt_wide), value = TRUE)
max_cov_cols <- grep("^max_cov_", names(dt_wide), value = TRUE)
setnames(dt_wide, base_ssu_cols, paste0("ssu", seq_along(base_ssu_cols)))
setnames(dt_wide, max_cov_cols, paste0("max_cov", seq_along(n_total_cols)))

# inverse variance weighting
dt_ssu_corrected <- dt_ssu[, .(theta = sum(logit_ssu / var_total) / sum(1 / var_total),
                                         var_theta = 1 / sum(1 / var_total)), 
                                         by = var_id]
dt_ssu_corrected[, ssu_est := plogis(theta)]

dt_ssu_corrected <- merge(dt_wide, dt_ssu_corrected, by = "var_id", all.x = TRUE)

# -- 7. shrinkage (empirical Bayes) -- #
# We now shrink variant estimates toward a global mean, exactly as DiMSum does for fitness.
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "7. shrinkage (empirical Bayes) ...")

mu_global <- mean(dt_ssu_corrected$theta)
tau2 <- var(dt_ssu_corrected$theta)

# shrinkage factor
# for each variant v: λ(v) = tau2 / (tau2 + var_theta(v))
dt_ssu_corrected[, shrinkage := tau2 / (tau2 + var_theta)]

# shrunk estimates
# θ(shrunk)​ = λ(v)​θ(v)​ + (1−λ(v)​) * mu_global
dt_ssu_corrected[, theta_shrunk := shrinkage * theta + (1 - shrinkage) * mu_global]
dt_ssu_corrected[, var_theta_shrunk := shrinkage^2 * var_theta]
dt_ssu_corrected[, ssu_corrected := plogis(theta_shrunk)]

# -- 8. confidence intervals -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "8. calculate confidence intervals ...")

# θ(hat​) = logit(SSU) ∼ N(θ,Var(θ(est)​))
# a standard normal variable Z ~ N(0,1), so P(|Z| <= 1.96) = 0.95
# then 95% CI = mean ± 1.96 × SD
z <- 1.96
dt_ssu_corrected[, ssu_corrected_lwr   := plogis(theta_shrunk - z * sqrt(var_theta_shrunk))]
dt_ssu_corrected[, ssu_corrected_upr   := plogis(theta_shrunk + z * sqrt(var_theta_shrunk))]

# -- 9. output -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "9. output ...")
num_cols <- names(dt_ssu_corrected)[sapply(dt_ssu_corrected, is.numeric)]
dt_ssu_corrected[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
output_file <- file.path(opt$output_dir, paste0(sample_prefix, ".details.tsv"))
fwrite(dt_ssu_corrected, file = output_file, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE)
