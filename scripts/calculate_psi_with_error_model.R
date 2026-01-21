#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "glue", "tidyverse", "data.table", "vroom", "gtools")
invisible(lapply(packages, quiet_library))

# -- Model description: PSI error model -- #
# -------------------------
# - PSI = Percent Spliced In
# - Observed counts: Inclusion, Skipping
# - Add small prior eps = 0.5 (Beta(0.5, 0.5), Jeffreys prior)
#   -> stabilizes low-count variants, avoids logit(0/1)
# - Variance model: multiplicative (Poisson/binomial) + additive (replicate noise)
# - Inverse variance weighting used to combine replicates
# - Empirical Bayes shrinkage reduces noise in low-confidence variants
# -------------------------
model_help_description <- glue(r"(
Calculate PSI with error model from splicing counts.

Reason: 
Read counts are finite and some variants may have low coverage, and replicates may differ due to technical noise.

Methods:
Target is PSI (Percent Spliced In) defined as the proportion of transcripts that include a given exon or splice junction.
And PSI should be in [0,1]

Asumptions: Binomial observation (matches the sequencing process)
Given N reads for a variant
Inclusion ~ Binomial(N, PSI)
Skipping  ~ Binomial(N, 1 - PSI)

1. Adding error term (eps) to avoid undefined PSI values
Raw PSI = Inclusion / (Inclusion + Skipping)
Breaks when Inclusion = 0 or Skipping = 0
To avoid this, we add a small error term (eps) to both Inclusion and Skipping
Adjusted PSI = (Inclusion + eps) / (Inclusion + Skipping + 2 * eps)
We use eps = 0.5, which corresponds to a Bayesian posterior mean of PSI under a Beta prior (Beta(0.5, 0.5), Jeffreys prior for binomial proportion)

2. Logit transformation
PSI in [0,1], may casue issue:
- Variance dpends on PSI
- Near 0 or 1 causes asymmetric noise
- 0 count causes undefined values
To stabilize variance and make the distribution more normal-like, we apply logit transformation to PSI
logit(PSI) = log(PSI / (1 - PSI)) = log((Inclusion + eps) / (Skipping + eps)
- Maps PSI from [0,1] to (-inf, +inf)
- Make noise appproximately Gaussaian for counts
- Allows linear modeling

3. Variance of logit(PSI): Poisson/Binomial noise model (Sampling)
Var(Inclusion) = N * PSI * (1 - PSI)
Var(PSI) = N * PSI * (1 - PSI) / N² = PSI * (1 - PSI) / N
Var(logit(PSI)) = 1 / (N * PSI * (1 - PSI))

4. Error model overview
For a variant v in replicate r
logit(PSI_vr) = logit(PSI_v0) + error(vr)[M] + error(vr)[A]
Where:
- logit(PSI_v0) is the true underlying logit(PSI) for variant v
- error(vr)[M] is a multiplicative error term
    - from sampling noise (finite read counts): Var(error(vr)[M]) = Var(logit(PSI_vr)) = 1 /(N_vr * PSI_vr * (1 - PSI_vr))
    - when N_vr is small, variance is large versus when N_vr is large, variance is small
    - error(vr)[M] ~ N(0, Var(logit(PSI_vr))) = N(0, 1 /(N_vr * PSI_vr * (1 - PSI_vr)))
- error(vr)[A] is an additive error term
    - from replicate (library prep, sequencing, PCR): Var(error(vr)[A]) = sigma_rep²
    - error(vr)[A] ~ N(0, Var(rep)) = N(0, sigma_rep²)
Total variance of error
    - Var(logit(PSI_v0)) = 1 /(N_vr * PSI_vr * (1 - PSI_vr)) + sigma_rep²

5. DiMSum approach - subset replicates
Issue: 
observed variance of variant v: Var(PSI_vr) = a1² + a2² + a3²
This is one equation with three unknowns. The likelihood is invariant under: (a1² + δ, a2² - δ, a3²)
Additive errors are not identifiable

How DiMSum solves this:
Fit the error model simultaneously on: R ∈ ((1,2,3),(1,2),(1,3),(2,3))
(1,2,3) -> a1² + a2² + a3²
(1,2)   -> a1² + a2²
(1,3)   -> a1² + a3²
(2,3)   -> a2² + a3²
This gives 4 equations with 3 unknowns, allowing estimation of all parameters

6. Estimation of replicate variance
- For each variant, compute variance across replicates
- Subtract expected sampling variance (from read counts)
- Pool across variants to get a robust estimate of replicate variance (robust median or regression vs depth)

7. Combining replicates: inverse variance weighting
Best estimate logit(PSI) for variant v:
logit(PSI_v) = sum_r (Y(v|r) / var_total(v|r)) / sum_r (1 / var_total(v|r))
where var_total(v|r) = Var(logit(PSI_vr)) + Var(rep)
- Maximal likelihood estimate under Gaussian noise
- Replicates with lower variance (higher coverage) get higher weight
- Low-coverage replicates contribute less

8. Shrinkage
Low-coverage varaints:
- High variance
- Small wieght
- Strong shrinkage to population mean
High-coverage variants:
- Low variance
- Large weight
- Little shrinkage to population mean
Overall, this model accounts for technical noise from finite read counts and replicate variability,
providing robust PSI estimates across a range of coverage levels.
)")

model_workflow <- glue(r"(
Raw counts per variant
 ├─ Inclusion counts
 └─ Skipping counts
        │
        ▼
Add Bayesian prior (eps = 0.5)
 └─ Adjusted PSI:
      PSI_adj = (Inclusion + eps) / (Inclusion + Skipping + 2*eps)
        │
        ▼
Logit transformation
 └─ logit(PSI) = log(PSI_adj / (1 - PSI_adj))
        │
        ▼
Compute multiplicative variance (sampling noise)
 └─ var_mult = 1 / (n_total * PSI_adj * (1 - PSI_adj))
        │
        ▼
Replicate subsets (DiMSum approach)
 └─ Build all combinations of replicates
        │
        ▼
Estimate additive replicate variance (σ²_rep)
 └─ Optimize Negative Log-Likelihood across subsets
        │
        ▼
Error-corrected PSI (MLE across replicates)
 └─ Inverse variance weighting:
      logit(PSI_v) = Σ(logit_psi / var_tot) / Σ(1 / var_tot)
      var_theta = 1 / Σ(1 / var_tot)
        │
        ▼
Empirical Bayes shrinkage
 └─ Shrink low-confidence variants toward global mean:
      theta_shrunk = λ * theta + (1 - λ) * mu_global
      λ = tau² / (tau² + var_theta)
        │
        ▼
Final PSI estimate:
 └─ psi_shrunk = plogis(theta_shrunk)
)")

# -- options -- #
option_list <- list(make_option(c("-r", "--rscript_dir"),     type = "character",    help = "directory path of R scripts", default = NULL),
                    make_option(c("-s", "--sample_id"),       type = "character",    help = "list of sample IDs",          default = NULL),
                    make_option(c("-d", "--splicing_counts"), type = "character",    help = "list of splicing counts",     default = NULL),
                    make_option(c("-o", "--output_dir"),      type = "character",    help = "output directory",            default = getwd()),
                    make_option(c("-p", "--prefix"),          type = "character",    help = "output prefix",               default = "sample"),
                    make_option(c("-m", "--model_help"),      action = "store_true", help = "print model description",     default=FALSE),
                    make_option(c("-w", "--model_workflow"),  action = "store_true", help = "print model workflow",        default=FALSE))

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

if(opt$model_workflow)
{
    print(model_workflow)
    quit(status = 0)
}

# -- check options -- #
if(is.null(opt$rscript_dir))          stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$sample_id))            stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$splicing_counts))      stop("-d, list of splicing counts is required!", call. = FALSE)

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))

# -- inputs -- #
sample_reps           <- unlist(strsplit(opt$sample_id, ","))
files_splicing_counts <- unlist(strsplit(opt$splicing_counts, ","))

sample_reps           <- mixedsort(sample_reps)
files_splicing_counts <- sort_paths_by_filename(files_splicing_counts)

# -- outputs -- #
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
setwd(opt$output_dir)

sample_prefix <- opt$prefix

# -- 1. reading files and formating -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "1. reading input files ...")
reformat_counts <- function(dt, inclusion_cols, skipping_cols)
{
    missing <- setdiff(c("var_id", inclusion_cols, skipping_cols), names(dt))
    if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

    other_cols <- setdiff(names(dt), c("var_id", inclusion_cols, skipping_cols))

    dt <- dt[, .( var_id    = var_id,
                  inclusion = rowSums(.SD[, inclusion_cols, with = FALSE]),
                  skipping  = rowSums(.SD[, skipping_cols,  with = FALSE]),
                  others    = rowSums(.SD[, other_cols,     with = FALSE]))]

    dt[, ratio_inclusion := inclusion / (inclusion + skipping + others)]
    dt[, ratio_skipping  := skipping  / (inclusion + skipping + others)]
    dt[, ratio_others    := others    / (inclusion + skipping + others)]
    dt[, ratio_splicing := { incl <- as.integer(ratio_inclusion * 100)
                             skip <- as.integer(ratio_skipping * 100)
                             other <- 100 - incl - skip
                             paste0(incl, ":", skip, ":", other) }]
    return(dt)
}

cols_canonical_inclusion <- c("canonical_inclusion")
cols_canonical_skipping  <- c("canonical_skipping")
splicing_counts <- list()
for(i in seq_along(sample_reps))
{    
    splicing_counts[[sample_reps[i]]] <- as.data.table(vroom(files_splicing_counts[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    splicing_counts[[sample_reps[i]]] <- reformat_counts(splicing_counts[[sample_reps[i]]], cols_canonical_inclusion, cols_canonical_skipping)
}
dt_splicing <- rbindlist(splicing_counts, idcol = "reps")


# -- 2. calculate PSI with eps -- #
# Note: For low-count variants, PSI is pulled toward 0.5. This is intentional Bayesian shrinkage.
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "2. calculate PSI with eps ...")
eps <- 0.5
dt_splicing[, n_total := inclusion + skipping]
dt_splicing <- dt_splicing[n_total > 0]
dt_splicing[, psi := (inclusion + eps) / (n_total + 2*eps)]

# -- 3. calculate logit(psi) and variance multiplicative -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "3. calculate logit(psi) and variance multiplicative ...")
dt_splicing[, logit_psi := qlogis(psi)] # qlogis(psi) = log(psi / (1 - psi))
dt_splicing[, var_mult := 1 /(n_total * psi * (1 - psi))]

# -- 4. build replicate subsets -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "4. build replicate subsets ...")
rep_subsets <- unlist(lapply(2:length(sample_reps), function(k) combn(sample_reps, k, simplify = FALSE)), recursive = FALSE)

# Negative log-likelihood (NLL)
# We usually minimize functions in optimization routines like optim() in R.
# So instead of maximizing L(θ), we minimize: NLL(θ) = -log(L(θ))
subset_nll <- function(a_log, dt_sub)
{
    # convert to actual additive variances
    a2 <- exp(a_log)
    dt_sub[, var_tot := var_mult + a2[reps]]

    # theta_i estimate
    theta_dt <- dt_sub[, .(theta = sum(logit_psi / var_tot) / sum(1 / var_tot)), by = var_id]
    dt_sub <- merge(dt_sub, theta_dt, by = "var_id")

    return(sum(0.5 * (log(dt_sub$var_tot) + (dt_sub$logit_psi - dt_sub$theta)^2 / dt_sub$var_tot)))
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

# initial values for log(a_r^2)
init <- rep(log(0.01), length(sample_reps))
names(init) <- sample_reps

# estimates
fit <- optim(par = init, fn = joint_nll, dt = dt_splicing, rep_subsets = rep_subsets, method = "BFGS")
a2_hat <- exp(fit$par)
a_hat  <- sqrt(a2_hat)

# -- 6. calculate error-corrected PSI -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "6. calculate error-corrected PSI ...")
dt_splicing[, a2 := a2_hat[reps]]
dt_splicing[, var_tot := var_mult + a2]

# inverse variance weighting
# PSI_v0 = logit(PSI_v)
# PSI_v0 = sum_r (Y(v|r) / var_total(v|r)) / sum_r (1 / var_total(v|r))
# Var(PSI_v0) = 1 / sum_r (1 / var_total(v|r))
ratio_wide <- dcast(dt_splicing, var_id ~ reps, value.var = "ratio_splicing")
setnames(ratio_wide, old = names(ratio_wide)[-1], new = paste0("ratio", seq_along(names(ratio_wide)[-1])))

dt_splicing_corrected <- dt_splicing[, .(theta = sum(logit_psi / var_tot) / sum(1 / var_tot),
                                         var_theta = 1 / sum(1 / var_tot)), 
                                         by = var_id]
dt_splicing_corrected[, psi_est := plogis(theta)]

dt_splicing_corrected <- merge(ratio_wide, dt_splicing_corrected, by = "var_id", all.x = TRUE)

# -- 7. shrinkage (empirical Bayes) -- #
# We now shrink variant estimates toward a global mean, exactly as DiMSum does for fitness.
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "7. shrinkage (empirical Bayes) ...")

mu_global <- mean(dt_splicing_corrected$theta)
tau2 <- var(dt_splicing_corrected$theta)

# shrinkage factor
# for each variant v: λ(v) = tau2 / (tau2 + var_theta(v))
dt_splicing_corrected[, shrinkage := tau2 / (tau2 + var_theta)]

# shrunk estimates
# θ(shrunk)​ = λ(v)​θ(v)​ + (1−λ(v)​) * mu_global
dt_splicing_corrected[, theta_shrunk := shrinkage * theta + (1 - shrinkage) * mu_global]
dt_splicing_corrected[, var_theta_shrunk := shrinkage^2 * var_theta]
dt_splicing_corrected[, psi_shrunk := plogis(theta_shrunk)]

# -- 8. confidence intervals -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "8. calculate confidence intervals ...")

# θ(hat​) = logit(PSI) ∼ N(θ,Var(θ(hat)​))
# a standard normal variable Z ~ N(0,1), so P(|Z| <= 1.96) = 0.95
# then 95% CI = mean ± 1.96 × SD
z <- 1.96
dt_splicing_corrected[, psi_shrunk_lwr   := plogis(theta_shrunk - z * sqrt(var_theta_shrunk))]
dt_splicing_corrected[, psi_shrunk_upr   := plogis(theta_shrunk + z * sqrt(var_theta_shrunk))]

# -- 9. output -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "9. output ...")
num_cols <- names(dt_splicing_corrected)[sapply(dt_splicing_corrected, is.numeric)]
dt_splicing_corrected[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
output_file <- file.path(opt$output_dir, paste0(sample_prefix, ".corrected_psi.tsv"))
fwrite(dt_splicing_corrected, file = output_file, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE)
