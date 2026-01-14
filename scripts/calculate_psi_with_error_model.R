#!/usr/bin/env Rscript
quiet_library <- function(pkg) { suppressMessages(suppressWarnings(library(pkg, character.only = TRUE))) }
packages <- c("optparse", "glue", "tidyverse", "data.table", "vroom", "gtools")
invisible(lapply(packages, quiet_library))

# -- description -- #
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

# -- options -- #
option_list <- list(make_option(c("-r", "--rscript_dir"),     type = "character",    help = "directory path of R scripts", default = NULL),
                    make_option(c("-s", "--sample_id"),       type = "character",    help = "list of sample IDs",          default = NULL),
                    make_option(c("-d", "--splicing_counts"), type = "character",    help = "list of splicing counts",     default = NULL),
                    make_option(c("-o", "--output_dir"),      type = "character",    help = "output directory",            default = getwd()),
                    make_option(c("-p", "--prefix"),          type = "character",    help = "output prefix",               default = "sample"),
                    make_option(c("-m", "--model_help"),      action = "store_true", help = "print model description"))

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
if(is.null(opt$rscript_dir))          stop("-r, directory path of R scripts is required!", call. = FALSE)
if(is.null(opt$sample_id))            stop("-s, list of sample IDs is required!", call. = FALSE)
if(is.null(opt$splicing_counts))      stop("-d, list of splicing counts is required!", call. = FALSE)

# -- modules -- #
source(file.path(opt$rscript_dir, "report_utils.R"))

# -- inputs -- #
sample_reps                <- unlist(strsplit(opt$sample_id, ","))
files_splicing_counts      <- unlist(strsplit(opt$splicing_counts, ","))

sample_reps                <- mixedsort(sample_reps)
files_splicing_counts      <- sort_paths_by_filename(files_splicing_counts)

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

    return(dt[, .( var_id    = var_id,
                   inclusion = rowSums(.SD[, ..inclusion_cols]),
                   skipping  = rowSums(.SD[, ..skipping_cols]))])
}

for(i in seq_along(sample_reps))
{    
    splicing_counts[[sample_reps[i]]] <- as.data.table(vroom(files_splicing_counts[i], delim = "\t", comment = "#", col_names = TRUE, show_col_types = FALSE))
    splicing_counts[[sample_reps[i]]] <- reformat_counts(splicing_counts[[sample_reps[i]]], "canonical_inclusion", "canonical_skipping")
}
dt_splicing <- rbindlist(splicing_counts, idcol = "reps")

# -- 2. calculate PSI with eps -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "2. calculate PSI with eps ...")
eps <- 0.5
dt_splicing[, n_total := inclusion + skipping]
dt_splicing[, psi := (inclusion + eps) / (n_total + 2*eps)]

# -- 3. calculate logit(psi) and variance multiplicative -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "3. calculate logit(psi) and variance multiplicative ...")
# qlogis(psi) = log(psi / (1 - psi))
dt_splicing[, logit_psi := qlogis(psi)]
dt_splicing[, var_mult := 1 /(n_total * psi * (1 - psi))]

# -- 4. build replicate subsets -- #
message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "4. build replicate subsets ...")
rep_subsets <- unlist(lapply(2:length(sample_reps), function(k) combn(sample_reps, k, simplify = FALSE)), recursive = FALSE)

# Negative log-likelihood (NLL)
# We usually minimize functions in optimization routines like optim() in R.
# So instead of maximizing L(θ), we minimize: NLL(θ) = -log(L(θ))
subset_variance <- function(a_log, dt_sub)
{
    # convert to actual additive variances
    a2 <- exp(a_log)
    names(a2) <- unique(dt_sub$reps)
    
    dt_sub[, var_tot := var_mult + a2[reps]]

    # theta_i estimate
    theta_dt <- dt_sub[, .(theta = sum(logit_psi / var_tot) / sum(1 / var_tot)), by = var_id]
    
    dt_sub <- merge(dt_sub, theta_dt, by = "var_id")
    sum(0.5 * (log(var_tot) + (logit_psi - theta)^2 / var_tot))
}
