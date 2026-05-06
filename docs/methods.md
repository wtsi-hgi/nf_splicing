## PSI calculation using an error model

### Rationale

Splicing measurements derived from RNA sequencing are affected by technical uncertainty due to finite read depth, low-coverage variants, and variability between replicate libraries. To obtain robust estimates of percent spliced in (PSI), we model both sampling noise and replicate-specific variability within a unified error framework.

---

### Observation model

For each variant (exon or splice junction), $\mathrm{PSI}$ is defined as the proportion of transcripts that include the variant and is constrained to the interval $[0,1]$.

Sequencing reads are assumed to arise from a binomial sampling process. For a variant with total read depth $\mathrm{N}$:

$$
\begin{aligned}
\mathrm{I} &\sim \mathrm{Binomial}(\mathrm{N}, \mathrm{PSI}) \\
\mathrm{S} &\sim \mathrm{Binomial}(\mathrm{N}, 1 - \mathrm{PSI})
\end{aligned}
$$

where $\mathrm{I}$ and $\mathrm{S}$ denote inclusion and skipping reads, respectively.

---

### Pseudocount adjustment

The raw PSI estimator,

$$
\mathrm{PSI} = \frac{\mathrm{I}}{\mathrm{I} + \mathrm{S}}
$$

is undefined when $\mathrm{I} = 0$ or $\mathrm{S} = 0$. To avoid this issue, we add a small pseudocount $\varepsilon$ to both inclusion and skipping reads:

$$
\mathrm{PSI}_{\text{adj}} =
\frac{\mathrm{I} + \varepsilon}{\mathrm{I} + \mathrm{S} + 2\varepsilon}
$$

We set $\varepsilon = 0.5$, corresponding to the posterior mean under a Jeffreys prior $\mathrm{Beta}(0.5, 0.5)$ for a binomial proportion. This provides a weakly informative Bayesian regularization and stabilizes estimates at low coverage.

---

### Logit transformation

Because PSI is bounded between 0 and 1, its variance depends on its mean and becomes highly asymmetric near the boundaries. To stabilize variance and enable Gaussian modeling, we apply a logit transformation:

$$
\theta = \mathrm{logit}(\mathrm{PSI}_{\text{adj}})
= \log\left(\frac{\mathrm{PSI}_{\text{adj}}}{1 - \mathrm{PSI}_{\text{adj}}}\right)
= \log\left(\frac{\mathrm{I} + \varepsilon}{\mathrm{S} + \varepsilon}\right)
$$

The logit transformation:

* Maps $\mathrm{PSI}$ from $[0,1]$ to $(-\infty, +\infty)$, stabilizes variance
* Reduces heteroscedasticity
* Makes sampling noise approximately Gaussian, allows Gaussian error modeling
* Enables linear and inverse-variance weighted modeling

---

### Sampling variance of logit-PSI

Under the binomial observation model, the number of reads supporting inclusion of a given exon follows a Binomial distribution:

$$
\mathrm{I} \sim \mathrm{Binomial}(\mathrm{N}, \mathrm{PSI})
$$

**Variance of inclusion counts:**

The variance of the number of inclusion reads $I$ is given by the binomial variance:

$$
\mathrm{Var(I)} = \mathrm{N} \cdot \mathrm{PSI} \cdot (1 - \mathrm{PSI})
$$

**Variance of PSI:**

The observed PSI is 

$$
\mathrm{PSI} = \frac{\mathrm{I}}{\mathrm{N}}
$$

Using standard variance scaling:

$$
{\mathrm{Var(PSI)}} = \frac{\mathrm{Var(I)}}{\mathrm{N^2}} = \frac{\mathrm{N} \cdot \mathrm{PSI} \cdot (1 - \mathrm{PSI})}{\mathrm{N}^2} = \frac{\mathrm{PSI} \cdot (1 - \mathrm{PSI})}{\mathrm{N}}
$$

**Variance of logit(PSI) via delta method:**

Let:

$$
\theta = \mathrm{logit(PSI)} = \log\frac{\mathrm{PSI}}{1-\mathrm{PSI}}\
$$ 

Using the delta method:

$$
\mathrm{Var}(\theta)
\approx
\left(\frac{d\theta}{d\mathrm{PSI}}\right)^2 \mathrm{Var}(\hat{\mathrm{PSI}})
$$

$$
\mathrm{Var}(\theta)
\approx
(\frac{1}{\mathrm{PSI} \cdot \mathrm{(1-PSI)}})^2 \cdot \frac{\mathrm{PSI} \cdot \mathrm{(1-PSI)}}{\mathrm{N}}
$$

$$
\mathrm{Var}(\theta)
\approx
\frac{1}{\mathrm{N} \cdot \mathrm{PSI} \cdot (1 - \mathrm{PSI})}
$$

This term represents sampling uncertainty due to finite read depth. It captures multiplicative (depth-dependent) uncertainty due to finite read counts. Variance decreases with increasing coverage and is largest when PSI is near 0 or 1.

---

### Error model

For variant $v$ in replicate $r$, the observed logit-PSI is modeled as:

$$
\theta_{vr} = \theta_{v} + \varepsilon^{(M)}_{vr} + \varepsilon^{(A)}_{r}
$$

where:
* $\theta_{v}$ is the true underlying logit-PSI of variant ${v}$
* $\varepsilon^{(M)}_{vr}$ is a multiplicative error term arising from sampling noise
    * variants with small $\mathrm{N}_{vr}$ show higher sampling variance, while variants with large $\mathrm{N}_{vr}$ have reduced variance.
* $\varepsilon^{(A)}_{r}$ is an additive replicate-specific error term.
    * it captures technical variation from library preparation, PCR amplification, and sequencing.

**Multiplicative error (sampling noise):**

$$
\varepsilon^{(M)}_{vr} \sim
\mathcal{N}\left(
0,\;
\frac{1}{\mathrm{N}_{vr} \cdot \mathrm{PSI}_{vr} \cdot (1-\mathrm{PSI}_{vr})}
\right)
$$

**Additive error (replicate variability):**

$$
\varepsilon^{(A)}_{r} \sim \mathcal{N}(0,\; \sigma_r^2)
$$

**The total variance of $\theta_{vr}$ is therefore:**

$$
\mathrm{Var}(\theta_{vr}) =
\frac{1}{\mathrm{N}_{vr} \cdot \mathrm{PSI}_{vr} \cdot (1-\mathrm{PSI}_{vr})} + \sigma_r^2
$$

---

### Identifiability via subset replicates

**Issue:**

For a given variant $v$, the observed variance across three replicates can be written as

$$
\mathrm{Var(PSI_{vr})} = \sigma_{1}^2 + \sigma_{2}^2 + \sigma_{3}^2
$$

where $\mathrm{\sigma_{r}^2}$ denotes the replicate-specific error variance.

This provides a single equation with three unknowns, making the model non-identifiable: the likelihood is invariant under transformations such as

$$
(\sigma_{1}^2 + \delta, \sigma_{2}^2 - \delta, \sigma_{3}^2)
$$

and therefore additive error components cannot be uniquely estimated.

When only full replicate sets are used, additive replicate variances are not identifiable because the observed variance is invariant under redistributions among replicates. To resolve this, we adopt the DiMSum strategy of jointly fitting the model across all subsets of replicates.

For three replicates, likelihoods are evaluated on:

$$
(1,2,3): \sigma_{1}^2 + \sigma_{2}^2 + \sigma_{3}^2
$$

$$
(1,2): \sigma_{1}^2 + \sigma_{2}^2
$$

$$
(2,3): \sigma_{2}^2 + \sigma_{3}^2
$$

$$
(1,3): \sigma_{1}^2 + \sigma_{3}^2
$$

providing multiple equations that allow unique estimation of replicate-specific additive variances $\sigma_r^2$

---

### Estimation of replicate variance

For each variant, the observed variance of PSI across replicates is decomposed into sampling variance and replicate-specific additive variance. The sampling variance is determined by read depth and the estimated PSI and is treated as known.

Replicate-specific additive variances $\sigma_r^2$ are assumed to be shared across variants and are estimated by pooling information across all variants. To ensure identifiability, the error model is fit jointly across all replicate subsets, and parameters are learned by maximizing the joint likelihood.

This pooling strategy provides robust estimates of replicate variance and reduces the influence of low-coverage or noisy variants.

To ensure identifiability and robust estimation of $\sigma_r^2$ , parameters are estimated jointly across all replicate subsets. For each subset $S$, the negative log-likelihood is

$$
L(S) = \sum_{r \in S} \sum_{v} \frac{1}{2} \left[ \log(\sigma_{vr}^2) + \frac{(\theta_{vr} - \theta_v)^2}{\sigma_{vr}^2} \right]
$$

$\theta_v$ is computed using only observations within subset S

The total objective function is defined as the sum over all subsets:

$$
L(total) = \sum_{S} L(S) 
$$

Then following the DiMSum framework, a joint negative log-likelihood was minimised over all replicate subsets simultaneously using the BFGS algorithm

$$
arg min_{\sigma_r^2} \sum_{S} \sum_{r \in S} \sum_v \frac{1}{2} \left[ \log(\sigma_{vr}^2) + \frac{(\theta_{vr} - \theta_v)^2}{\sigma_{vr}^2} \right]
$$

Then under the Gaussian error model

$$
\theta_{vr} \sim N(\theta_v, \sigma_{vr}^2)
$$

For a given set of variances, the maximum likelihood estimate of $\theta_v$ is given by the inverse-variance weighted mean.

$$
\theta_v = \frac{ \sum_r \frac{\theta_{vr}}{\sigma_{vr}^2}}{ \sum_r \frac{1}{\sigma_{vr}^2}}
$$


---

### Combining replicates by inverse-variance weighting

The maximum-likelihood estimate of the true logit-PSI for variant $v$ is computed by inverse-variance weighting:

$$
\hat{\theta}_v =
\frac{\sum_r \theta_{vr} / \mathrm{Var}(\theta_{vr})}
     {\sum_r 1 / \mathrm{Var}(\theta_{vr})}
$$

$$
\mathrm{Var}(\hat{\theta}_v) =
\frac{1}{\sum_r 1 / \mathrm{Var}(\theta_{vr})}
$$

Replicates with higher coverage and lower technical variance contribute more weight, while low-coverage replicates contribute less.

---

### Empirical Bayes shrinkage

To stabilize estimates for low-coverage variants, logit-PSI estimates are shrunk toward the global mean using an empirical Bayes framework. Variants with large uncertainty exhibit stronger shrinkage, whereas high-coverage variants remain largely unchanged. Final PSI estimates are obtained by inverse logit transformation.

Assuming a normal prior
$$
\theta_v \sim N(\mu, \tau^2)
$$

$\mu$ and $\tau^2$ are the mean and variance of logit-PSI across all variants, the shrinkage factor for each variant is:

$$
\lambda_v = \frac{\tau^2}{\tau^2 + \sigma_v^2}
$$

The posterior (shrunk) estimate and variance are given by:

$$
\bar{\theta_v} = \lambda_v \cdot \theta_v + (1 - \lambda_v) \cdot \mu
$$

$$
\bar{\sigma_v} = \lambda_v^2 \cdot \sigma_v^2
$$

The final PSI estimate and 95% confidence interval were obtained by back-transforming through the logistic function. Only variants with at least one replicate meeting a minimum coverage of 10 reads were retained for downstream analysis.

---

### Summary

This error model explicitly accounts for both sampling noise from finite read counts and replicate-specific technical variability. By combining logit transformation, variance modeling, subset-based variance estimation, and empirical Bayes shrinkage, the approach yields robust PSI estimates across a wide range of sequencing depths.

---