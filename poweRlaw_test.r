library("poweRlaw")
library("dplyr")
# https://cran.r-project.org/web/packages/poweRlaw/vignettes/c_comparing_distributions.pdf

# Safe wrapper for compare_distributions that handles NA values
safe_compare <- function(dist1, dist2) {
  result <- tryCatch({
    compare_distributions(dist1, dist2)
  }, error = function(e) {
    list(test_statistic = NA, p_two_sided = NA)
  })

  # If the result contains NA, return a list with NA values
  if (is.na(result$test_statistic) || is.na(result$p_two_sided)) {
    return(list(test_statistic = NA, p_two_sided = NA))
  }
  return(result)
}

# Test how many samples you need to distinguish between distributions
# For three true distributions: power law, binomial, and lognormal
# Run 10 replicates for each combination

sample_sizes <- c(50, 100, 200, 500, 1000, 10000, 100000, 1000000)
n_replicates <- 10

# Set significance level
alpha <- 0.05

results <- data.frame(
  replicate = integer(),
  true_distribution = character(),
  n = integer(),
  R_plaw_vs_pois = numeric(),
  p_plaw_vs_pois = numeric(),
  R_plaw_vs_lnorm = numeric(),
  p_plaw_vs_lnorm = numeric(),
  R_plaw_vs_exp = numeric(),
  p_plaw_vs_exp = numeric(),
  R_pois_vs_exp = numeric(),
  p_pois_vs_exp = numeric(),
  R_pois_vs_lnorm = numeric(),
  p_pois_vs_lnorm = numeric(),
  concluded_distribution = character(),
  stringsAsFactors = FALSE
)

for (rep in 1:n_replicates) {
  cat("\n\n########## REPLICATE", rep, "##########\n")

  for (n in sample_sizes) {

    # Test 1: Data from power law distribution
    cat("\n=== Rep", rep, "- Testing n =", n, "from power law ===\n")
    x_plaw <- rpldis(n, xmin = 2, alpha = 3)
    x_plaw <- ceiling(x_plaw)

    # Fit models
    powerlawfit <- displ$new(x_plaw)
    powerlawfit$setPars(estimate_pars(powerlawfit))

    poissonfit <- dispois$new(x_plaw)
    poissonfit$setPars(estimate_pars(poissonfit))

    lnormfit <- dislnorm$new(x_plaw)
    lnormfit$setPars(estimate_pars(lnormfit))

    expfit <- disexp$new(x_plaw)
    expfit$setPars(estimate_pars(expfit))

    # Compare distributions
    plaw_pois <- safe_compare(powerlawfit, poissonfit)
    plaw_lnorm <- safe_compare(powerlawfit, lnormfit)
    plaw_exp <- safe_compare(powerlawfit, expfit)
    pois_exp <- safe_compare(poissonfit, expfit)
    pois_lnorm <- safe_compare(poissonfit, lnormfit)

    # Determine conclusion: power law should win (positive R, p < alpha)
    conclusion <- "ns"
    if (!is.na(plaw_pois$p_two_sided) && !is.na(plaw_lnorm$p_two_sided) && !is.na(plaw_exp$p_two_sided) &&
        plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic > 0 &&
        plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic > 0 &&
        plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic > 0) {
      conclusion <- "power_law"
    } else if (!is.na(plaw_pois$p_two_sided) && plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic < 0) {
      conclusion <- "poisson"
    } else if (!is.na(plaw_lnorm$p_two_sided) && plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic < 0) {
      conclusion <- "lognormal"
    } else if (!is.na(plaw_exp$p_two_sided) && plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic < 0) {
      conclusion <- "exponential"
    }

    results <- rbind(results, data.frame(
      replicate = rep,
      true_distribution = "power_law",
      n = n,
      R_plaw_vs_pois = plaw_pois$test_statistic,
      p_plaw_vs_pois = plaw_pois$p_two_sided,
      R_plaw_vs_lnorm = plaw_lnorm$test_statistic,
      p_plaw_vs_lnorm = plaw_lnorm$p_two_sided,
      R_plaw_vs_exp = plaw_exp$test_statistic,
      p_plaw_vs_exp = plaw_exp$p_two_sided,
      R_pois_vs_exp = pois_exp$test_statistic,
      p_pois_vs_exp = pois_exp$p_two_sided,
      R_pois_vs_lnorm = pois_lnorm$test_statistic,
      p_pois_vs_lnorm = pois_lnorm$p_two_sided,
      concluded_distribution = conclusion,
      stringsAsFactors = FALSE
    ))

    # Test 2: Data from binomial distribution
    cat("=== Rep", rep, "- Testing n =", n, "from binomial ===\n")
    x_binom <- rbinom(n, size = 100, prob = 0.2)  # reduced size for more reasonable values
    x_binom <- ceiling(x_binom)

    # Fit models
    powerlawfit <- displ$new(x_binom)
    powerlawfit$setPars(estimate_pars(powerlawfit))

    poissonfit <- dispois$new(x_binom)
    poissonfit$setPars(estimate_pars(poissonfit))

    lnormfit <- dislnorm$new(x_binom)
    lnormfit$setPars(estimate_pars(lnormfit))

    expfit <- disexp$new(x_binom)
    expfit$setPars(estimate_pars(expfit))

    # Compare distributions
    plaw_pois <- safe_compare(powerlawfit, poissonfit)
    plaw_lnorm <- safe_compare(powerlawfit, lnormfit)
    plaw_exp <- safe_compare(powerlawfit, expfit)
    pois_exp <- safe_compare(poissonfit, expfit)
    pois_lnorm <- safe_compare(poissonfit, lnormfit)

    # Determine conclusion: poisson should win (negative R for plaw_pois, p < alpha)
    conclusion <- "ns"
    if (!is.na(plaw_pois$p_two_sided) && plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic < 0) {
      # Poisson beats power law, now check if it beats others
      if (!is.na(pois_exp$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
          pois_exp$p_two_sided < alpha && pois_exp$test_statistic > 0 &&
          pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic > 0) {
        conclusion <- "poisson"
      } else {
        conclusion <- "ns"  # Poisson doesn't clearly beat all
      }
    } else if (!is.na(plaw_lnorm$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
               plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic < 0 &&
               pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic < 0) {
      conclusion <- "lognormal"
    } else if (!is.na(plaw_exp$p_two_sided) && !is.na(pois_exp$p_two_sided) &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic < 0 &&
               pois_exp$p_two_sided < alpha && pois_exp$test_statistic < 0) {
      conclusion <- "exponential"
    } else if (!is.na(plaw_pois$p_two_sided) && !is.na(plaw_lnorm$p_two_sided) && !is.na(plaw_exp$p_two_sided) &&
               plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic > 0 &&
               plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic > 0 &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic > 0) {
      conclusion <- "power_law"
    }

    results <- rbind(results, data.frame(
      replicate = rep,
      true_distribution = "binomial",
      n = n,
      R_plaw_vs_pois = plaw_pois$test_statistic,
      p_plaw_vs_pois = plaw_pois$p_two_sided,
      R_plaw_vs_lnorm = plaw_lnorm$test_statistic,
      p_plaw_vs_lnorm = plaw_lnorm$p_two_sided,
      R_plaw_vs_exp = plaw_exp$test_statistic,
      p_plaw_vs_exp = plaw_exp$p_two_sided,
      R_pois_vs_exp = pois_exp$test_statistic,
      p_pois_vs_exp = pois_exp$p_two_sided,
      R_pois_vs_lnorm = pois_lnorm$test_statistic,
      p_pois_vs_lnorm = pois_lnorm$p_two_sided,
      concluded_distribution = conclusion,
      stringsAsFactors = FALSE
    ))

    # Test 3: Data from lognormal distribution
    cat("=== Rep", rep, "- Testing n =", n, "from lognormal ===\n")
    x_lnorm <- rlnorm(n, meanlog = 0, sdlog = 3)
    x_lnorm <- ceiling(x_lnorm)

    # Fit models
    powerlawfit <- displ$new(x_lnorm)
    powerlawfit$setPars(estimate_pars(powerlawfit))

    poissonfit <- dispois$new(x_lnorm)
    poissonfit$setPars(estimate_pars(poissonfit))

    lnormfit <- dislnorm$new(x_lnorm)
    lnormfit$setPars(estimate_pars(lnormfit))

    expfit <- disexp$new(x_lnorm)
    expfit$setPars(estimate_pars(expfit))

    # Compare distributions
    plaw_pois <- safe_compare(powerlawfit, poissonfit)
    plaw_lnorm <- safe_compare(powerlawfit, lnormfit)
    plaw_exp <- safe_compare(powerlawfit, expfit)
    pois_exp <- safe_compare(poissonfit, expfit)
    pois_lnorm <- safe_compare(poissonfit, lnormfit)

    # Determine conclusion: lognormal should win (negative R for plaw_lnorm, p < alpha)
    conclusion <- "ns"
    if (!is.na(plaw_lnorm$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
        plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic < 0 &&
        pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic < 0) {
      conclusion <- "lognormal"
    } else if (!is.na(plaw_pois$p_two_sided) && !is.na(pois_exp$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
               plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic < 0 &&
               pois_exp$p_two_sided < alpha && pois_exp$test_statistic > 0 &&
               pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic > 0) {
      conclusion <- "poisson"
    } else if (!is.na(plaw_exp$p_two_sided) && !is.na(pois_exp$p_two_sided) &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic < 0 &&
               pois_exp$p_two_sided < alpha && pois_exp$test_statistic < 0) {
      conclusion <- "exponential"
    } else if (!is.na(plaw_pois$p_two_sided) && !is.na(plaw_lnorm$p_two_sided) && !is.na(plaw_exp$p_two_sided) &&
               plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic > 0 &&
               plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic > 0 &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic > 0) {
      conclusion <- "power_law"
    }

    results <- rbind(results, data.frame(
      replicate = rep,
      true_distribution = "lognormal",
      n = n,
      R_plaw_vs_pois = plaw_pois$test_statistic,
      p_plaw_vs_pois = plaw_pois$p_two_sided,
      R_plaw_vs_lnorm = plaw_lnorm$test_statistic,
      p_plaw_vs_lnorm = plaw_lnorm$p_two_sided,
      R_plaw_vs_exp = plaw_exp$test_statistic,
      p_plaw_vs_exp = plaw_exp$p_two_sided,
      R_pois_vs_exp = pois_exp$test_statistic,
      p_pois_vs_exp = pois_exp$p_two_sided,
      R_pois_vs_lnorm = pois_lnorm$test_statistic,
      p_pois_vs_lnorm = pois_lnorm$p_two_sided,
      concluded_distribution = conclusion,
      stringsAsFactors = FALSE
    ))

    # Test 4: Data from Poisson distribution
    cat("=== Rep", rep, "- Testing n =", n, "from Poisson ===\n")
    x_pois <- rpois(n, lambda = 20)
    x_pois <- ceiling(x_pois)

    # Fit models
    powerlawfit <- displ$new(x_pois)
    powerlawfit$setPars(estimate_pars(powerlawfit))

    poissonfit <- dispois$new(x_pois)
    poissonfit$setPars(estimate_pars(poissonfit))

    lnormfit <- dislnorm$new(x_pois)
    lnormfit$setPars(estimate_pars(lnormfit))

    expfit <- disexp$new(x_pois)
    expfit$setPars(estimate_pars(expfit))

    # Compare distributions
    plaw_pois <- safe_compare(powerlawfit, poissonfit)
    plaw_lnorm <- safe_compare(powerlawfit, lnormfit)
    plaw_exp <- safe_compare(powerlawfit, expfit)
    pois_exp <- safe_compare(poissonfit, expfit)
    pois_lnorm <- safe_compare(poissonfit, lnormfit)

    # Determine conclusion: poisson should win (negative R for plaw_pois, p < alpha)
    conclusion <- "ns"
    if (!is.na(plaw_pois$p_two_sided) && plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic < 0) {
      # Poisson beats power law, now check if it beats others
      if (!is.na(pois_exp$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
          pois_exp$p_two_sided < alpha && pois_exp$test_statistic > 0 &&
          pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic > 0) {
        conclusion <- "poisson"
      } else {
        conclusion <- "ns"  # Poisson doesn't clearly beat all
      }
    } else if (!is.na(plaw_lnorm$p_two_sided) && !is.na(pois_lnorm$p_two_sided) &&
               plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic < 0 &&
               pois_lnorm$p_two_sided < alpha && pois_lnorm$test_statistic < 0) {
      conclusion <- "lognormal"
    } else if (!is.na(plaw_exp$p_two_sided) && !is.na(pois_exp$p_two_sided) &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic < 0 &&
               pois_exp$p_two_sided < alpha && pois_exp$test_statistic < 0) {
      conclusion <- "exponential"
    } else if (!is.na(plaw_pois$p_two_sided) && !is.na(plaw_lnorm$p_two_sided) && !is.na(plaw_exp$p_two_sided) &&
               plaw_pois$p_two_sided < alpha && plaw_pois$test_statistic > 0 &&
               plaw_lnorm$p_two_sided < alpha && plaw_lnorm$test_statistic > 0 &&
               plaw_exp$p_two_sided < alpha && plaw_exp$test_statistic > 0) {
      conclusion <- "power_law"
    }

    results <- rbind(results, data.frame(
      replicate = rep,
      true_distribution = "poisson",
      n = n,
      R_plaw_vs_pois = plaw_pois$test_statistic,
      p_plaw_vs_pois = plaw_pois$p_two_sided,
      R_plaw_vs_lnorm = plaw_lnorm$test_statistic,
      p_plaw_vs_lnorm = plaw_lnorm$p_two_sided,
      R_plaw_vs_exp = plaw_exp$test_statistic,
      p_plaw_vs_exp = plaw_exp$p_two_sided,
      R_pois_vs_exp = pois_exp$test_statistic,
      p_pois_vs_exp = pois_exp$p_two_sided,
      R_pois_vs_lnorm = pois_lnorm$test_statistic,
      p_pois_vs_lnorm = pois_lnorm$p_two_sided,
      concluded_distribution = conclusion,
      stringsAsFactors = FALSE
    ))
  }
}

# Save all detailed results
write.csv(results, file = "distribution_test_results_all_replicates.csv", row.names = FALSE)
cat("\nAll results saved to distribution_test_results_all_replicates.csv\n")

# Create summary table with proportion correct
summary_table <- results %>%
  mutate(correct = case_when(
    true_distribution == "power_law" & concluded_distribution == "power_law" ~ 1,
    true_distribution == "binomial" & concluded_distribution == "poisson" ~ 1,
    true_distribution == "lognormal" & concluded_distribution == "lognormal" ~ 1,
    true_distribution == "poisson" & concluded_distribution == "poisson" ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(true_distribution, n) %>%
  summarise(
    n_replicates = n(),
    mean_R_plaw_vs_pois = mean(R_plaw_vs_pois),
    mean_p_plaw_vs_pois = mean(p_plaw_vs_pois),
    mean_R_plaw_vs_lnorm = mean(R_plaw_vs_lnorm),
    mean_p_plaw_vs_lnorm = mean(p_plaw_vs_lnorm),
    mean_R_plaw_vs_exp = mean(R_plaw_vs_exp),
    mean_p_plaw_vs_exp = mean(p_plaw_vs_exp),
    proportion_correct = mean(correct),
    .groups = "drop"
  )

cat("\n\n=== SUMMARY TABLE WITH PROPORTION CORRECT ===\n\n")
print(summary_table)

# Save summary table
write.csv(summary_table, file = "distribution_test_summary.csv", row.names = FALSE)
cat("\nSummary table saved to distribution_test_summary.csv\n")

# Create a simplified view
simple_summary <- summary_table %>%
  select(true_distribution, n, proportion_correct)

cat("\n\n=== SIMPLIFIED SUMMARY: PROPORTION CORRECT BY SAMPLE SIZE ===\n\n")
print(simple_summary, n = Inf)
