library(tidyverse); library(data.table); library(glue); library(brms); library(broom); library(bayesplot)
source("helpfuncs.R")

prior_informed_cohensd <- 0.28 # cohen's d
nchains <- 20
samples <- 2000

ddm <- fread("Data/ddm.csv")
ddm[condition == "control", conditionEC := -0.5]
ddm[condition == "deplete", conditionEC := 0.5]
ddm[congruency == "congruent", congruentEC := -0.5]
ddm[congruency == "incongruent", congruentEC := 0.5]




prior_coef <- expectedBeta(expected_d = -prior_informed_cohensd,
                           sd1 = ddm[condition == "control", sd(a, na.rm = T)], # REMEMBER TO CHANGE VARIALBE!
                           sd2 = ddm[condition == "deplete", sd(a, na.rm = T)])
prior_coef
get_prior(a ~ conditionEC + (1 | study/pNo), ddm[study == 1])
priors <- c(set_prior("normal(0, 1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"),
            set_prior(glue("normal({prior_coef}, {abs(prior_coef/2)})"), class = "b", coef = "conditionEC"),
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 1)", class = "sigma")) %>% print()

# fit model for each study (2-level model)
mbayes_bound_condition_congruency <- vector(mode = "list", length = 5)
for (i in 1:4) {
    mbayes_bound_condition_congruency[[i]] <- brm(a ~ conditionEC + congruentEC + (1 | pNo), data = ddm[study == i],
                             # family = lognormal(),
                             cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples,
                             file = glue("brms_models/mbayes_bound_condition_congruency_study{i}"),
                             control = list(adapt_delta = 0.99))
}

# fit 3-level model
mbayes_bound_condition_congruency[[5]] <- brm(a ~ conditionEC + congruentEC + (1 | study/pNo), data = ddm,
                             # family = lognormal(),
                             cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples * 3,
                             file = "brms_models/mbayes_bound_condition_congruency_study_all",
                             control = list(adapt_delta = 0.99))

# summarize model results
mbayes_bound_condition_congruency_results <- lapply(1:5, function(x) summarizebrms(mbayes_bound_condition_congruency[[x]], conf.method = "HPDinterval", effect = "conditionEC"))
manuscriptformat <- data.table(results = sapply(1:5, function(x) mbayes_bound_condition_congruency_results[[x]][effect == "manuscriptformat", result]))
manuscriptformat
tableformat <- lapply(1:5, function(x) formattable(mbayes_bound_condition_congruency_results[[x]]))
tableformat

mbayes_bound_condition_congruency_results <- lapply(1:5, function(x) summarizebrms(mbayes_bound_condition_congruency[[x]], conf.method = "HPDinterval", effect = "congruentEC"))
manuscriptformat <- data.table(results = sapply(1:5, function(x) mbayes_bound_condition_congruency_results[[x]][effect == "manuscriptformat", result]))
manuscriptformat

tableformat <- lapply(1:5, function(x) formattable(mbayes_bound_condition_congruency_results[[x]]))
tableformat





# NULL MODEL (no condition)
priors <- c(set_prior("normal(0, 1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"),
            # set_prior(glue("normal({prior_coef}, {abs(prior_coef/2)})"), class = "b", coef = "conditionEC"),
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 1)", class = "sigma")) %>% print()

# fit model for each study (2-level model)
mbayes_bound_condition_congruency_nocondition <- vector(mode = "list", length = 5)
for (i in 1:4) {
    mbayes_bound_condition_congruency_nocondition[[i]] <- brm(a ~ congruentEC + (1 | pNo), data = ddm[study == i],
                                                  # family = lognormal(),
                                                  cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples,
                                                  file = glue("brms_models/mbayes_bound_condition_congruency_study{i}_nocondition"),
                                                  control = list(adapt_delta = 0.99))
}

# fit 3-level model
mbayes_bound_condition_congruency_nocondition[[5]] <- brm(a ~ congruentEC + (1 | study/pNo), data = ddm,
                                              # family = lognormal(),
                                              cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples * 3,
                                              file = "brms_models/mbayes_bound_condition_congruency_study_all_nocondition",
                                              control = list(adapt_delta = 0.99))

# bridge sampling bayes factors
compute_bfs(mbayes_bound_condition_congruency, mbayes_bound_condition_congruency_nocondition)






# NULL MODEL (no congruency)
priors <- c(set_prior("normal(0, 1)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"),
            set_prior(glue("normal({prior_coef}, {abs(prior_coef/2)})"), class = "b", coef = "conditionEC"),
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 1)", class = "sigma")) %>% print()

# fit model for each study (2-level model)
mbayes_bound_condition_congruency_nocongruency <- vector(mode = "list", length = 5)
for (i in 1:4) {
    mbayes_bound_condition_congruency_nocongruency[[i]] <- brm(a ~ conditionEC + (1 | pNo), data = ddm[study == i],
                                                               # family = lognormal(),
                                                               cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples,
                                                               file = glue("brms_models/mbayes_bound_condition_congruency_study{i}_nocongruency"),
                                                               control = list(adapt_delta = 0.99))
}

# fit 3-level model
mbayes_bound_condition_congruency_nocongruency[[5]] <- brm(a ~ conditionEC + (1 | study/pNo), data = ddm,
                                                           # family = lognormal(),
                                                           cores = nchains, chains = nchains, sample_prior = TRUE, save_pars = save_pars(all = TRUE), prior = priors, iter = samples * 3,
                                                           file = "brms_models/mbayes_bound_condition_congruency_study_all_nocongruency",
                                                           control = list(adapt_delta = 0.99))

# bridge sampling bayes factors
compute_bfs(mbayes_bound_condition_congruency, mbayes_bound_condition_congruency_nocongruency)

# ===== Figure: Posterior Distribution of Condition Effect on Boundary =====

library(posterior)

# --- Extract posterior samples from the 3-level combined model (index 5) ---
posterior_samples <- as_draws_matrix(
    mbayes_bound_condition_congruency[[5]],
    variable = "b_conditionEC"
)[, "b_conditionEC"]

# --- Calculate posterior density ---
posterior_density <- density(posterior_samples)
posterior_df <- data.frame(
    estimate = posterior_density$x,
    density = posterior_density$y
)

# --- Calculate prior density (informed normal) ---
prior_mean <- prior_coef
prior_sd <- abs(prior_coef / 2)
prior_x <- seq(min(posterior_density$x), max(posterior_density$x), length.out = 512)
prior_df <- data.frame(
    estimate = prior_x,
    density = dnorm(prior_x, mean = prior_mean, sd = prior_sd)
)

# --- Extract BF from 3-level model ---
hypo_3level <- brms::hypothesis(
    mbayes_bound_condition_congruency[[5]],
    "conditionEC = 0"
)
bf_value <- round(1 / hypo_3level$hypothesis$Evid.Ratio, 2)

# --- Create the plot ---
p_boundary <- ggplot() +
    # Prior distribution (light blue, broader)
    geom_area(data = prior_df, aes(x = estimate, y = density),
              fill = "#A6CEE3", alpha = 0.6, color = "#A6CEE3", linewidth = 0.3) +
    # Posterior distribution (dark blue, narrower)
    geom_area(data = posterior_df, aes(x = estimate, y = density),
              fill = "#1F78B4", alpha = 0.8, color = "#1F78B4", linewidth = 0.3) +
    # Reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    # BF annotation
    annotate("text", x = min(posterior_df$estimate) * 0.8,
             y = max(posterior_df$density) * 0.9,
             label = paste0("BF = ", bf_value),
             size = 5, hjust = 0) +
    # Labels and theme
    labs(title = "Effect of Condition on Boundary",
         x = "Estimate", y = "Density") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Save the figure ---
dir.create("Figures", showWarnings = FALSE)
ggsave("Figures/effect_condition_on_boundary.png", p_boundary,
       width = 8, height = 6, dpi = 300)
ggsave("Figures/effect_condition_on_boundary.pdf", p_boundary,
       width = 8, height = 6)

# Display
print(p_boundary)
