## =============================================================================
## Diagnostic and additional code (not part of the SI document)
##
## Materials and participant descriptives, parametric (full-random-effects)
## linear models, the maximal-exclusion robustness check, and diagnostic
## comparison figures.
## =============================================================================

## Set working directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}
setwd(get_script_dir())

## Libraries, plotting theme, helper functions, `crit`, the by-item data frames
## and the shared plotting helpers all live in setup.R.
source("setup.R")

## Not part of the SI document. Materials & participant descriptives, parametric
## (full-random-effects) linear models, the maximal-exclusion robustness check,
## and diagnostic comparison figures.
## #############################################################################

## ---- Materials Information ----
### ---- Sentence Length Across Languages ----
lengths <- crit %>%
  select(language, item, n_words) %>%
  distinct() %>%
  mutate(
    language = case_when(
      language == "da" ~ "Danish",
      language == "de" ~ "German",
      language == "en" ~ "English",
      language == "hi" ~ "Hindi",
      language == "fr" ~ "French",
      language == "ko" ~ "Korean",
      language == "ru" ~ "Russian",
      language == "zh" ~ "Mandarin",
      language == "en_exp" ~ "English\n(controlled)"
    ),
    language = factor(
      language,
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
      )
  )


lengths %>%
  group_by(language) %>% 
  summarize(
    median_length = median(n_words),
    min_length = min(n_words),
    max_length = max(n_words),
    sd_length = sd(n_words),
  )

ggplot(lengths) +
  geom_histogram(
    aes(x=n_words)
  ) + 
  facet_wrap(~language, scales="free_y") + 
  scale_y_continuous(breaks = scales::breaks_width(5))

ggsave("figures/sentence_length_histograms.png", width=10, height=6)

### ---- Comprehension Accuracy by Edit Distance ----
## Descriptive check: how accurately did participants answer the comprehension
## question that followed each sentence, as a function of how corrupted that
## sentence was? Restricted to the three error conditions in the design (0, 1, and
## 3 adjacent-token exchanges); the rare intermediate DL distance of 2 is dropped
## as it is only sparsely populated. Scored the same way as the exclusion
## criterion above (correct = literal_response == response), on the same retained
## sample (participants below the 80% threshold already excluded). Not part of the
## SI document.
comprehension_accuracy <- comprehension %>%
  filter(!participant_id %in% exclusions,
         dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(language, dl_distance) %>%
  summarize(n = n(), k = sum(correct, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    accuracy = k / n,
    ## Wilson 95% confidence interval for a binomial proportion (well behaved
    ## near the ceiling, where a normal-approximation interval would exceed 1)
    z = 1.96,
    wilson_center = (accuracy + z^2 / (2 * n)) / (1 + z^2 / n),
    wilson_half   = z * sqrt(accuracy * (1 - accuracy) / n + z^2 / (4 * n^2)) /
                    (1 + z^2 / n),
    lower = wilson_center - wilson_half,
    upper = wilson_center + wilson_half,
    language = label_language(language)
  )

## Pooled accuracy per DL distance, printed for reference
comprehension %>%
  filter(!participant_id %in% exclusions, dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(dl_distance) %>%
  summarize(accuracy = mean(correct, na.rm = TRUE), n = n(), .groups = "drop")

ggplot(comprehension_accuracy, aes(x = dl_distance, y = accuracy)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_line(color = "grey60") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(size = 2) +
  facet_wrap(~language) +
  scale_x_continuous(breaks = c(0, 1, 3), name = "Damerau-Levenshtein distance") +
  scale_y_continuous(name = "Comprehension accuracy", limits = c(0.5, 1)) +
  theme(panel.grid.major.y = element_line(color = "grey90"))

ggsave("figures/comprehension_accuracy.png", width = 10, height = 6)

### ---- Count Damerau-Levenshtein distance for each language ----
distance_counts <- crit %>%
  select(language, item, condition, dl_distance, n_words) %>%
  distinct() %>% mutate(
    language = case_when(
      language == "da" ~ "Danish",
      language == "de" ~ "German",
      language == "en" ~ "English",
      language == "hi" ~ "Hindi",
      language == "fr" ~ "French",
      language == "ko" ~ "Korean",
      language == "ru" ~ "Russian",
      language == "zh" ~ "Mandarin",
      language == "en_exp" ~ "English\n(controlled)"
    ),
    language = factor(
      language,
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
    )
  )

table(distance_counts$language)
table(distance_counts$item, distance_counts$language)

sink("outputs/dl_distance_counts.txt")
print("Original Sentences")
table(
  distance_counts$dl_distance[distance_counts$condition == "original"], 
  distance_counts$language[distance_counts$condition == "original"]
) %>%
  print(., n=10^4)

print("1 Exchange")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_1"], 
  distance_counts$language[distance_counts$condition == "swap_1"]
) %>%
  print(., n=10^4)

print("3 Exchanges")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_3"], 
  distance_counts$language[distance_counts$condition == "swap_3"]
) %>%
  print(., n=10^4)

print("Shuffle")
table(
  distance_counts$dl_distance[distance_counts$condition == "random"], 
  distance_counts$language[distance_counts$condition == "random"]
) %>%
  print(., n=10^4)

sink()

## ---- Demographic Information ----
demographics <- dat %>%
  select(
    participant_id, country_birth, country_birth_text, first_language, 
    first_language_text,second_language , second_language_text, 
    language_contexts, english_first_language, country_of_birth, 
    french_first_language,  french_contexts, language
  ) %>%
  distinct

## Danish
danish_demographics <- demographics %>% filter(language == "da")
table(danish_demographics$first_language, danish_demographics$country_birth)

danish_demographics %>%
  filter(country_birth == "Andet")

## English
english_demographics <- demographics %>% filter(language == "en")

table(english_demographics$english_first_language, 
      english_demographics$country_of_birth)

## English (Length-Controlled Experiment)
english_exp_demographics <- demographics %>% filter(language == "en_exp")

table(english_exp_demographics$english_first_language, 
      english_exp_demographics$country_of_birth,
      useNA = "always"
)

## French
french_demographics <- demographics %>% filter(language == "fr")

table(
  french_demographics$french_first_language,
  french_demographics$country_of_birth
)

french_demographics %>%
  filter(country_of_birth == "Ailleurs") %>%
  pull(french_contexts)

## German
german_demographics <- demographics %>% filter(language == "de")

table(
  german_demographics$first_language,
  german_demographics$country_birth
)

german_demographics %>% 
  filter(country_birth == "Anderes") %>%
  pull(country_birth_text)

## Hindi
hindi_demographics <- demographics %>% filter(language == "hi")

table(
  hindi_demographics$first_language,
  hindi_demographics$country_birth
)

table(hindi_demographics$first_language_text)
table(hindi_demographics$country_birth_text, hindi_demographics$first_language)

## Korean
korean_demographics <- demographics %>% filter(language == "ko")

table(
  korean_demographics$first_language,
  korean_demographics$country_of_birth
)

## Chinese
chinese_demographics <- demographics %>% filter(language == "zh")

table(chinese_demographics$english_first_language)
table(chinese_demographics$country_of_birth)

table(
  chinese_demographics$english_first_language,
  chinese_demographics$country_of_birth
)

## Russian
russian_demographics <- demographics %>% filter(language == "ru")

table(
  russian_demographics$first_language,
  russian_demographics$country_birth
)

table(russian_demographics$country_birth,
      russian_demographics$country_birth_text)



## ---- Linear Models ----
### ---- All Materials ----
### ---- Collapsing Across Language ----
if (file.exists("models/linear_model_fits.Rdata")) {
  load("models/linear_model_fits.Rdata")
} else {
  
  model_grammatical_linear <- lmer(
    response_num ~ grammatical + 
      (grammatical|participant_id) + (grammatical|item) + (grammatical|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_distance_linear <- lmer(
    response_num ~ dl_distance + 
      (dl_distance|participant_id) + (dl_distance|item) + (dl_distance|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_weighted_distance_linear <- lmer(
    response_num ~ weighted_dl_distance + n_words + 
      (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item) +
      (weighted_dl_distance + n_words|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  save(
    model_grammatical_linear, 
    model_distance_linear, 
    model_weighted_distance_linear,
    file="models/linear_model_fits.Rdata"
  )
}

## Extract log likelihoods
model_grammatical_linear
model_distance_linear
model_weighted_distance_linear

logLik(model_grammatical_linear) %>% as.numeric()
logLik(model_distance_linear) %>% as.numeric()
logLik(model_weighted_distance_linear) %>% as.numeric()

AIC(model_grammatical_linear)
AIC(model_distance_linear)
AIC(model_weighted_distance_linear)

BIC(model_grammatical_linear)
BIC(model_distance_linear)
BIC(model_weighted_distance_linear)


main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    logLik(model_grammatical_linear) %>% as.numeric(),
    logLik(model_distance_linear) %>% as.numeric(),
    logLik(model_weighted_distance_linear) %>% as.numeric()
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("grammatical", "distance", "weighted_distance"))
  )

p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Languages",
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40"),
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods


### ---- Analysis By Language ----
if (file.exists("models/models_by_language_linear.RData")) {
  load("models/models_by_language_linear.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit$language)) {
    ## Get specific data
    crit_language_specific <- crit %>%
      filter(language == language_iter)
    
    ## Collapse to the 72 base items for English (controlled): the item id
    ## encodes the length version (en_experimental_62m), so strip the trailing
    ## s/m/l. Length then varies within item.
    if (language_iter == "en_exp") {
      crit_language_specific <- crit_language_specific %>%
        mutate(
          item = as.factor(gsub("[sml]$", "", item))
        )
      
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + 
          (weighted_dl_distance + n_words|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    } else {
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    }
    
    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )
    
    stats_models_by_language <- stats_models_by_language %>% 
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = grammatical %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = distance %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = weighted_distance %>% logLik() %>% as.numeric()
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
    mutate(
      language = case_when(
        language == "da" ~ "Danish",
        language == "de" ~ "German",
        language == "en" ~ "English",
        language == "hi" ~ "Hindi",
        language == "fr" ~ "French",
        language == "ko" ~ "Korean",
        language == "ru" ~ "Russian",
        language == "zh" ~ "Mandarin",
        language == "en_exp" ~ "English\n(controlled)"
      ),
      language = factor(
        language,
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_linear.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40")
  ) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_linear.png", width=8, height=8)

### ---- DL <= 3 ----
if (file.exists("models/linear_model_fits_dl3_interactions.RData")) {
  load("models/linear_model_fits_dl3_interactions.RData")
} else {
  fit_lmer_interaction <- lmer(
    response_num ~ dl_distance * n_words +
      (dl_distance * n_words|participant_id) + (dl_distance|item) +
      (dl_distance * n_words|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ) %>%
      mutate(
        n_words = scale(n_words)[,1],
        dl_distance = scale(dl_distance)[,1]
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  fit_lmers_by_language <- lapply(
    levels(crit$language),
    function(lang) {
      print(paste0("Fitting for language: ", lang))
      lmer(
        response_num ~ dl_distance * n_words +
          (dl_distance * n_words|participant_id) + (dl_distance|item),
        data=crit %>%
          filter(
            language == lang,
            dl_distance <= 3
          ) %>%
          mutate(
            n_words = scale(n_words)[,1],
            dl_distance = scale(dl_distance)[,1]
          ),
        control=lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 1e8)
        )
      )
    }
  )
  
  names(fit_lmers_by_language) <- levels(crit$language)
  
  ## Collapse the item id to the 72 base items (see Section 3): length varies
  ## within item only after the trailing s/m/l is stripped.
  fit_lmers_by_language$en_exp <- lmer(
    response_num ~ dl_distance * n_words +
      (dl_distance * n_words|participant_id) + (dl_distance * n_words|item),
    data=crit %>%
      filter(
        language == "en_exp",
        dl_distance <= 3
      ) %>%
      mutate(
        item = as.factor(gsub("[sml]$", "", item)),
        n_words = scale(n_words)[,1],
        dl_distance = scale(dl_distance)[,1]
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e8)
    )
  ) 
  
  lapply(fit_lmers_by_language, function(langfit)
    summary(langfit)$coef
         )
  
  save(
    fit_lmer_interaction,
    fit_lmers_by_language,
    file="models/linear_model_fits_dl3_interactions.RData"
  )  
}

summary(fit_lmer_interaction)


if (file.exists("models/linear_model_fits_dl3.Rdata")) {
  load("models/linear_model_fits_dl3.Rdata")
} else {
  
  model_grammatical_linear_dl3 <- lmer(
    response_num ~ grammatical + 
      (grammatical|participant_id) + (grammatical|item) + (grammatical|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_distance_linear_dl3 <- lmer(
    response_num ~ dl_distance + 
      (dl_distance|participant_id) + (dl_distance|item) + (dl_distance|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_weighted_distance_linear_dl3 <- lmer(
    response_num ~ weighted_dl_distance + n_words + 
      (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item) +
      (weighted_dl_distance + n_words|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  save(
    model_grammatical_linear_dl3,
    model_distance_linear_dl3,
    model_weighted_distance_linear_dl3,
    file="models/linear_model_fits_dl3.Rdata"
  )
}

## Extract log likelihoods
model_grammatical_linear_dl3
model_distance_linear_dl3
model_weighted_distance_linear_dl3

logLik(model_grammatical_linear_dl3) %>% as.numeric()
logLik(model_distance_linear_dl3) %>% as.numeric()
logLik(model_weighted_distance_linear_dl3) %>% as.numeric()

AIC(model_grammatical_linear_dl3)
AIC(model_distance_linear_dl3)
AIC(model_weighted_distance_linear_dl3)

BIC(model_grammatical_linear_dl3)
BIC(model_distance_linear_dl3)
BIC(model_weighted_distance_linear_dl3)


main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    logLik(model_grammatical_linear_dl3) %>% as.numeric(),
    logLik(model_distance_linear_dl3) %>% as.numeric(),
    logLik(model_weighted_distance_linear_dl3) %>% as.numeric()
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("grammatical", "distance", "weighted_distance"))
  )

p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Languages",
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40"),
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods


### ---- Analysis By Language ----
if (file.exists("models/models_by_language_linear_dl3.RData")) {
  load("models/models_by_language_linear_dl3.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit_by_items_full$language)) {
    ## Get specific data
    crit_language_specific <- crit %>%
      filter(language == language_iter,
             dl_distance <= 3)
    
    ## Collapse to the 72 base items for English (controlled): the item id
    ## encodes the length version (en_experimental_62m), so strip the trailing
    ## s/m/l. Length then varies within item.
    if (language_iter == "en_exp") {
      crit_language_specific <- crit_language_specific %>%
        mutate(
          item = as.factor(gsub("[sml]$", "", item))
        )
      
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + 
          (weighted_dl_distance + n_words|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    } else {
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    }
    
    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )
    
    stats_models_by_language <- stats_models_by_language %>% 
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = grammatical %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = distance %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = weighted_distance %>% logLik() %>% as.numeric()
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
    mutate(
      language = case_when(
        language == "da" ~ "Danish",
        language == "de" ~ "German",
        language == "en" ~ "English",
        language == "hi" ~ "Hindi",
        language == "fr" ~ "French",
        language == "ko" ~ "Korean",
        language == "ru" ~ "Russian",
        language == "zh" ~ "Mandarin",
        language == "en_exp" ~ "English\n(controlled)"
      ),
      language = factor(
        language,
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_linear_dl3.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40")
  ) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_linear_dl3.png", width=8, height=8)


## ---- Main Analysis with Maximal Exclusions ----
crit_exclusions <- crit %>%
  filter(
    first_language %in% c(
      "Dansk",
      "Deutsch",
      "हिंदी",
      "한국어",
      "Русский"
    ) |
      english_first_language %in% c("Yes", "是") |
      french_first_language == "Oui"
  )

table(crit_exclusions$language)

crit_exclusions_by_items <- crit_exclusions %>% 
  mutate(
    weighted_dl_distance = crit_exclusions$dl_distance / crit_exclusions$n_words,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response),
    response_num = as.numeric(as.character(response)),
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance, 
           grammatical, sum_surprisal, condition, weighted_dl_distance_surprisal) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  )

### ---- Across Languages ----
if (file.exists("models/main_model_fits_maximal_exclusions.RData")) {
  load("models/main_model_fits_maximal_exclusions.RData")
} else {
  
  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +# + # Fixed effects
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"),# RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects
      s(language, bs="re") + s(language, weighted_dl_distance, bs="re") + s(language, n_words, bs="re"), # RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  save(
    mod_grammatical, 
    mod_distance, 
    mod_weighted_distance,
    file="models/main_model_fits_maximal_exclusions.RData"
  )
}

## Inspect models
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv

## Extract the marginal log-likelihood of each model (see the note above: scam's
## own logLik()/AIC()/BIC() are conditional on the fitted random effects).
marginal_logLik(mod_grammatical)
marginal_logLik(mod_distance)
marginal_logLik(mod_weighted_distance)

main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    marginal_logLik(mod_grammatical),
    marginal_logLik(mod_distance),
    marginal_logLik(mod_weighted_distance)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("grammatical", "distance", "weighted_distance"))
  )

p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Languages",
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40"),
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods

### ---- Analysis By Language: Length ----
if (file.exists("models/models_by_language_maximal_exclusions.RData")) {
  load("models/models_by_language_maximal_exclusions.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_language_specific <- crit_exclusions_by_items %>%
      filter(language == language_iter)
    
    ## Fit models
    grammatical <- scam(
      response_num ~ grammatical,
      data=crit_language_specific
    )
    
    distance <- scam(
      response_num ~ s(dl_distance, bs="mpd"),
      data=crit_language_specific
    )
    
    weighted_distance <- scam(
      response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd"),
      data=crit_language_specific
    )

    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )

    ## Unlike the other by-language blocks, this tibble is NOT rebuilt after the
    ## cache guard, so it must carry the marginal value directly. These models
    ## carry no random effects, so marginal_logLik() returns the same number as
    ## logLik() -- but using it keeps the figure correct if REs are ever added.
    stats_models_by_language <- stats_models_by_language %>%
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = marginal_logLik(grammatical)
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = marginal_logLik(distance)
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = marginal_logLik(weighted_distance)
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
    mutate(
      language = case_when(
        language == "da" ~ "Danish",
        language == "de" ~ "German",
        language == "en" ~ "English",
        language == "hi" ~ "Hindi",
        language == "fr" ~ "French",
        language == "ko" ~ "Korean",
        language == "ru" ~ "Russian",
        language == "zh" ~ "Mandarin",
        language == "en_exp" ~ "English\n(controlled)"
      ),
      language = factor(
        language,
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_maximal_exclusions.RData"
  )
}

stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40")
  ) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_max_exclusions.png", width=8, height=8)


## ---- Figures ----
#### ---- Helper Functions 
YLIM <- c(1, 7)

DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")
DL_BUCKET_LEVELS <- c("0", "1-2", "3", "4+")

bucket_dl <- function(dl_distance) {
  dplyr::case_when(
    dl_distance == 0 ~ "0",
    dl_distance %in% 1:2 ~ "1-2",
    dl_distance == 3 ~ "3",
    dl_distance >= 4 ~ "4+"
  ) %>%
    factor(levels = DL_BUCKET_LEVELS)
}

scale_dl_color <- function() {
  ggplot2::scale_color_manual(
    name   = "Damerau-Levenshtein Distance",
    values = DL_COLORS,
    breaks = DL_BUCKET_LEVELS,
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 3)
    )
  )
}


get_prediction_grid_scam <- function(
    lang,
    model_type,
    n_words=4:43, 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = n_words,
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(models_by_language[[lang]][[model_type]]$model$item[1], levels=levels(models_by_language[[lang]][[model_type]]$model$item))
    )
}

get_prediction_grid_lmer <- function(
    lang,
    model_type,
    n_words=4:43, 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = n_words,
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(models_by_language[[lang]][[model_type]]@flist[1]$item[1], levels=levels(models_by_language[[lang]][[model_type]]@flist[1]$item)),
      participant_id=factor(models_by_language[[lang]][[model_type]]@flist[2]$participant_id[1], levels=levels(models_by_language[[lang]][[model_type]]@flist[2]$participant_id))
    )
}

### ---- Baseline data for plotting (item means) ----
plot_data_by_item <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response),
    response_num = as.numeric(as.character(response)),
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance, 
           grammatical, condition) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  dplyr::group_by(language, item, n_words, dl_distance) %>%
  dplyr::summarise(response = mean(response_num), .groups = "drop") %>%
  dplyr::mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )


p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = n_words, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence Length in Words", y="Rating") +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    
    # put facet strips on the outside + rotate text
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  scale_dl_color()


### ---- Model-based plots ----
#### ---- lmer predictions ----
##### ----- Communication, full data ----
load("models/models_by_language_linear.RData")

## Communication
fulldata_communication_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$weighted_distance,
      newdata = get_prediction_grid_lmer(lang, "weighted_distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

fulldata_communication_plot_data <- bind_rows(
  lapply(names(fulldata_communication_predictions), function(lang) {
    preds <- fulldata_communication_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "weighted_distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_full_communication <- p_baseline +
  geom_ribbon(
    data=fulldata_communication_plot_data %>%
      mutate(
        lower = response - 2*se,
        upper = response + 2*se,
        upper = ifelse(upper > 7.5, 7.5, upper),
      ),
    aes(
      x = n_words,
      ymin = lower,
      ymax = upper,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=fulldata_communication_plot_data %>%
      mutate(
        response = ifelse(response > 7.5, 7.5, response)
      ),
    aes(
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
    y = ""
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    
    # put facet strips on the outside + rotate text
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  ggtitle("Communication Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(YLIM[1], YLIM[2], by=2),
    labels = seq(YLIM[1], YLIM[2], by=2)
  ) +
  ggplot2::guides(color = "none")

p_full_communication


## Gradient
fulldata_gradient_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$distance,
      newdata = get_prediction_grid_lmer(lang, "distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

fulldata_gradient_predictions_plot_data <- bind_rows(
  lapply(names(fulldata_gradient_predictions), function(lang) {
    preds <- fulldata_gradient_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_full_gradient <- p_baseline +
  geom_ribbon(
    data=fulldata_gradient_predictions_plot_data %>%
      mutate(
        lower = response - 2*se,
        upper = response + 2*se,
        upper = ifelse(upper > 7.5, 7.5, upper),
      ),
    aes(
      x = n_words,
      ymin = lower,
      ymax = upper,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=fulldata_gradient_predictions_plot_data %>%
      mutate(
        response = ifelse(response > 7.5, 7.5, response)
      ),
    aes(
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Gradient Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_full_gradient

## with_shuffled
load("models/models_by_language_linear_dl3.RData")
with_shuffled_communication_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$weighted_distance,
      newdata = get_prediction_grid_lmer(lang, "weighted_distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

with_shuffled_communication_plot_data <- bind_rows(
  lapply(names(with_shuffled_communication_predictions), function(lang) {
    preds <- with_shuffled_communication_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "weighted_distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_with_shuffled_communication <- p_baseline +
  geom_ribbon(
    data=with_shuffled_communication_plot_data,
    aes(
      x = n_words,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=with_shuffled_communication_plot_data,
    aes(
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Communication Model (DL ≤ 3)") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_with_shuffled_communication

## Gradient
with_shuffled_gradient_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$distance,
      newdata = get_prediction_grid_lmer(lang, "distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

with_shuffled_gradient_plot_data <- bind_rows(
  lapply(names(with_shuffled_gradient_predictions), function(lang) {
    preds <- with_shuffled_gradient_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_with_shuffled_gradient <- p_baseline +
  geom_ribbon(
    data=with_shuffled_gradient_plot_data,
    aes(
      x = n_words,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=with_shuffled_gradient_plot_data,
    aes(
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Gradient Model (DL ≤ 3)") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_with_shuffled_gradient

## Combine and save linear figures
linear_combined_plot <- (p_full_communication |
  p_full_gradient |
  p_with_shuffled_communication |
  p_with_shuffled_gradient) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )

## Save as pdf
ggsave(
  filename = "figures/main_figure_linear_models.pdf",
  plot     = linear_combined_plot,
  width    = 25,
  height   = 35,
  units    = "in"
)
