# Wheat field-trial analysis: single-environment and multi-environment
# Required packages: lme4, emmeans, multcompView

library(lme4)
library(emmeans)
library(multcompView)

# Create output folders if needed
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Import and clean the downloaded data
wheat <- read.csv("data/George_et_al_2019.csv", stringsAsFactors = FALSE)
wheat$YIELD[wheat$YIELD == "#VALUE!"] <- NA
wheat$YIELD <- as.numeric(wheat$YIELD)
wheat <- subset(wheat, !is.na(YIELD))
wheat$Entry_Number <- factor(wheat$Entry_Number)
wheat$BLOCK <- factor(wheat$BLOCK)
wheat$Environment <- interaction(wheat$YEAR, wheat$LOCATION, sep = ":", drop = TRUE)

cat("Usable plot records:", nrow(wheat), "\n")
cat("Entries:", nlevels(wheat$Entry_Number), "\n")
cat("Year-location environments:", nlevels(wheat$Environment), "\n")

# Part 1: one balanced trial, for a clear beginner-level comparison
trial <- subset(wheat, YEAR == 2010 & LOCATION == "COLUSA")
cat("\nCOLUSA 2010 records per block:\n")
print(table(trial$BLOCK))

single_model <- aov(YIELD ~ Entry_Number + BLOCK, data = trial)
cat("\nCOLUSA 2010 ANOVA:\n")
print(summary(single_model))


# Save model diagnostics and flag observations with |standardized residual| >= 2
pdf("figures/model_diagnostics.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
plot(single_model)
dev.off()

trial$Standardized_Residual <- rstandard(single_model)
observations_to_review <- subset(trial, abs(Standardized_Residual) >= 2)
observations_to_review <- observations_to_review[c(
  "Entry_Number", "BLOCK", "YIELD", "Standardized_Residual"
)]
write.csv(observations_to_review, "results/observations_to_review.csv", row.names = FALSE)

single_tukey <- TukeyHSD(single_model, "Entry_Number")$Entry_Number
tukey_table <- data.frame(
  Comparison = rownames(single_tukey),
  single_tukey,
  row.names = NULL,
  check.names = FALSE
)
write.csv(tukey_table, "results/tukey_colusa_2010.csv", row.names = FALSE)

tukey_p <- single_tukey[, "p adj"]
names(tukey_p) <- rownames(single_tukey)
tukey_letters <- multcompLetters(tukey_p, threshold = 0.05)$Letters

single_means <- aggregate(YIELD ~ Entry_Number, data = trial, FUN = mean)
single_means <- single_means[order(-single_means$YIELD), ]
residual_mse <- summary(single_model)[[1]]["Residuals", "Mean Sq"]
replicates <- as.numeric(table(trial$Entry_Number)[as.character(single_means$Entry_Number)])
single_means$SE <- sqrt(residual_mse / replicates)
t_critical <- qt(0.975, df.residual(single_model))
single_means$CI_lower <- single_means$YIELD - t_critical * single_means$SE
single_means$CI_upper <- single_means$YIELD + t_critical * single_means$SE
single_means$Tukey_group <- tukey_letters[as.character(single_means$Entry_Number)]
write.csv(single_means, "results/colusa_2010_entry_means.csv", row.names = FALSE)

# Plot all 29 entries: dot estimates, 95% confidence intervals, and Tukey groups
plot_data <- single_means[rev(seq_len(nrow(single_means))), ]
y <- seq_len(nrow(plot_data))
letter_space <- max(diff(range(plot_data$CI_lower, plot_data$CI_upper)) * 0.03, 100)
plot_right <- max(plot_data$CI_upper) + letter_space * 4
png("figures/colusa_2010_all_entries.png", width = 1500, height = 1100)
par(mar = c(5, 6, 4, 2))
plot(
  plot_data$YIELD, y,
  type = "n",
  xlim = c(min(plot_data$CI_lower) - 100, plot_right),
  ylim = c(0.5, nrow(plot_data) + 0.5),
  yaxt = "n",
  xlab = "Mean yield (kg/ha)",
  ylab = "Entry number",
  main = "COLUSA 2010: all entries with 95% CIs and Tukey groups"
)
axis(2, at = y, labels = plot_data$Entry_Number, las = 1)
segments(plot_data$CI_lower, y, plot_data$CI_upper, y)
points(plot_data$YIELD, y, pch = 19)
text(plot_data$CI_upper + letter_space, y, labels = plot_data$Tukey_group, adj = 0)
text(plot_right, 0.1, labels = "Shared letters = not significantly different (Tukey, 5%)", adj = 1, xpd = NA, cex = 0.8)
dev.off()

# Part 2: full multi-environment mixed model
# Entry is fixed; environment, blocks within environment, and entry-by-environment
# deviations are random effects.
multi_model <- lmer(
  YIELD ~ Entry_Number +
    (1 | Environment) +
    (1 | Environment:BLOCK) +
    (1 | Entry_Number:Environment),
  data = wheat,
  REML = TRUE
)

entry_emmeans <- as.data.frame(emmeans(multi_model, ~ Entry_Number, lmer.df = "asymptotic"))
entry_emmeans <- entry_emmeans[order(-entry_emmeans$emmean), ]

# Summarize observed within-environment rank variation as a simple GxE view
entry_environment_means <- aggregate(
  YIELD ~ Entry_Number + Environment,
  data = wheat,
  FUN = mean
)
entry_environment_means$Within_environment_rank <- ave(
  entry_environment_means$YIELD,
  entry_environment_means$Environment,
  FUN = function(x) rank(-x, ties.method = "average")
)
rank_summaries <- lapply(split(entry_environment_means, entry_environment_means$Entry_Number), function(x) {
  data.frame(
    Entry_Number = as.character(x$Entry_Number[1]),
    Environments_tested = nrow(x),
    Mean_within_environment_rank = mean(x$Within_environment_rank),
    SD_within_environment_rank = if (nrow(x) > 1) sd(x$Within_environment_rank) else NA_real_
  )
})
rank_summary <- do.call(rbind, rank_summaries)
entry_emmeans$Entry_Number <- as.character(entry_emmeans$Entry_Number)
entry_emmeans <- merge(entry_emmeans, rank_summary, by = "Entry_Number", all.x = TRUE, sort = FALSE)
entry_emmeans <- entry_emmeans[order(-entry_emmeans$emmean), ]
write.csv(entry_emmeans, "results/multi_environment_entry_means.csv", row.names = FALSE)
write.csv(entry_environment_means, "results/entry_environment_means.csv", row.names = FALSE)

variance_components <- as.data.frame(VarCorr(multi_model))
write.csv(variance_components, "results/multi_environment_variance_components.csv", row.names = FALSE)

cat("\nMulti-environment mixed-model variance components:\n")
print(variance_components[, c("grp", "vcov", "sdcor")])
cat("\nTop 10 model-adjusted entry means across environments:\n")
print(head(entry_emmeans[, c("Entry_Number", "emmean", "SE", "asymp.LCL", "asymp.UCL", "Environments_tested", "SD_within_environment_rank")], 10))
cat("\nMixed model singular fit:", isSingular(multi_model), "\n")

# Save diagnostic plots and the software session details
pdf("figures/multi_environment_diagnostics.pdf", width = 8, height = 6)
plot(fitted(multi_model), resid(multi_model),
     xlab = "Fitted yield", ylab = "Conditional residual",
     main = "Multi-environment model: residuals vs fitted")
abline(h = 0, lty = 2, col = "red")
dev.off()

capture.output(sessionInfo(), file = "results/session_info.txt")

# Remove the earlier top-ten-only chart if it remains from a previous run
if (file.exists("figures/top_trial.png")) file.remove("figures/top_trial.png")
