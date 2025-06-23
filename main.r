# Quality Control Analysis of
# Laser Welded Steel-copper Lap Joints
#
#
source("utils.r")
get_requirements("./requirements.txt")

# Load all required libraries
library(qcc)
library(dplyr)


# Create analysis directory if it doesn't exist
analysis_dir <- "analysis"
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir)
  cat("Created analysis directory:", analysis_dir, "\n")
}

data_file <- "data.csv"

if (!file.exists(data_file)) {
  stop("Data file 'data.csv' not found.")
}
data <- read.csv(data_file, sep = ",", dec = ".", stringsAsFactors = FALSE)
cat("=== COLUMN NAMES IN DATA ===\n")
print(colnames(data))
cat("\n")

#### Processing and cleaning the data
#
#


# Create binary cracking variable (1 = yes, 0 = no)
data$cracking_binary <- ifelse(
  data$cracking.in.the.weld.metal..yes.no. == "yes", 1, 0
)

numeric_cols <- c(
  "power..W.", "welding.speed..m.min.", "gas.flow.rate..l.min.",
  "focal.position..mm.", "angular.position....",
  "material.strenght.steel..mm.", "weld.depth.copper..µm.",
  "count.of.cracks"
)

for (col in numeric_cols) {
  data[[col]] <- as.numeric(data[[col]])
  data[[col]][is.na(data[[col]])] <- 0
}

df <- data[data$weld.depth.copper..µm. != 0, ]

cat("=== DATA SUMMARY ===\n")
cat("Original data rows:", nrow(data), "\n")
cat("Clean data rows:", nrow(df), "\n")
cat("Cracking cases:", sum(df$cracking_binary), "\n")
cat("No cracking cases:", sum(df$cracking_binary == 0), "\n\n")

#
#
###

### Proportion of Cracking Cases
#
#
cat("=== PROPORTION OF CRACKING CASES ===\n")

# Grouping by weldnumber ensures we summarize correctly per group of 4 cuts × 5 repeats
grouped <- df |>
  group_by(weldnumber) |>
  summarise(
    n_inspected = n(), # total number of cross-sections per weld
    cracked = sum(cracking_binary), # number of cross-sections that cracked
    crack_counts = sum(count.of.cracks, na.rm = TRUE), # total number of cracks (not just cracked yes/no)
    .groups = "drop"
  ) |>
  arrange(weldnumber)

# Extracting vectors to use in qcc()
n_inspected <- grouped$n_inspected
cracked <- grouped$cracked
crack_counts <- grouped$crack_counts

# Logging to console
cat("Total welds inspected:", sum(n_inspected), "\n")
cat("Total welds with cracks:", sum(cracked), "\n")
cat("Total cracks:", sum(crack_counts), "\n")

# Plotting the p-chart and c-chart
png(file.path(analysis_dir, "proportion_of_cracking_cases.png"), width = 800, height = 600)
p_cracking_cases <- qcc(cracked, sizes = n_inspected, type = "p", title = "p Chart: Cracking Proportion")
dev.off()
summary(p_cracking_cases)

png(file.path(analysis_dir, "qty_of_cracking_cases.png"), width = 800, height = 600)
np_cracking_cases <- qcc(cracked, sizes = n_inspected, type = "np", title = "np Chart: Cracking Count")
dev.off()
summary(np_cracking_cases)

#
#
###

### Logistic Regression Analysis
#
#
cat("=== LOGISTIC REGRESSION ANALYSIS ===\n")
logistic_model <- glm(
  cracking_binary ~
    power..W. +
    welding.speed..m.min. +
    gas.flow.rate..l.min. +
    focal.position..mm. +
    angular.position.... +
    material.strenght.steel..mm. +
    weld.depth.copper..µm.,
  data = df,
  family = binomial
)

print(summary(logistic_model))

df$predicted_prob <- predict(logistic_model, type = "response")

png(file.path(analysis_dir, "crack_prob_vs_weld_depth.png"), width = 800, height = 600)
plot(df$weld.depth.copper..µm., df$predicted_prob,
  xlab = "Weld Depth (µm)", ylab = "Cracking Probability",
  main = "Cracking Risk vs Weld Depth",
  pch = 16, col = ifelse(df$cracking_binary == 1, "red", "blue")
)
abline(h = 0.5, col = "orange", lty = 2, lwd = 2)
legend("bottomright", c("Cracked", "Not Cracked", "50% Risk Line"),
  col = c("red", "blue", "orange"), pch = c(16, 16, NA), lty = c(NA, NA, 2)
)
dev.off()
#
#
###

### Monitoring Weld Depth Copper
#
#
grouped_depth <- df |>
  group_by(weldnumber) |>
  summarise(
    values = list(weld.depth.copper..µm.),
    .groups = "drop"
  )

depth_matrix <- do.call(rbind, lapply(grouped_depth$values, function(x) as.numeric(x)))

png(file.path(analysis_dir, "weld_depth_copper.png"), width = 800, height = 600)
x_chart_weld_depth_copper <- qcc(depth_matrix, type = "xbar", title = "x̄ Chart: Weld Depth in Copper")
dev.off()

summary(x_chart_weld_depth_copper)

png(file.path(analysis_dir, "weld_depth_copper_individual.png"), width = 800, height = 600)
qcc_depth <- qcc(df$weld.depth.copper..µm.,
  type = "xbar.one",
  plot = TRUE,
  title = "Individual Control Chart: Weld Depth"
)
dev.off()
summary(qcc_depth)
#
#
###

### Iterative Outlier Removal (Smaller is Better)
#
#
cat("=== ITERATIVE OUTLIER REMOVAL (SMALLER IS BETTER) ===\n")

df_conforme <- df
max_iterations <- 10
tolerance <- 1e-6

for (i in 1:max_iterations) {
  # Create control chart
  qcc_depth <- qcc(df_conforme$weld.depth.copper..µm.,
    type = "xbar.one",
    plot = FALSE, # Don't plot during iterations
    title = paste("Iteration", i)
  )

  outside_upper_limit <- sum(df_conforme$weld.depth.copper..µm. > qcc_depth$limits[2])

  cat(
    "Iteration", i, "- Points above upper limit:", outside_upper_limit,
    "- Sample size:", nrow(df_conforme), "\n"
  )

  if (outside_upper_limit == 0) {
    cat("Convergence achieved at iteration", i, "- No points above upper control limit\n")
    break
  }

  df_conforme_new <- df_conforme %>%
    filter(weld.depth.copper..µm. <= qcc_depth$limits[2])

  if (abs(nrow(df_conforme_new) - nrow(df_conforme)) < tolerance) {
    cat("No significant change in sample size. Stopping.\n")
    break
  }

  df_conforme <- df_conforme_new
}

# Final control chart with plot
png(file.path(analysis_dir, "weld_depth_final_control.png"), width = 800, height = 600)
qcc_depth_final <- qcc(df_conforme$weld.depth.copper..µm.,
  type = "xbar.one",
  plot = TRUE,
  title = "Individual Control Chart - Weld Depth (After Outlier Removal)"
)
dev.off()

# Process capability analysis
png(file.path(analysis_dir, "process_capability.png"), width = 800, height = 600)
capability_analysis <- process.capability(qcc_depth_final, spec.limits = c(0, 50))
dev.off()

cat("\n=== FINAL RESULTS ===\n")
cat("Original sample size:", nrow(df), "\n")
cat("Final conforming sample size:", nrow(df_conforme), "\n")
cat("Data reduction:", round((1 - nrow(df_conforme) / nrow(df)) * 100, 2), "%\n")

summary(qcc_depth_final)
summary(capability_analysis)
