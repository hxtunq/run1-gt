#!/usr/bin/env Rscript
#===============================================================================
# Variant Calling Benchmarking Visualization
# Tạo các biểu đồ so sánh hiệu suất của 4 variant callers
#===============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║        VARIANT CALLING BENCHMARK VISUALIZATION                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

#-------------------------------------------------------------------------------
# 1. Install and load packages
#-------------------------------------------------------------------------------
cat("[1/5] Loading packages...\n")

packages <- c("ggplot2", "dplyr", "tidyr", "readr", "scales", "patchwork", "gridExtra")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(gridExtra)
})

#-------------------------------------------------------------------------------
# 2. Configuration
#-------------------------------------------------------------------------------
cat("[2/5] Setting up configuration...\n")

args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) > 0, args[1], ".")
metrics_dir <- file.path(project_dir, "results", "final_metrics")
output_dir <- file.path(project_dir, "results", "figures")

cat(sprintf("  Project:  %s\n", project_dir))
cat(sprintf("  Output:   %s\n", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Color scheme
COLORS <- c(
  "gatk" = "#E41A1C",
  "deepvariant" = "#377EB8",
  "strelka2" = "#4DAF4A",
  "freebayes" = "#984EA3"
)

LABELS <- c(
  "gatk" = "GATK",
  "deepvariant" = "DeepVariant",
  "strelka2" = "Strelka2",
  "freebayes" = "FreeBayes"
)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot. title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95")
    )
}

#-------------------------------------------------------------------------------
# 3. Load data
#-------------------------------------------------------------------------------
cat("[3/5] Loading data...\n")

bench_file <- file.path(metrics_dir, "benchmark_summary. tsv")

if (file.exists(bench_file)) {
  data <- read_tsv(bench_file, show_col_types = FALSE)
  cat(sprintf("  Loaded %d records from %s\n", nrow(data), bench_file))
} else {
  cat("  [WARN] benchmark_summary.tsv not found, using demo data\n")
  
  # Demo data
  data <- tibble(
    Caller = rep(c("gatk", "deepvariant", "strelka2", "freebayes"), each = 3),
    VariantType = rep(c("ALL", "SNP", "INDEL"), 4),
    TP = c(47500, 42500, 5000, 48500, 43500, 5000, 47000, 42000, 5000, 46000, 41000, 5000),
    FP = c(250, 150, 100, 100, 50, 50, 300, 200, 100, 500, 350, 150),
    FN = c(2500, 2250, 250, 1500, 1250, 250, 3000, 2750, 250, 4000, 3750, 250),
    Precision = c(0.9948, 0.9965, 0.9804, 0.9979, 0.9989, 0.9901, 
                  0.9937, 0.9953, 0.9804, 0.9892, 0.9915, 0.9709),
    Recall = c(0.9500, 0.9497, 0.9524, 0.9700, 0.9720, 0.9524,
               0.9400, 0.9387, 0.9524, 0.9200, 0.9163, 0.9524),
    F1 = c(0.9719, 0.9726, 0.9662, 0.9838, 0.9852, 0.9709,
           0.9662, 0.9663, 0.9662, 0.9534, 0.9527, 0.9615)
  )
}

# Standardize
data <- data %>%
  rename_with(tolower) %>%
  mutate(
    caller = tolower(caller),
    varianttype = toupper(varianttype)
  )

#-------------------------------------------------------------------------------
# 4. Create plots
#-------------------------------------------------------------------------------
cat("[4/5] Creating plots...\n")

# Helper function
save_plot <- function(plot, name, width = 10, height = 6) {
  png_file <- file.path(output_dir, paste0(name, ". png"))
  pdf_file <- file.path(output_dir, paste0(name, ".pdf"))
  
  ggsave(png_file, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height)
  
  cat(sprintf("  Saved:  %s\n", basename(png_file)))
}

#--- Plot 1: F1 Score Comparison ---#
cat("  Creating F1 score plot...\n")

p1 <- data %>%
  filter(varianttype %in% c("ALL", "SNP", "INDEL")) %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    varianttype = factor(varianttype, levels = c("ALL", "SNP", "INDEL"))
  ) %>%
  ggplot(aes(x = caller_label, y = f1, fill = caller)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", f1)), vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ varianttype, nrow = 1) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2), labels = percent) +
  labs(
    title = "F1 Score Comparison Across Variant Callers",
    subtitle = "Stratified by Variant Type (ALL, SNP, INDEL)",
    x = NULL,
    y = "F1 Score"
  ) +
  theme_pub()

save_plot(p1, "01_f1_score_comparison", width = 12, height = 5)

#--- Plot 2: Precision vs Recall ---#
cat("  Creating precision-recall plot...\n")

p2 <- data %>%
  filter(varianttype %in% c("SNP", "INDEL")) %>%
  mutate(caller_label = factor(caller, levels = names(COLORS), labels = LABELS)) %>%
  ggplot(aes(x = recall, y = precision, color = caller, shape = varianttype)) +
  geom_point(size = 6, stroke = 1.2) +
  scale_color_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_shape_manual(values = c("SNP" = 16, "INDEL" = 17), name = "Variant Type") +
  scale_x_continuous(limits = c(0. 90, 1.0), labels = percent, breaks = seq(0.90, 1, 0.02)) +
  scale_y_continuous(limits = c(0.96, 1.0), labels = percent, breaks = seq(0.96, 1, 0.01)) +
  labs(
    title = "Precision vs Recall Trade-off",
    subtitle = "Ideal performance:  top-right corner",
    x = "Recall (Sensitivity)",
    y = "Precision (PPV)"
  ) +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

save_plot(p2, "02_precision_recall", width = 9, height = 6)

#--- Plot 3: TP/FP/FN Stacked Bar ---#
cat("  Creating classification plot...\n")

class_colors <- c(
  "True Positive" = "#27AE60",
  "False Positive" = "#E74C3C",
  "False Negative" = "#F39C12"
)

p3 <- data %>%
  filter(varianttype == "ALL") %>%
  select(caller, tp, fp, fn) %>%
  pivot_longer(c(tp, fp, fn), names_to = "category", values_to = "count") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    category = factor(category, 
                      levels = c("tp", "fp", "fn"),
                      labels = c("True Positive", "False Positive", "False Negative"))
  ) %>%
  ggplot(aes(x = caller_label, y = count, fill = category)) +
  geom_col(position = "stack", color = "black", width = 0.7, linewidth = 0.3) +
  geom_text(
    aes(label = comma(count)),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = class_colors, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Variant Classification by Caller",
    subtitle = "Distribution of TP, FP, and FN (All Variants)",
    x = NULL,
    y = "Number of Variants"
  ) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(p3, "03_variant_classification", width = 9, height = 6)

#--- Plot 4: Metrics Side-by-Side ---#
cat("  Creating metrics comparison plot...\n")

p4 <- data %>%
  filter(varianttype == "ALL") %>%
  select(caller, precision, recall, f1) %>%
  pivot_longer(c(precision, recall, f1), names_to = "metric", values_to = "value") %>%
  mutate(
    caller_label = factor(caller, levels = names(COLORS), labels = LABELS),
    metric = factor(metric, 
                    levels = c("precision", "recall", "f1"),
                    labels = c("Precision", "Recall", "F1 Score"))
  ) %>%
  ggplot(aes(x = metric, y = value, fill = caller)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2), labels = percent) +
  labs(
    title = "Performance Metrics Comparison",
    subtitle = "Precision, Recall, and F1 Score (All Variants)",
    x = NULL,
    y = "Value"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_plot(p4, "04_metrics_comparison", width = 10, height = 6)

#--- Plot 5: SNP vs INDEL ---#
cat("  Creating SNP vs INDEL comparison.. .\n")

p5 <- data %>%
  filter(varianttype %in% c("SNP", "INDEL")) %>%
  select(caller, varianttype, f1) %>%
  pivot_wider(names_from = varianttype, values_from = f1) %>%
  mutate(caller_label = factor(caller, levels = names(COLORS), labels = LABELS)) %>%
  ggplot(aes(x = SNP, y = INDEL, color = caller)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 8, alpha = 0.8) +
  geom_text(aes(label = caller_label), vjust = -1. 2, size = 4, fontface = "bold",
            show.legend = FALSE) +
  scale_color_manual(values = COLORS, labels = LABELS, name = "Caller") +
  scale_x_continuous(limits = c(0.94, 1.0), labels = percent) +
  scale_y_continuous(limits = c(0.94, 1.0), labels = percent) +
  labs(
    title = "SNP vs INDEL Calling Performance",
    subtitle = "F1 Score (diagonal = equal performance)",
    x = "F1 Score (SNPs)",
    y = "F1 Score (INDELs)"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  coord_fixed()

save_plot(p5, "05_snp_vs_indel", width = 8, height = 8)

#--- Plot 6: Summary Dashboard ---#
cat("  Creating summary dashboard...\n")

dashboard <- (p1 / (p2 | p3)) +
  plot_annotation(
    title = "Variant Calling Benchmarking Dashboard",
    subtitle = paste("GATK | DeepVariant | Strelka2 | FreeBayes —", format(Sys.Date(), "%B %d, %Y")),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

save_plot(dashboard, "06_summary_dashboard", width = 14, height = 12)

#-------------------------------------------------------------------------------
# 5. Generate summary table
#-------------------------------------------------------------------------------
cat("[5/5] Generating summary table...\n")

summary_table <- data %>%
  filter(varianttype == "ALL") %>%
  mutate(
    Caller = factor(caller, levels = names(COLORS), labels = LABELS)
  ) %>%
  select(
    Caller,
    TP = tp,
    FP = fp,
    FN = fn,
    Precision = precision,
    Recall = recall,
    F1 = f1
  ) %>%
  arrange(desc(F1))

# Save as CSV
write_csv(summary_table, file.path(output_dir, "benchmark_summary_table.csv"))
cat(sprintf("  Saved:  benchmark_summary_table.csv\n"))

#-------------------------------------------------------------------------------
# Print summary
#-------------------------------------------------------------------------------
cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("BENCHMARK SUMMARY (Ranked by F1 Score)\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# Print formatted table
print(
  summary_table %>%
    mutate(
      Precision = sprintf("%.4f", Precision),
      Recall = sprintf("%.4f", Recall),
      F1 = sprintf("%.4f", as.numeric(F1))
    ),
  n = Inf
)

cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║              VISUALIZATION COMPLETE!                            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nGenerated files:\n")
for (f in list.files(output_dir, pattern = "\\.(png|pdf|csv)$")) {
  cat(sprintf("  • %s\n", f))
}
cat("\n")
