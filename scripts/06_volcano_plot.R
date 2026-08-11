# Volcano plot visualizing meta-analysis results (colorectal cancer)
# Load packages
library(tidyverse)
library(ggrepel)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

dir.create("results/figures/meta-analysis", showWarnings = FALSE, recursive = TRUE)

# Load Random Effect Model result
meta_result <- import("results/tables/meta-analysis/random_effect_model.csv") |>
  select(Gene_Symbol, randomSummary, randomP) |>
  rename(log2FC = randomSummary, P.Value = randomP)

# Significance thresholds: padj < 0.05 and |log2FC| > 1
meta_result <- meta_result |>
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC >  1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE                         ~ "NoSignificant"
  ))

# Top genes to label (annotated, significant, largest |log2FC|)
top_genes <- meta_result |>
  mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |>
  filter(!is.na(Gene_Symbol), P.Value < 0.05, abs(log2FC) > 1) |>
  slice_max(order_by = abs(log2FC), n = 500)

# Symmetric x-axis limit from the data (avoids clipping large fold changes)
x_lim   <- ceiling(max(abs(meta_result$log2FC), na.rm = TRUE))
x_break <- if (x_lim > 10) 2 else 1

# Volcano plot
volcano <- ggplot(meta_result, aes(x = log2FC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Down" = "#2c7fb8",
                                "NoSignificant" = "#636363",
                                "Up" = "#e34a33")) +
  theme_minimal() +
  labs(title = "", x = "log2 (Fold Change)", y = "-log10 (adjusted P-value)") +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "right",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(colour = "black", hjust = 1, size = 12),
    axis.text.y  = element_text(colour = "black", hjust = 1, size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  geom_text_repel(
    data = top_genes, aes(label = Gene_Symbol),
    colour = "black", size = 3, max.overlaps = 10,
    direction = "both", max.time = 5,
    force = 5, force_pull = 5, point.padding = 0.5, seed = 40
  ) +
  scale_x_continuous(limits = c(-x_lim, x_lim),
                     breaks = seq(-x_lim, x_lim, by = x_break))

# Save
ggsave("results/figures/meta-analysis/volcano_plot.png",
       plot = volcano, width = 12, height = 12, dpi = 600, bg = "white")

# ============================================================
# KIDNEY CANCER META-ANALYSIS
# Combining Method + Random-Effects Model Volcano Plots
# ============================================================

library(dplyr)
library(ggplot2)
library(rio)

# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_dir <- "D:/kidney_cancer_meta_analysis"

meta_dir <- file.path(
  project_dir,
  "results/tables/meta-analysis"
)

figure_dir <- file.path(
  project_dir,
  "results/figures/meta-analysis"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# PART A: COMBINING-METHOD META-VOLCANO
# ============================================================

# ------------------------------------------------------------
# 2. Load the combining-method result
# ------------------------------------------------------------

# Use the existing meta-analysis summary table
combining_file <- file.path(
  meta_dir,
  "meta_degs_summary_stats.csv"
)

combining_result <- rio::import(combining_file)

# Check columns
print(names(combining_result))


# ------------------------------------------------------------
# 3. Identify the correct combining-method columns
# ------------------------------------------------------------

# Your previous combining-method table contains:
# Gene_ID
# Gene_Symbol
# meta_log2FC
# meta_SE
# meta_pvalue
# meta_padj
# ci_lower
# ci_upper
# tau2
# I2
# Q
# Q_pvalue

combining_result <- combining_result %>%
  mutate(
    Regulation = case_when(
      meta_padj < 0.05 & meta_log2FC > 0 ~ "Upregulated",
      meta_padj < 0.05 & meta_log2FC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )


# ------------------------------------------------------------
# 4. Prepare volcano data
# ------------------------------------------------------------

combining_volcano_data <- combining_result %>%
  mutate(
    neg_log10_padj =
      -log10(pmax(meta_padj, .Machine$double.xmin)),
    
    Regulation = factor(
      Regulation,
      levels = c(
        "Downregulated",
        "Not significant",
        "Upregulated"
      )
    )
  )


# ------------------------------------------------------------
# 5. Combining-method volcano plot
# ------------------------------------------------------------

p_combining_volcano <- ggplot(
  combining_volcano_data,
  aes(
    x = meta_log2FC,
    y = neg_log10_padj
  )
) +
  
  geom_point(
    aes(shape = Regulation),
    size = 1.8,
    alpha = 0.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  labs(
    title = "Kidney Cancer Meta-analysis",
    subtitle = "Combining method",
    x = "Meta-analysis log2 Fold Change",
    y = "-log10 Adjusted P-value",
    shape = "Regulation"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.position = "right"
  )


# ------------------------------------------------------------
# 6. Save combining-method volcano
# ------------------------------------------------------------

ggsave(
  file.path(
    figure_dir,
    "Kidney_Cancer_Combining_Method_Meta_Volcano.png"
  ),
  plot = p_combining_volcano,
  width = 10,
  height = 8,
  dpi = 400
)

ggsave(
  file.path(
    figure_dir,
    "Kidney_Cancer_Combining_Method_Meta_Volcano.pdf"
  ),
  plot = p_combining_volcano,
  width = 10,
  height = 8,
  units = "in"
)


# ============================================================
# PART B: RANDOM-EFFECTS MODEL META-VOLCANO
# ============================================================

# ------------------------------------------------------------
# 7. Load corrected random-effects model
# ------------------------------------------------------------

random_file <- file.path(
  meta_dir,
  "random_effect_model_CORRECTED.csv"
)

random_result <- rio::import(random_file)

print(names(random_result))


# ------------------------------------------------------------
# 8. Create regulation classification
# ------------------------------------------------------------

random_result <- random_result %>%
  mutate(
    Regulation = case_when(
      meta_padj < 0.05 & meta_log2FC > 0 ~ "Upregulated",
      meta_padj < 0.05 & meta_log2FC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )


# ------------------------------------------------------------
# 9. Prepare random-effects volcano data
# ------------------------------------------------------------

random_volcano_data <- random_result %>%
  mutate(
    neg_log10_padj =
      -log10(pmax(meta_padj, .Machine$double.xmin)),
    
    Regulation = factor(
      Regulation,
      levels = c(
        "Downregulated",
        "Not significant",
        "Upregulated"
      )
    )
  )


# ------------------------------------------------------------
# 10. Random-effects model volcano
# ------------------------------------------------------------

p_random_volcano <- ggplot(
  random_volcano_data,
  aes(
    x = meta_log2FC,
    y = neg_log10_padj
  )
) +
  
  geom_point(
    aes(shape = Regulation),
    size = 1.8,
    alpha = 0.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  labs(
    title = "Kidney Cancer Meta-analysis",
    subtitle = "Random-effects model (REML)",
    x = "Meta-analysis log2 Fold Change",
    y = "-log10 Adjusted P-value",
    shape = "Regulation"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.position = "right"
  )


# ------------------------------------------------------------
# 11. Save random-effects volcano
# ------------------------------------------------------------

ggsave(
  file.path(
    figure_dir,
    "Kidney_Cancer_Random_Effects_Meta_Volcano.png"
  ),
  plot = p_random_volcano,
  width = 10,
  height = 8,
  dpi = 400
)

ggsave(
  file.path(
    figure_dir,
    "Kidney_Cancer_Random_Effects_Meta_Volcano.pdf"
  ),
  plot = p_random_volcano,
  width = 10,
  height = 8,
  units = "in"
)


# ============================================================
# PART C: CHECK RESULTS
# ============================================================

cat("\nCombining-method significant genes:\n")

print(
  sum(
    combining_result$meta_padj < 0.05,
    na.rm = TRUE
  )
)


cat("\nRandom-effects significant genes:\n")

print(
  sum(
    random_result$meta_padj < 0.05,
    na.rm = TRUE
  )
)


cat("\nCombining-method regulation:\n")

print(
  table(
    combining_result$Regulation
  )
)


cat("\nRandom-effects regulation:\n")

print(
  table(
    random_result$Regulation
  )
)


# ------------------------------------------------------------
# 12. Check output files
# ------------------------------------------------------------

list.files(
  figure_dir,
  pattern = "Volcano",
  full.names = TRUE
)
