# Summary of meta-analysis results (Kidney cancer,tumor vs normal)
# Load packages
library(tidyverse)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

meta_dir <- "results/tables/meta-analysis"

# Up / down regulation summary from the Random Effect Model result
meta_result <- import(file.path(meta_dir, "random_effect_model.csv")) |>
  select(Gene_Symbol, randomSummary, randomP, Gene_Description) |>
  rename(log2FC = randomSummary, P.Value = randomP)

# Significance call: padj < 0.05 and |log2FC| > 1
meta_result <- meta_result |>
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC >  1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE                         ~ "NS"
  ))

# Counts and percentage of up- / down-regulated genes
gene_stats <- meta_result |>
  filter(Significance != "NS") |>
  group_by(Significance) |>
  summarise(
    Count      = n(),
    Percentage = n() / nrow(meta_result) * 100,
    .groups    = "drop"
  )

print(gene_stats)
export(gene_stats, file.path(meta_dir, "meta_degs_summary_stats.csv"))

# Significant DEGs that carry a gene symbol (annotated only)
annotated_genes <- meta_result |>
  mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |>
  filter(!is.na(Gene_Symbol), Significance != "NS")

export(annotated_genes, file.path(meta_dir, "filtered_meta_degs_annotated_only.csv"))

# Regulation table from the filtered meta DEGs
meta_key_results <- import(file.path(meta_dir, "filtered_meta_degs.csv"), na.strings = "") |>
  select(Gene_ID, Gene_Symbol, Gene_Description,
         log2FoldChange = randomSummary, adjusted.P.Value = randomP) |>
  drop_na(Gene_Symbol) |>
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

export(meta_key_results, file.path(meta_dir, "meta_degs_regulation.csv"))

# ------------------------------------------------------------------------------
# 04. Cross-reference Venn intersection genes with meta-analysis effect sizes
venn_common <- "results/tables/venn/common_genes_all_datasets.csv"
comb_file   <- file.path(meta_dir, "meta_combining_mean.csv")

if (file.exists(venn_common) && file.exists(comb_file)) {
  intersect_ids <- import(venn_common) |> pull(Gene_ID) |> as.character()
  
  combined_meta <- import(comb_file) |>
    mutate(Gene_ID = as.character(Gene_ID))
  
  intersects <- combined_meta |> filter(Gene_ID %in% intersect_ids)
  
  export(intersects, "results/tables/venn/intersect_genes_meta_analysis.csv")
  message("Venn-intersect genes with meta stats: ", nrow(intersects))
} else {
  message("Missing Venn or combining-approach file; skipping cross-reference. Run 03_venn.R and 04_metavolcano.R first.")
}

meta_result <- rio::import(
  file.path(meta_dir, "random_effect_model.csv")
)


head(meta_result)
meta_result <- meta_result %>%
  dplyr::select(
    Gene_ID,
    Gene_Symbol,
    meta_log2FC,
    meta_pvalue,
    meta_padj,
    ci_lower,
    ci_upper,
    tau2,
    I2,
    Q,
    Q_pvalue
  ) %>%
  dplyr::rename(
    log2FC = meta_log2FC,
    P.Value = meta_pvalue,
    FDR = meta_padj
  )
head(meta_result)
names(meta_result)
meta_result <- meta_result %>%
  mutate(
    Regulation = case_when(
      FDR < 0.05 & log2FC > 0 ~ "Upregulated",
      FDR < 0.05 & log2FC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

table(meta_result$Regulation)
meta_result <- rio::import(
  file.path(meta_dir, "random_effect_model.csv")
)

names(meta_result)

head(meta_result)

library(metafor)

# Check Gene_ID 2 manually
test_gene <- meta_input[
  meta_input$Gene_ID == 2,
]

test_gene
exists("gse24455")
exists("gse63183")
ls()
list.files(
  "D:/kidney_cancer_meta_analysis",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
meta_result <- rio::import(
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/random_effect_model.csv"
)

dim(meta_result)
names(meta_result)

table(
  ifelse(
    meta_result$meta_padj < 0.05 & meta_result$meta_log2FC > 0,
    "Upregulated",
    ifelse(
      meta_result$meta_padj < 0.05 & meta_result$meta_log2FC < 0,
      "Downregulated",
      "Not significant"
    )
  )
)

# Load the two original DESeq2 result files
gse24455 <- rio::import(
  "D:/kidney_cancer_meta_analysis/results/tables/DESeq2/GSE24455.csv"
)

gse63183 <- rio::import(
  "D:/kidney_cancer_meta_analysis/results/tables/DESeq2/GSE63183.csv"
)

# Check their column names
names(gse24455)
names(gse63183)
# Check Gene ID 2 in both datasets
gse24455[gse24455$Gene_ID == 2, ]
gse63183[gse63183$Gene_ID == 2, ]
meta_result[meta_result$Gene_ID == 2, ]
library(dplyr)
library(metafor)

s1 <- gse24455 %>%
  select(Gene_ID, log2FoldChange, lfcSE, pvalue, padj) %>%
  rename(
    logFC_24455 = log2FoldChange,
    SE_24455 = lfcSE,
    p_24455 = pvalue,
    padj_24455 = padj
  )

s2 <- gse63183 %>%
  select(Gene_ID, log2FoldChange, lfcSE, pvalue, padj) %>%
  rename(
    logFC_63183 = log2FoldChange,
    SE_63183 = lfcSE,
    p_63183 = pvalue,
    padj_63183 = padj
  )

meta_input <- inner_join(
  s1,
  s2,
  by = "Gene_ID"
)

dim(meta_input)
head(meta_input)
test_manual <- metafor::rma.uni(
  yi = c(
    meta_input$logFC_24455[meta_input$Gene_ID == 2],
    meta_input$logFC_63183[meta_input$Gene_ID == 2]
  ),
  sei = c(
    meta_input$SE_24455[meta_input$Gene_ID == 2],
    meta_input$SE_63183[meta_input$Gene_ID == 2]
  ),
  method = "REML",
  test = "z"
)

test_manual

run_meta_REML <- function(fc1, se1, fc2, se2) {
  
  fit <- tryCatch(
    metafor::rma.uni(
      yi = c(fc1, fc2),
      sei = c(se1, se2),
      method = "REML",
      test = "z"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(c(
      meta_log2FC = NA_real_,
      meta_SE = NA_real_,
      meta_pvalue = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      Q = NA_real_,
      Q_pvalue = NA_real_
    ))
  }
  
  c(
    meta_log2FC = as.numeric(coef(fit)),
    meta_SE = fit$se,
    meta_pvalue = fit$pval,
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
    tau2 = fit$tau2,
    I2 = fit$I2,
    Q = fit$QE,
    Q_pvalue = fit$QEp
  )
}

reml_results <- t(
  vapply(
    seq_len(nrow(meta_input)),
    function(i) {
      run_meta_REML(
        meta_input$logFC_24455[i],
        meta_input$SE_24455[i],
        meta_input$logFC_63183[i],
        meta_input$SE_63183[i]
      )
    },
    FUN.VALUE = c(
      meta_log2FC = 0,
      meta_SE = 0,
      meta_pvalue = 0,
      ci_lower = 0,
      ci_upper = 0,
      tau2 = 0,
      I2 = 0,
      Q = 0,
      Q_pvalue = 0
    )
  )
)

reml_results <- as.data.frame(reml_results)

reml_results$Gene_ID <- meta_input$Gene_ID

reml_results <- reml_results[, c(
  "Gene_ID",
  "meta_log2FC",
  "meta_SE",
  "meta_pvalue",
  "ci_lower",
  "ci_upper",
  "tau2",
  "I2",
  "Q",
  "Q_pvalue"
)]

reml_results$meta_padj <- p.adjust(
  reml_results$meta_pvalue,
  method = "BH"
)

reml_results[reml_results$Gene_ID == 2, ]
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(reml_results$Gene_ID),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL")
)

gene_annotation <- gene_annotation[
  !duplicated(gene_annotation$ENTREZID),
  c("ENTREZID", "SYMBOL")
]

colnames(gene_annotation) <- c(
  "Gene_ID",
  "Gene_Symbol"
)

reml_results_annotated <- merge(
  reml_results,
  gene_annotation,
  by = "Gene_ID",
  all.x = TRUE,
  sort = FALSE
)

dir.create(
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  reml_results_annotated,
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/random_effect_model_CORRECTED.csv",
  row.names = FALSE
)

sum(
  reml_results_annotated$meta_padj < 0.05,
  na.rm = TRUE
)

table(
  ifelse(
    reml_results_annotated$meta_padj < 0.05 &
      reml_results_annotated$meta_log2FC > 0,
    "Upregulated",
    ifelse(
      reml_results_annotated$meta_padj < 0.05 &
        reml_results_annotated$meta_log2FC < 0,
      "Downregulated",
      "Not significant"
    )
  )
)
# Corrected meta-analysis result
corrected_meta <- reml_results_annotated

# Significant meta-DEGs
sig_meta_corrected <- corrected_meta %>%
  filter(!is.na(meta_padj), meta_padj < 0.05)

# Upregulated
up_meta_corrected <- sig_meta_corrected %>%
  filter(meta_log2FC > 0)

# Downregulated
down_meta_corrected <- sig_meta_corrected %>%
  filter(meta_log2FC < 0)

# Save
write.csv(
  sig_meta_corrected,
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/filtered_meta_degs_CORRECTED.csv",
  row.names = FALSE
)

write.csv(
  up_meta_corrected,
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/upregulated_meta_degs_CORRECTED.csv",
  row.names = FALSE
)

write.csv(
  down_meta_corrected,
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/downregulated_meta_degs_CORRECTED.csv",
  row.names = FALSE
)

# Check
dim(sig_meta_corrected)
dim(up_meta_corrected)
dim(down_meta_corrected)
meta_summary <- corrected_meta %>%
  mutate(
    Regulation = case_when(
      meta_padj < 0.05 & meta_log2FC > 0 ~ "Upregulated",
      meta_padj < 0.05 & meta_log2FC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  select(
    Gene_ID,
    Gene_Symbol,
    meta_log2FC,
    meta_SE,
    meta_pvalue,
    meta_padj,
    ci_lower,
    ci_upper,
    tau2,
    I2,
    Q,
    Q_pvalue,
    Regulation
  )

write.csv(
  meta_summary,
  "D:/kidney_cancer_meta_analysis/results/tables/meta-analysis/final_meta_analysis_results.csv",
  row.names = FALSE
)

library(ggplot2)
library(ggrepel)
library(dplyr)

volcano_data <- meta_summary %>%
  mutate(
    neg_log10_padj = -log10(pmax(meta_padj, .Machine$double.xmin)),
    Regulation = factor(
      Regulation,
      levels = c(
        "Downregulated",
        "Not significant",
        "Upregulated"
      )
    )
  )

p_volcano_corrected <- ggplot(
  volcano_data,
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
    subtitle = "Random-effects model",
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
ggsave(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Volcano_CORRECTED.png",
  plot = p_volcano_corrected,
  width = 10,
  height = 8,
  dpi = 400
)
ggsave(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Volcano_CORRECTED.pdf",
  plot = p_volcano_corrected,
  width = 10,
  height = 8,
  units = "in"
)
file.exists(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Volcano_CORRECTED.png"
)
library(dplyr)
library(ggplot2)

forest_genes_corrected <- corrected_meta %>%
  filter(
    !is.na(meta_padj),
    meta_padj < 0.05
  ) %>%
  arrange(meta_padj) %>%
  slice_head(n = 20) %>%
  mutate(
    Gene_Label = paste0(
      Gene_Symbol,
      " (",
      Gene_ID,
      ")"
    )
  )

forest_genes_corrected <- forest_genes_corrected %>%
  arrange(meta_log2FC) %>%
  mutate(
    Gene_Label = factor(
      Gene_Label,
      levels = Gene_Label
    )
  )

p_forest_corrected <- ggplot(
  forest_genes_corrected,
  aes(
    x = meta_log2FC,
    y = Gene_Label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    width = 0.15
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "Kidney Cancer Meta-analysis",
    subtitle = "Top 20 significant genes — Random-effects model",
    x = "Meta-analysis log2 Fold Change (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

p_forest_corrected
ggsave(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Forest_CORRECTED.png",
  plot = p_forest_corrected,
  width = 8,
  height = 10,
  dpi = 400
)
ggsave(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Forest_CORRECTED.pdf",
  plot = p_forest_corrected,
  width = 8,
  height = 10,
  units = "in"
)
file.exists(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Forest_CORRECTED.png"
)

file.exists(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_Meta_Forest_CORRECTED.pdf"
)
corrected_meta

library(dplyr)
library(VennDiagram)
library(grid)
library(rio)

gse24455 <- rio::import(
  "D:/kidney_cancer_meta_analysis/results/tables/DESeq2/GSE24455.csv"
)

gse63183 <- rio::import(
  "D:/kidney_cancer_meta_analysis/results/tables/DESeq2/GSE63183.csv"
)

deg_24455 <- gse24455 %>%
  filter(!is.na(padj), padj < 0.05)

deg_63183 <- gse63183 %>%
  filter(!is.na(padj), padj < 0.05)

length(unique(deg_24455$Gene_ID))
length(unique(deg_63183$Gene_ID))
venn_24455 <- unique(deg_24455$Gene_ID)
venn_63183 <- unique(deg_63183$Gene_ID)

venn_common <- intersect(
  venn_24455,
  venn_63183
)

venn_24455_only <- setdiff(
  venn_24455,
  venn_63183
)

venn_63183_only <- setdiff(
  venn_63183,
  venn_24455
)

length(venn_24455_only)
length(venn_63183_only)
length(venn_common)
deg_overlap_summary <- data.frame(
  Category = c(
    "GSE24455 significant DEGs",
    "GSE63183 significant DEGs",
    "Common DEGs",
    "GSE24455-specific DEGs",
    "GSE63183-specific DEGs"
  ),
  Number_of_DEGs = c(
    length(venn_24455),
    length(venn_63183),
    length(venn_common),
    length(venn_24455_only),
    length(venn_63183_only)
  )
)

deg_overlap_summary

write.csv(
  deg_overlap_summary,
  "D:/kidney_cancer_meta_analysis/results/tables/venn/DEG_overlap_summary.csv",
  row.names = FALSE
)
common_deg_table <- gse24455 %>%
  filter(Gene_ID %in% venn_common) %>%
  select(
    Gene_ID,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  ) %>%
  rename(
    log2FC_24455 = log2FoldChange,
    SE_24455 = lfcSE,
    p_24455 = pvalue,
    padj_24455 = padj
  ) %>%
  inner_join(
    gse63183 %>%
      filter(Gene_ID %in% venn_common) %>%
      select(
        Gene_ID,
        log2FoldChange,
        lfcSE,
        pvalue,
        padj
      ) %>%
      rename(
        log2FC_63183 = log2FoldChange,
        SE_63183 = lfcSE,
        p_63183 = pvalue,
        padj_63183 = padj
      ),
    by = "Gene_ID"
  )

dim(common_deg_table)
head(common_deg_table)
write.csv(
  common_deg_table,
  "D:/kidney_cancer_meta_analysis/results/tables/venn/common_DEGs_904.csv",
  row.names = FALSE
)
deg_frequency <- data.frame(
  Gene_ID = unique(c(venn_24455, venn_63183))
) %>%
  mutate(
    GSE24455 = ifelse(Gene_ID %in% venn_24455, "Yes", "No"),
    GSE63183 = ifelse(Gene_ID %in% venn_63183, "Yes", "No"),
    Dataset_Count =
      as.numeric(GSE24455 == "Yes") +
      as.numeric(GSE63183 == "Yes")
  ) %>%
  arrange(desc(Dataset_Count), Gene_ID)
write.csv(
  deg_frequency,
  "D:/kidney_cancer_meta_analysis/results/tables/venn/deg_frequency.csv",
  row.names = FALSE
)
venn_plot <- venn.diagram(
  x = list(
    GSE24455 = unique(venn_24455),
    GSE63183 = unique(venn_63183)
  ),
  filename = NULL,
  fill = c("grey70", "grey40"),
  alpha = 0.55,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  margin = 0.1,
  main = "Kidney Cancer DEG Overlap",
  main.cex = 1.5
)
png(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis/Kidney_Cancer_DEG_Overlap_Venn.png",
  width = 2400,
  height = 2200,
  res = 400
)

grid.newpage()
grid.draw(venn_plot)

dev.off()

list.files(
  "D:/kidney_cancer_meta_analysis/results/tables/venn",
  full.names = TRUE
)

list.files(
  "D:/kidney_cancer_meta_analysis/results/figures/meta-analysis",
  pattern = "Venn|Overlap",
  full.names = TRUE
)
