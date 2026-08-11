# Functional enrichment (ORA + GSEA) of colorectal cancer meta-analysis DEGs
# Load packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
library(rio)

# org.Hs.eg.db (AnnotationDbi) masks dplyr verbs -- pin them back
select <- dplyr::select
filter <- dplyr::filter

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter

# Thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# Input / output
META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"
OUT_CSV   <- "results/tables/enrichment"
OUT_FIG   <- "results/figures/enrichment"
PREFIX    <- "colorectal"
dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# Helper functions
# Ranked list for GSEA. Names = Entrez (Gene_ID), score = sign(log2FC)*-log10(p).
build_ranked_list <- function(meta_df) {
  df <- meta_df |>
    filter(!is.na(randomP), !is.na(randomSummary), randomP > 0, !is.na(Gene_ID)) |>
    mutate(rank_score = sign(randomSummary) * -log10(randomP),
           Gene_ID    = as.character(Gene_ID)) |>
    filter(!duplicated(Gene_ID)) |>
    arrange(desc(rank_score))
  setNames(df$rank_score, df$Gene_ID)
}

# ORA: GO (BP/MF/CC) + KEGG + Reactome
run_ora <- function(entrez_ids, universe_entrez) {
  if (length(entrez_ids) < 5) {
    message("    Too few genes for ORA (n=", length(entrez_ids), "), skipping.")
    return(NULL)
  }
  results <- list()
  
  for (ont in c("BP", "MF", "CC")) {
    results[[paste0("GO_", ont)]] <- tryCatch(
      enrichGO(gene = entrez_ids, universe = universe_entrez, OrgDb = org.Hs.eg.db,
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { message("    GO-", ont, " error: ", e$message); NULL }
    )
  }
  
  results[["KEGG"]] <- tryCatch(
    enrichKEGG(gene = entrez_ids, universe = universe_entrez, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) { message("    KEGG error: ", e$message); NULL }
  )
  
  results[["Reactome"]] <- tryCatch(
    enrichPathway(gene = entrez_ids, universe = universe_entrez, organism = "human",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) { message("    Reactome error: ", e$message); NULL }
  )
  
  results
}

# GSEA: GO-BP + KEGG + Hallmark
run_gsea <- function(ranked_list) {
  if (length(ranked_list) < 100) {
    message("    Too few ranked genes (n=", length(ranked_list), "), skipping GSEA.")
    return(NULL)
  }
  results <- list()
  
  results[["GSEA_GO_BP"]] <- tryCatch(
    gseGO(geneList = ranked_list, OrgDb = org.Hs.eg.db, ont = "BP",
          minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA GO error: ", e$message); NULL }
  )
  
  results[["GSEA_KEGG"]] <- tryCatch(
    gseKEGG(geneList = ranked_list, organism = "hsa",
            minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA KEGG error: ", e$message); NULL }
  )
  
  # msigdbr >= 10: use collection= (category deprecated); Entrez col is ncbi_gene
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
    select(gs_name, ncbi_gene) |>
    mutate(ncbi_gene = as.character(ncbi_gene))
  
  results[["GSEA_Hallmark"]] <- tryCatch(
    GSEA(geneList = ranked_list, TERM2GENE = hallmark,
         minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA Hallmark error: ", e$message); NULL }
  )
  
  results
}

# Export enrichment results to CSV
export_enrichment <- function(result_list, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (is.null(res)) next
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    export(df, file.path(out_dir, paste0(prefix, "_", name, ".csv")))
  }
}

# Tidy long MSigDB / GO set names for plotting.
.pretty_set_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x); x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  acronyms <- c("DNA","RNA","TNFA","NFKB","IL1","IL2","IL6","JAK","STAT1","STAT3",
                "STAT5","KRAS","MTORC1","TGF","E2F","UV","WNT","MYC")
  vapply(strsplit(x, " ", fixed = TRUE), function(words) {
    hit <- toupper(words) %in% acronyms
    words[hit] <- toupper(words[hit])
    paste(words, collapse = " ")
  }, character(1))
}

# Dotplot. GSEA drawn on NES axis; ORA uses clusterProfiler dotplot.
save_dotplot <- function(res, title, path, width = 10, height = 8, show = 20) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  is_gsea <- all(c("NES", "setSize") %in% colnames(df))
  if (is_gsea) {
    df <- df[order(df$p.adjust), , drop = FALSE]
    df <- head(df, show)
    df$Label <- stringr::str_wrap(.pretty_set_label(df$Description), 34)
    p <- ggplot(df, aes(x = NES, y = reorder(Label, NES))) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = NES, yend = reorder(Label, NES)),
                   color = "grey80", linewidth = 0.5) +
      geom_point(aes(size = setSize, fill = p.adjust), shape = 21, color = "grey30") +
      scale_fill_gradient(low = "#b2182b", high = "#2166ac", trans = "log10",
                          name = "p.adjust",
                          labels = function(v) format(v, scientific = TRUE, digits = 1)) +
      scale_size_continuous(range = c(3, 10), name = "Set size") +
      labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(color = "black"),
            panel.grid.minor = element_blank())
  } else {
    p <- dotplot(res, showCategory = show) + ggtitle(title)
  }
  ggsave(path, plot = p, width = width, height = height, dpi = 600, bg = "white")
  invisible(p)
}

# GSEA enrichment plot (top n pathways)
save_gsea_plot <- function(res, title, path, n = 3) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  ids <- head(df$ID, n)
  p <- gseaplot2(res, geneSetID = ids, title = title)
  ggsave(path, plot = p, width = 12, height = 6 * ceiling(n / 2), dpi = 300)
  invisible(p)
}

# ------------------------------------------------------------------------------
# Driver (single cancer: kidney)
# ------------------------------------------------------------------------------
message("\n=== ", toupper(PREFIX), " ENRICHMENT ===")

meta <- import(META_FILE) |>
  mutate(Gene_ID = as.character(Gene_ID))
message("  Genes in meta: ", nrow(meta))

# REM output has no adjusted p-value; derive BH-adjusted padj from randomP.
meta <- meta |>
  mutate(rem_padj = p.adjust(randomP, method = "BH"))

# Universe: all tested genes with an Entrez ID
universe_entrez <- meta |> filter(!is.na(Gene_ID)) |> pull(Gene_ID) |> unique()

# Significant DEGs for ORA
sig <- meta |> filter(rem_padj < PADJ_CUTOFF, abs(randomSummary) >= LFC_CUTOFF)
message("  Significant DEGs (padj<", PADJ_CUTOFF, ", |LFC|>=", LFC_CUTOFF, "): ", nrow(sig))

entrez_up   <- sig |> filter(randomSummary > 0) |> pull(Gene_ID) |> unique()
entrez_down <- sig |> filter(randomSummary < 0) |> pull(Gene_ID) |> unique()
entrez_all  <- sig |> pull(Gene_ID) |> unique()
message("  Up: ", length(entrez_up), "  Down: ", length(entrez_down))

# ----- ORA -----
message("  Running ORA...")
ora_up   <- run_ora(entrez_up,   universe_entrez)
ora_down <- run_ora(entrez_down, universe_entrez)
ora_all  <- run_ora(entrez_all,  universe_entrez)

export_enrichment(ora_up,   OUT_CSV, paste0(PREFIX, "_ora_up"))
export_enrichment(ora_down, OUT_CSV, paste0(PREFIX, "_ora_down"))
export_enrichment(ora_all,  OUT_CSV, paste0(PREFIX, "_ora_all"))

for (name in names(ora_all)) {
  save_dotplot(ora_all[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_ora_all_", name, "_dotplot.png")))
}

# ----- GSEA -----
message("  Running GSEA...")
ranked <- build_ranked_list(meta)
message("  Ranked genes for GSEA: ", length(ranked))

gsea_res <- run_gsea(ranked)
export_enrichment(gsea_res, OUT_CSV, paste0(PREFIX, "_gsea"))

for (name in names(gsea_res)) {
  save_dotplot(gsea_res[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_", name, "_dotplot.png")))
  save_gsea_plot(gsea_res[[name]], paste(PREFIX, name, "- Top Pathways"),
                 file.path(OUT_FIG, paste0(PREFIX, "_", name, "_enrichplot.png")))
}
names(meta_input)
ls()
exists("gse24455")
exists("gse63183")
exists("combined_results")
exists("reml_results_annotated")





library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(rio)
project_dir <- "D:/kidney_cancer_meta_analysis"

meta_dir <- file.path(
  project_dir,
  "results/tables/meta-analysis"
)

figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment"
)

table_dir <- file.path(
  project_dir,
  "results/tables/enrichment"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

random_file <- file.path(
  meta_dir,
  "random_effect_model_CORRECTED.csv"
)

random_result <- rio::import(random_file)

names(random_result)

randomP
randomSummary
build_ranked_list <- function(meta_df) {
  
  df <- meta_df %>%
    filter(
      !is.na(meta_pvalue),
      !is.na(meta_log2FC),
      meta_pvalue > 0,
      !is.na(Gene_ID)
    ) %>%
    mutate(
      rank_score = sign(meta_log2FC) * -log10(meta_pvalue),
      Gene_ID = as.character(Gene_ID)
    ) %>%
    filter(!duplicated(Gene_ID)) %>%
    arrange(desc(rank_score))
  
  ranked_list <- df$rank_score
  names(ranked_list) <- df$Gene_ID
  
  ranked_list
}
ranked_genes <- build_ranked_list(random_result)

length(ranked_genes)
head(ranked_genes)

universe_entrez <- random_result %>%
  filter(!is.na(Gene_ID)) %>%
  pull(Gene_ID) %>%
  as.character() %>%
  unique()
length(universe_entrez)
head(universe_entrez)
run_ora <- function(entrez_ids, universe_entrez) {
  
  if (length(entrez_ids) < 5) {
    
    message(
      "Too few genes for ORA (n=",
      length(entrez_ids),
      "), skipping."
    )
    
    return(NULL)
  }
  
  results <- list()
  
  for (ont in c("BP", "MF", "CC")) {
    
    results[[paste0("GO_", ont)]] <- tryCatch(
      
      enrichGO(
        gene = entrez_ids,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
      ),
      
      error = function(e) {
        
        message(
          "GO-",
          ont,
          " error: ",
          e$message
        )
        
        NULL
      }
    )
  }
  
  return(results)
}
up_genes <- random_result %>%
  filter(
    !is.na(meta_padj),
    meta_padj < 0.05,
    meta_log2FC > 0
  ) %>%
  pull(Gene_ID) %>%
  as.character() %>%
  unique()
down_genes <- random_result %>%
  filter(
    !is.na(meta_padj),
    meta_padj < 0.05,
    meta_log2FC < 0
  ) %>%
  pull(Gene_ID) %>%
  as.character() %>%
  unique()
length(up_genes)
length(down_genes)
ora_up <- run_ora(
  entrez_ids = up_genes,
  universe_entrez = universe_entrez
)

ora_down <- run_ora(
  entrez_ids = down_genes,
  universe_entrez = universe_entrez
)
ora_up
ora_down
if (!is.null(ora_up$GO_BP)) {
  
  write.csv(
    as.data.frame(ora_up$GO_BP),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_UP_GO_BP.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_up$GO_MF)) {
  
  write.csv(
    as.data.frame(ora_up$GO_MF),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_UP_GO_MF.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_up$GO_CC)) {
  
  write.csv(
    as.data.frame(ora_up$GO_CC),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_UP_GO_CC.csv"
    ),
    row.names = FALSE
  )
}
if (!is.null(ora_down$GO_BP)) {
  
  write.csv(
    as.data.frame(ora_down$GO_BP),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_DOWN_GO_BP.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_down$GO_MF)) {
  
  write.csv(
    as.data.frame(ora_down$GO_MF),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_DOWN_GO_MF.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_down$GO_CC)) {
  
  write.csv(
    as.data.frame(ora_down$GO_CC),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_DOWN_GO_CC.csv"
    ),
    row.names = FALSE
  )
}
build_ranked_list()
ranked_genes <- build_ranked_list(random_result)
length(ranked_genes)
head(ranked_genes)
tail(ranked_genes)
up_genes
down_genes
summary(ranked_genes)

head(ranked_genes, 10)

tail(ranked_genes, 10)
rank_score = sign(meta_log2FC) * -log10(meta_pvalue)
build_ranked_list <- function(meta_df) {
  
  df <- meta_df %>%
    dplyr::filter(
      !is.na(meta_pvalue),
      !is.na(meta_log2FC),
      meta_pvalue > 0,
      !is.na(Gene_ID)
    ) %>%
    dplyr::mutate(
      rank_score = sign(meta_log2FC) * -log10(meta_pvalue),
      Gene_ID = as.character(Gene_ID)
    ) %>%
    dplyr::filter(!duplicated(Gene_ID)) %>%
    dplyr::arrange(desc(rank_score))
  
  setNames(df$rank_score, df$Gene_ID)
}
ranked_genes <- build_ranked_list(random_result)
length(ranked_genes)
head(ranked_genes)
tail(ranked_genes)

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rio)
library(ReactomePA)
library(msigdbr)
project_dir <- "D:/kidney_cancer_meta_analysis"

meta_dir <- file.path(
  project_dir,
  "results/tables/meta-analysis"
)

figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment"
)

table_dir <- file.path(
  project_dir,
  "results/tables/enrichment"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
gsea_GO_BP <- tryCatch({
  
  gseGO(
    geneList = ranked_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
}, error = function(e) {
  
  message("GO-BP GSEA error: ", e$message)
  NULL
  
})
if (!is.null(gsea_GO_BP)) {
  
  gsea_GO_BP_df <- as.data.frame(gsea_GO_BP)
  
  write.csv(
    gsea_GO_BP_df,
    file.path(
      table_dir,
      "kidney_random_REML_GSEA_GO_BP.csv"
    ),
    row.names = FALSE
  )
  
  message(
    "GO-BP GSEA pathways: ",
    nrow(gsea_GO_BP_df)
  )
}
if (!is.null(gsea_GO_BP) &&
    nrow(as.data.frame(gsea_GO_BP)) > 0) {
  
  p_gsea_GO_BP_dot <- dotplot(
    gsea_GO_BP,
    showCategory = 20,
    font.size = 10
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    labs(
      title = "Kidney Cancer — GSEA GO Biological Process",
      subtitle = "Random-effects model (REML)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_GO_BP_dotplot.png"
    ),
    p_gsea_GO_BP_dot,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_GO_BP_dotplot.pdf"
    ),
    p_gsea_GO_BP_dot,
    width = 10,
    height = 8,
    units = "in"
  )
}
if (!is.null(gsea_GO_BP) &&
    nrow(as.data.frame(gsea_GO_BP)) > 0) {
  
  top_GO_BP <- as.data.frame(gsea_GO_BP) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1)
  
  if (nrow(top_GO_BP) > 0) {
    
    pathway_GO_BP <- top_GO_BP$ID[1]
    
    p_gsea_GO_BP_enrich <- gseaplot2(
      gsea_GO_BP,
      geneSetID = pathway_GO_BP,
      title = top_GO_BP$Description[1],
      color = "#2166AC",
      base_size = 12
    )
    
    ggsave(
      file.path(
        figure_dir,
        "kidney_random_REML_GSEA_GO_BP_enrichplot.png"
      ),
      p_gsea_GO_BP_enrich,
      width = 10,
      height = 8,
      dpi = 400
    )
    
    ggsave(
      file.path(
        figure_dir,
        "kidney_random_REML_GSEA_GO_BP_enrichplot.pdf"
      ),
      p_gsea_GO_BP_enrich,
      width = 10,
      height = 8,
      units = "in"
    )
  }
}
hallmark_df <- tryCatch({
  
  msigdbr(
    species = "Homo sapiens",
    category = "H"
  )
  
}, error = function(e) {
  
  message("Old msigdbr syntax failed; trying collection='H'")
  
  tryCatch(
    
    msigdbr(
      species = "Homo sapiens",
      collection = "H"
    ),
    
    error = function(e2) {
      message(
        "Hallmark msigdbr error: ",
        e2$message
      )
      NULL
    }
  )
})

if (!is.null(hallmark_df)) {
  
  hallmark_term2gene <- hallmark_df %>%
    filter(
      !is.na(ncbi_gene),
      ncbi_gene != ""
    ) %>%
    select(
      gs_name,
      ncbi_gene
    ) %>%
    distinct()
  
  hallmark_term2name <- hallmark_df %>%
    select(
      gs_name,
      gs_description
    ) %>%
    distinct()
  
  hallmark_term2gene$ncbi_gene <-
    as.character(hallmark_term2gene$ncbi_gene)
}

if (!is.null(hallmark_df)) {
  
  gsea_Hallmark <- tryCatch({
    
    GSEA(
      geneList = ranked_genes,
      TERM2GENE = hallmark_term2gene,
      TERM2NAME = hallmark_term2name,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
  }, error = function(e) {
    
    message(
      "Hallmark GSEA error: ",
      e$message
    )
    
    NULL
  })
}
if (exists("gsea_Hallmark") &&
    !is.null(gsea_Hallmark)) {
  
  hallmark_gsea_df <-
    as.data.frame(gsea_Hallmark)
  
  write.csv(
    hallmark_gsea_df,
    file.path(
      table_dir,
      "kidney_random_REML_GSEA_Hallmark.csv"
    ),
    row.names = FALSE
  )
  
  message(
    "Hallmark GSEA pathways: ",
    nrow(hallmark_gsea_df)
  )
}

if (exists("gsea_Hallmark") &&
    !is.null(gsea_Hallmark) &&
    nrow(as.data.frame(gsea_Hallmark)) > 0) {
  
  p_Hallmark_dot <- dotplot(
    gsea_Hallmark,
    showCategory = 20,
    font.size = 10
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    labs(
      title = "Kidney Cancer — Hallmark GSEA",
      subtitle = "Random-effects model (REML)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_Hallmark_dotplot.png"
    ),
    p_Hallmark_dot,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_Hallmark_dotplot.pdf"
    ),
    p_Hallmark_dot,
    width = 10,
    height = 8,
    units = "in"
  )
}
if (exists("gsea_Hallmark") &&
    !is.null(gsea_Hallmark) &&
    nrow(as.data.frame(gsea_Hallmark)) > 0) {
  
  top_Hallmark <-
    as.data.frame(gsea_Hallmark) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1)
  
  pathway_Hallmark <-
    top_Hallmark$ID[1]
  
  p_Hallmark_enrich <- gseaplot2(
    gsea_Hallmark,
    geneSetID = pathway_Hallmark,
    title = top_Hallmark$Description[1],
    color = "#2166AC",
    base_size = 12
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_Hallmark_enrichplot.png"
    ),
    p_Hallmark_enrich,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_Hallmark_enrichplot.pdf"
    ),
    p_Hallmark_enrich,
    width = 10,
    height = 8,
    units = "in"
  )
}
gsea_KEGG <- tryCatch({
  
  gseKEGG(
    geneList = ranked_genes,
    organism = "hsa",
    keyType = "ncbi-geneid",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
}, error = function(e) {
  
  message(
    "KEGG GSEA error: ",
    e$message
  )
  
  NULL
})

if (!is.null(gsea_KEGG)) {
  
  kegg_gsea_df <-
    as.data.frame(gsea_KEGG)
  
  write.csv(
    kegg_gsea_df,
    file.path(
      table_dir,
      "kidney_random_REML_GSEA_KEGG.csv"
    ),
    row.names = FALSE
  )
  
  message(
    "KEGG GSEA pathways: ",
    nrow(kegg_gsea_df)
  )
}

if (!is.null(gsea_KEGG) &&
    nrow(as.data.frame(gsea_KEGG)) > 0) {
  
  p_gsea_KEGG_dot <- dotplot(
    gsea_KEGG,
    showCategory = 20,
    font.size = 10
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    labs(
      title = "Kidney Cancer — KEGG GSEA",
      subtitle = "Random-effects model (REML)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_KEGG_dotplot.png"
    ),
    p_gsea_KEGG_dot,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_KEGG_dotplot.pdf"
    ),
    p_gsea_KEGG_dot,
    width = 10,
    height = 8,
    units = "in"
  )
}

if (!is.null(gsea_KEGG) &&
    nrow(as.data.frame(gsea_KEGG)) > 0) {
  
  top_KEGG <-
    as.data.frame(gsea_KEGG) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1)
  
  pathway_KEGG <- top_KEGG$ID[1]
  
  p_gsea_KEGG_enrich <- gseaplot2(
    gsea_KEGG,
    geneSetID = pathway_KEGG,
    title = top_KEGG$Description[1],
    color = "#2166AC",
    base_size = 12
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_KEGG_enrichplot.png"
    ),
    p_gsea_KEGG_enrich,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      "kidney_random_REML_GSEA_KEGG_enrichplot.pdf"
    ),
    p_gsea_KEGG_enrich,
    width = 10,
    height = 8,
    units = "in"
  )
}
ora_up_KEGG <- tryCatch({
  
  enrichKEGG(
    gene = up_genes,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20
  )
  
}, error = function(e) {
  
  message(
    "UP KEGG ORA error: ",
    e$message
  )
  
  NULL
})

ora_down_KEGG <- tryCatch({
  
  enrichKEGG(
    gene = down_genes,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20
  )
  
}, error = function(e) {
  
  message(
    "DOWN KEGG ORA error: ",
    e$message
  )
  
  NULL
})
if (!is.null(ora_up_KEGG)) {
  
  write.csv(
    as.data.frame(ora_up_KEGG),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_UP_KEGG.csv"
    ),
    row.names = FALSE
  )
}

ora_up_Reactome <- tryCatch({
  
  enrichPathway(
    gene = up_genes,
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
}, error = function(e) {
  
  message(
    "UP Reactome ORA error: ",
    e$message
  )
  
  NULL
})



if (!is.null(ora_down_KEGG)) {
  
  write.csv(
    as.data.frame(ora_down_KEGG),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_DOWN_KEGG.csv"
    ),
    row.names = FALSE
  )
}
ora_up_Reactome <- tryCatch({
  
  enrichPathway(
    gene = up_genes,
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
}, error = function(e) {
  
  message(
    "UP Reactome ORA error: ",
    e$message
  )
  
  NULL
})
ora_down_Reactome <- tryCatch({
  
  enrichPathway(
    gene = down_genes,
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
}, error = function(e) {
  
  message(
    "DOWN Reactome ORA error: ",
    e$message
  )
  
  NULL
})
if (!is.null(ora_up_Reactome)) {
  
  write.csv(
    as.data.frame(ora_up_Reactome),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_UP_Reactome.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_down_Reactome)) {
  
  write.csv(
    as.data.frame(ora_down_Reactome),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_DOWN_Reactome.csv"
    ),
    row.names = FALSE
  )
}

save_ora_dotplot <- function(
    enrichment_object,
    filename,
    title_text) {
  
  if (is.null(enrichment_object)) {
    return(NULL)
  }
  
  if (nrow(as.data.frame(enrichment_object)) == 0) {
    return(NULL)
  }
  
  p <- dotplot(
    enrichment_object,
    showCategory = 20,
    font.size = 10
  ) +
    scale_color_gradient(
      low = "#2166AC",
      high = "#B2182B"
    ) +
    labs(
      title = title_text
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
  
  ggsave(
    file.path(
      figure_dir,
      paste0(filename, ".png")
    ),
    p,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  ggsave(
    file.path(
      figure_dir,
      paste0(filename, ".pdf")
    ),
    p,
    width = 10,
    height = 8,
    units = "in"
  )
  
  return(p)
}

all_sig_genes <- random_result %>%
  filter(
    !is.na(meta_padj),
    meta_padj < 0.05,
    !is.na(Gene_ID)
  ) %>%
  pull(Gene_ID) %>%
  as.character() %>%
ora_all_GO_CC <- enrichGO(
  gene = all_sig_genes,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

  unique()

length(all_sig_genes)
ora_all_GO_BP <- enrichGO(
  gene = all_sig_genes,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ora_all_GO_MF <- enrichGO(
  gene = all_sig_genes,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

write.csv(
  as.data.frame(ora_all_GO_BP),
  file.path(
    table_dir,
    "kidney_random_REML_ORA_ALL_GO_BP.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(ora_all_GO_CC),
  file.path(
    table_dir,
    "kidney_random_REML_ORA_ALL_GO_CC.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(ora_all_GO_MF),
  file.path(
    table_dir,
    "kidney_random_REML_ORA_ALL_GO_MF.csv"
  ),
  row.names = FALSE
)

ora_all_KEGG <- tryCatch({
  
  enrichKEGG(
    gene = all_sig_genes,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20
  )
  
}, error = function(e) {
  
  message(
    "ALL KEGG ORA error: ",
    e$message
  )
  
  NULL
})

ora_all_Reactome <- tryCatch({
  
  enrichPathway(
    gene = all_sig_genes,
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
}, error = function(e) {
  
  message(
    "ALL Reactome ORA error: ",
    e$message
  )
  
  NULL
})

if (!is.null(ora_all_KEGG)) {
  
  write.csv(
    as.data.frame(ora_all_KEGG),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_ALL_KEGG.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(ora_all_Reactome)) {
  
  write.csv(
    as.data.frame(ora_all_Reactome),
    file.path(
      table_dir,
      "kidney_random_REML_ORA_ALL_Reactome.csv"
    ),
    row.names = FALSE
  )
}

save_ora_dotplot(
  ora_all_GO_BP,
  "kidney_random_REML_ORA_ALL_GO_BP_dotplot",
  "Kidney Cancer — GO Biological Process ORA"
)
save_ora_dotplot(
  ora_all_GO_CC,
  "kidney_random_REML_ORA_ALL_GO_CC_dotplot",
  "Kidney Cancer — GO Cellular Component ORA"
)
save_ora_dotplot(
  ora_all_KEGG,
  "kidney_random_REML_ORA_ALL_KEGG_dotplot",
  "Kidney Cancer — KEGG ORA"
)
save_ora_dotplot(
  ora_all_Reactome,
  "kidney_random_REML_ORA_ALL_Reactome_dotplot",
  "Kidney Cancer — Reactome ORA"
)
# UP KEGG
save_ora_dotplot(
  ora_up_KEGG,
  "kidney_random_REML_ORA_UP_KEGG_dotplot",
  "Kidney Cancer — Upregulated KEGG ORA"
)
# DOWN KEGG
save_ora_dotplot(
  ora_down_KEGG,
  "kidney_random_REML_ORA_DOWN_KEGG_dotplot",
  "Kidney Cancer — Downregulated KEGG ORA"
)
save_ora_dotplot(
  ora_up_Reactome,
  "kidney_random_REML_ORA_UP_Reactome_dotplot",
  "Kidney Cancer — Upregulated Reactome ORA"
)
save_ora_dotplot(
  ora_down_Reactome,
  "kidney_random_REML_ORA_DOWN_Reactome_dotplot",
  "Kidney Cancer — Downregulated Reactome ORA"
)
cat("TABLE OUTPUTS:\n")

print(
  list.files(
    table_dir,
    full.names = FALSE
  )
)
cat("\nFIGURE OUTPUTS:\n")

print(
  list.files(
    figure_dir,
    full.names = FALSE
  )
)
cat("NUMBER OF TABLE FILES: ",
    length(list.files(table_dir)),
    "\n")
cat("NUMBER OF FIGURE FILES: ",
    length(list.files(figure_dir)),
    "\n")
