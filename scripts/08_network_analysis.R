
# Output file
##########################################################

output_file <- file.path(
  
  network_figure_dir,
  
  paste0(
    prefix,
    "_module_",
    plot_number,
    "_network.png"
  )
)


png(
  
  filename =
    output_file,
  
  width = 2400,
  
  height = 2000,
  
  res = 300
)


set.seed(
  5000 + plot_number
)


module_layout <-
  
  if (vcount(module_graph) > 1) {
    
    layout_with_fr(
      module_graph,
      weights =
        E(module_graph)$score
    )
    
  } else {
    
    matrix(
      c(0, 0),
      ncol = 2
    )
  }


plot(
  
  module_graph,
  
  layout =
    module_layout,
  
  vertex.color =
    V(module_graph)$color,
  
  vertex.size =
    V(module_graph)$size,
  
  vertex.frame.color =
    V(module_graph)$frame.color,
  
  vertex.frame.width =
    V(module_graph)$frame.width,
  
  vertex.label =
    
    ifelse(
      
      V(module_graph)$ModuleHub,
      
      V(module_graph)$Gene_Symbol,
      
      NA
    ),
  
  vertex.label.cex =
    0.85,
  
  vertex.label.color =
    "#111111",
  
  edge.color =
    adjustcolor(
      "#777777",
      alpha.f = 0.35
    ),
  
  edge.width =
    0.9,
  
  main =
    paste(
      "STRING PPI Module",
      plot_number
    )
)


legend(
  
  "topleft",
  
  legend =
    c(
      "Upregulated",
      "Downregulated",
      "Module hub"
    ),
  
  pch =
    21,
  
  pt.bg =
    c(
      "#E41A1C",
      "#377EB8",
      "#F7F7F7"
    ),
  
  col =
    "#111111",
  
  pt.cex =
    1.5,
  
  bty = "n"
)


dev.off()


message(
  "Created: ",
  basename(output_file)
)
}
# STRING-PPI figure directory
network_figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment/STRING_PPI"
)

dir.create(
  network_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
project_dir <- "D:/kidney_cancer_meta_analysis"
network_figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment/STRING_PPI"
)

project_dir
network_figure_dir
dir.exists(network_figure_dir)
figure_dir

file.path(
  network_figure_dir,
  paste0(
    "STRING_PPI_Module_",
    plot_number,
    "_Module_",
    module_number,
    ".png"
  )
)
project_dir <- "D:/kidney_cancer_meta_analysis"

network_figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment/STRING_PPI"
)

dir.create(
  network_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
create_module_plot <- function(
    module_number,
    plot_number,
    graph
) {
  
  # Check graph
  if (!igraph::is_igraph(graph)) {
    stop("ERROR: graph must be an igraph object.")
  }
  
  # Check Module attribute
  if (!"Module" %in% igraph::vertex_attr_names(graph)) {
    stop("ERROR: Module attribute is missing from the graph.")
  }
  
  # Select module nodes
  module_values <- igraph::V(graph)$Module
  
  module_nodes <- igraph::V(graph)[
    !is.na(module_values) &
      module_values == module_number
  ]
  
  # Check module
  if (length(module_nodes) == 0) {
    message(
      "Module ",
      module_number,
      " contains no nodes. Skipping."
    )
    return(NULL)
  }
  
  # Create subgraph
  subgraph <- igraph::induced_subgraph(
    graph,
    vids = module_nodes
  )
  
  # Layout
  set.seed(123)
  
  module_layout <- igraph::layout_with_fr(
    subgraph
  )
  
  # --------------------------------
  # Node colors
  # --------------------------------
  
  node_colors <- ifelse(
    is.na(V(subgraph)$Regulation),
    "#BDBDBD",
    ifelse(
      V(subgraph)$Regulation == "Upregulated",
      "#E41A1C",
      ifelse(
        V(subgraph)$Regulation == "Downregulated",
        "#377EB8",
        "#BDBDBD"
      )
    )
  )
  
  # --------------------------------
  # Node size
  # --------------------------------
  
  node_sizes <- ifelse(
    !is.na(V(subgraph)$Hub) &
      V(subgraph)$Hub == TRUE,
    10,
    6
  )
  
  # --------------------------------
  # Labels
  # --------------------------------
  
  node_labels <- ifelse(
    !is.na(V(subgraph)$Hub) &
      V(subgraph)$Hub == TRUE &
      !is.na(V(subgraph)$Gene_Symbol),
    V(subgraph)$Gene_Symbol,
    NA
  )
  
  # --------------------------------
  # PNG output
  # --------------------------------
  
  png_file <- file.path(
    network_figure_dir,
    paste0(
      "STRING_PPI_Module_",
      plot_number,
      "_Module_",
      module_number,
      ".png"
    )
  )
  
  png(
    filename = png_file,
    width = 3000,
    height = 3000,
    res = 400
  )
  
  plot(
    subgraph,
    layout = module_layout,
    vertex.label = node_labels,
    vertex.label.cex = 0.8,
    vertex.label.color = "black",
    vertex.color = node_colors,
    vertex.size = node_sizes,
    vertex.frame.color = "black",
    vertex.frame.width = 1,
    edge.color = "#BDBDBD",
    edge.width = 0.8,
    main = paste0(
      "STRING PPI Module ",
      module_number
    )
  )
  
  dev.off()
  
  # --------------------------------
  # PDF output
  # --------------------------------
  
  pdf_file <- file.path(
    network_figure_dir,
    paste0(
      "STRING_PPI_Module_",
      plot_number,
      "_Module_",
      module_number,
      ".pdf"
    )
  )
  
  pdf(
    file = pdf_file,
    width = 8,
    height = 8
  )
  
  plot(
    subgraph,
    layout = module_layout,
    vertex.label = node_labels,
    vertex.label.cex = 0.8,
    vertex.label.color = "black",
    vertex.color = node_colors,
    vertex.size = node_sizes,
    vertex.frame.color = "black",
    vertex.frame.width = 1,
    edge.color = "#BDBDBD",
    edge.width = 0.8,
    main = paste0(
      "STRING PPI Module ",
      module_number
    )
  )
  
  dev.off()
  
  message(
    "Module ",
    module_number,
    " completed successfully."
  )
  
  return(subgraph)
}
test_module <- create_module_plot(
  module_number = selected_modules$Module[1],
  plot_number = selected_modules$Plot_Number[1],
  graph = ppi_graph
)
############################################################
# RECOVER STRING-PPI NETWORK OBJECTS
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

############################################################
# 1. PROJECT DIRECTORIES
############################################################

project_dir <- "D:/kidney_cancer_meta_analysis"

figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment"
)

table_dir <- file.path(
  project_dir,
  "results/tables/enrichment"
)

network_figure_dir <- file.path(
  figure_dir,
  "STRING_PPI"
)

dir.create(
  network_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. REBUILD PATHWAY-GENE TABLE
############################################################

# This requires gsea_bp_df to exist.
# If it exists, run:

pathway_gene_table <- gsea_bp_df %>%
  dplyr::filter(
    !is.na(ID),
    !is.na(Description),
    !is.na(core_enrichment),
    core_enrichment != ""
  ) %>%
  dplyr::select(
    Pathway_ID = ID,
    Pathway = Description,
    core_enrichment
  ) %>%
  tidyr::separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  dplyr::rename(
    Gene_ID = core_enrichment
  ) %>%
  dplyr::distinct()

pathway_gene_table <- pathway_gene_table %>%
  mutate(
    Gene_ID = as.character(Gene_ID),
    Pathway = as.character(Pathway)
  )

############################################################
# 3. CREATE NETWORK EDGES
############################################################

network_edges <- pathway_gene_table %>%
  dplyr::select(
    Pathway,
    Gene_ID
  ) %>%
  dplyr::distinct()

############################################################
# 4. CREATE PATHWAY NODES
############################################################

pathway_nodes <- data.frame(
  name = unique(network_edges$Pathway),
  type = "Pathway",
  stringsAsFactors = FALSE
)

############################################################
# 5. CREATE GENE NODES
############################################################

gene_nodes <- data.frame(
  name = unique(network_edges$Gene_ID),
  type = "Gene",
  stringsAsFactors = FALSE
)

############################################################
# 6. COMBINE NODES
############################################################

network_nodes <- dplyr::bind_rows(
  pathway_nodes,
  gene_nodes
)

############################################################
# 7. CREATE IGRAPH
############################################################

ppi_graph <- igraph::graph_from_data_frame(
  d = network_edges,
  vertices = network_nodes,
  directed = FALSE
)

############################################################
# CHECK
############################################################

class(ppi_graph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rio)

project_dir <- "D:/kidney_cancer_meta_analysis"

meta_dir <- file.path(
  project_dir,
  "results/tables/meta-analysis"
)

table_dir <- file.path(
  project_dir,
  "results/tables/enrichment"
)

network_figure_dir <- file.path(
  project_dir,
  "results/figures/enrichment/STRING_PPI"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_figure_dir, recursive = TRUE, showWarnings = FALSE)

igraph::is_igraph(ppi_graph)

igraph::vcount(ppi_graph)

igraph::ecount(ppi_graph)
objects_to_check <- c(
  "random_result",
  "ranked_genes",
  "gsea_GO_BP",
  "gsea_bp_df",
  "pathway_gene_table",
  "ppi_graph",
  "metric_df",
  "module_size_df",
  "selected_modules"
)

data.frame(
  Object = objects_to_check,
  Exists = sapply(objects_to_check, exists)
)
pathway_gene_table <- gsea_bp_df %>%
  dplyr::filter(
    !is.na(ID),
    !is.na(Description),
    !is.na(core_enrichment),
    core_enrichment != ""
  ) %>%
  dplyr::select(
    Pathway_ID = ID,
    Pathway = Description,
    core_enrichment
  ) %>%
  tidyr::separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  dplyr::rename(
    Gene_ID = core_enrichment
  ) %>%
  dplyr::distinct()
gsea_GO_BP <- ...
gsea_bp_df <- as.data.frame(gsea_GO_BP)
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_GO_BP <- tryCatch({
  
  clusterProfiler::gseGO(
    geneList = ranked_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0
  )
  
}, error = function(e) {
  
  message("GO-BP GSEA error: ", e$message)
  
  NULL
})
library(dplyr)
library(rio)

project_dir <- "D:/kidney_cancer_meta_analysis"

meta_dir <- file.path(
  project_dir,
  "results/tables/meta-analysis"
)

random_file <- file.path(
  meta_dir,
  "random_effect_model_CORRECTED.csv"
)

random_result <- rio::import(random_file)

dim(random_result)
names(random_result)
build_ranked_list <- function(meta_df) {
  
  df <- meta_df %>%
    dplyr::filter(
      !is.na(meta_pvalue),
      !is.na(meta_log2FC),
      meta_pvalue > 0,
      !is.na(Gene_ID)
    ) %>%
    dplyr::mutate(
      Gene_ID = as.character(Gene_ID),
      rank_score = sign(meta_log2FC) *
        -log10(meta_pvalue)
    ) %>%
    dplyr::filter(!duplicated(Gene_ID)) %>%
    dplyr::arrange(desc(rank_score))
  
  ranked_list <- df$rank_score
  names(ranked_list) <- df$Gene_ID
  
  return(ranked_list)
}

ranked_genes <- build_ranked_list(random_result)
length(ranked_genes)

head(ranked_genes)

tail(ranked_genes)
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_GO_BP <- tryCatch({
  
  clusterProfiler::gseGO(
    geneList = ranked_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0
  )
  
}, error = function(e) {
  
  message("GO-BP GSEA error: ", e$message)
  
  NULL
})