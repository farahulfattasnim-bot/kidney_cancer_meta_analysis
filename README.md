# Renal_Cell_Carcinoma_meta_analysis
Integrated Transcriptomic Meta-analysis Reveals 
Recurrent Gene Expression Alterations in Clear Cell 
Renal Cell Carcinoma
# Abstract
Background: Clear cell renal cell carcinoma (ccRCC) is the most common histological subtype 
of renal cell carcinoma and is characterized by substantial molecular and transcriptional 
alterations. Although individual transcriptomic studies have identified numerous genes and 
biological pathways associated with ccRCC, findings from individual cohorts may be influenced 
by study-specific sample characteristics, experimental platforms, and biological heterogeneity. 
Integrative analysis of independent transcriptomic datasets can therefore help identify expression 
changes that are consistently observed across studies.

 # Introduction
Renal cell carcinoma (RCC) is a major malignancy of the urinary system, with clear cell renal 
cell carcinoma (ccRCC) being its most common histological subtype. ccRCC exhibits extensive 
molecular alterations involving cellular metabolism, hypoxia signaling, angiogenesis, epigenetic 
regulation, and tumor-cell proliferation [1,2]. In particular, dysregulation of the VHL–hypoxia￾inducible factor (HIF) axis represents a central molecular feature of ccRCC development and 
progression [1–3]. Characterizing the resulting transcriptional alterations is therefore important 
for understanding the molecular biology of ccRCC.
High-throughput RNA sequencing (RNA-seq) enables genome-wide characterization of gene￾expression changes, while public repositories such as the NCBI Gene Expression Omnibus 
(GEO) provide opportunities to reanalyze and integrate independent transcriptomic datasets [4]. 
However, differences in sample characteristics, sequencing platforms, and experimental 
conditions can produce study-specific expression signatures, limiting the reproducibility of 
individual findings.
Transcriptomic meta-analysis provides a framework for integrating evidence across independent 
studies while accounting for variation in study-specific effects. Rather than directly pooling raw 
expression data generated under different experimental conditions, study-wise differential 
expression followed by random-effects meta-analysis can provide a more robust estimate of 
shared gene-expression effects [5]. Appropriate multiple-testing correction is essential for 
reliable genome-wide inference [6].
In this study, two human ccRCC transcriptomic datasets, GSE24455 and GSE63183, were 
integrated. GSE24455 comprises 10 paired ccRCC tumor and adjacent normal tissues, whereas 
GSE63183 contains two ccRCC patients with matched tumor and normal tissues [7,8]. The study 
aimed to identify recurrent gene-expression alterations in ccRCC by performing study-wise 
differential expression, gene-identifier harmonization, and random-effects meta-analysis, 
followed by functional enrichment and protein-interaction network analysis. Given the limited 
sample size of GSE63183, the findings were considered exploratory and interpreted with 
appropriate caution regarding study-specific effects and heterogeneity

# Discussion
This integrated meta-analysis demonstrates that kidney cancer is characterized by coordinated 
disruption of metabolic, immune-inflammatory, renal epithelial, and oncogenic signaling 
programs. The random-effects analysis identified 3,386 genes at FDR <0.05, including 2,239 
with |log2FC| ≥1, indicating a substantial and reproducible transcriptional shift across datasets.
A major finding was metabolic reprogramming, with enrichment of amino-acid, organic-acid, 
glycolytic, TCA-cycle, lipid, and mitochondrial pathways. Network modules containing GPT, 
EPHX2, AGXT2, PKM, TPI1, ALDOA, PGK1, ACLY, and PCK1 further supported 
disruption of central metabolism. The simultaneous representation of reduced 
oxidative/metabolic functions and increased glycolytic/hypoxia-associated programs suggests 
adaptation of tumor cells toward altered energy utilization.
The analysis also revealed strong immune and inflammatory remodeling, including interferon, 
TNF, phagocytosis, cytokine, and antigen-processing signatures. The large immune/signaling 
module containing CCL5, CD8A, FCGR3A, B2M, CXCR4, EGFR, and ERBB2 supports 
substantial interaction between immune responses and tumor-associated signaling. Infection- and 
transplant-rejection-related gene sets should be interpreted as shared immune/host-response 
signatures rather than direct evidence of infection or rejection.
Kidney-specific functional alterations were additionally evident through enrichment of ion 
transport, epithelial membrane, and acid–base regulatory processes, supported by Module 1 
genes such as KCNJ10, SLC26A7, SLC4A4, SCNN1A, and CA2. Network analysis further 
identified ALB, EGFR, TPI1, ERBB2, PKM, ACLY, PGK1, CD8A, B2M, CXCR4, and 
CAV1 as prominent hubs, linking metabolic, immune, signaling, and epithelial processes.
Overall, the findings support a multidimensional molecular phenotype of kidney cancer 
involving metabolic reprogramming, immune activation, altered renal epithelial function, 
hypoxia/glycolytic adaptation, and oncogenic signaling. The identified hub genes represent 
promising candidates for further investigation, although their network centrality does not 
establish causal function and requires experimental validation.

# References
1. Cancer Genome Atlas Research Network. Comprehensive molecular characterization of 
clear cell renal cell carcinoma. Nature. 2013;499:43–49. doi:10.1038/nature12222. 
2. Linehan WM, Ricketts CJ. The metabolic basis of kidney cancer. Seminars in Cancer 
Biology. 2013;23(1):46–55. doi:10.1016/j.semcancer.2012.06.002. 
3. Kaelin WG Jr. The von Hippel-Lindau tumour suppressor gene: insights into oxygen 
sensing and cancer. Nature Reviews Cancer. 2008;8:865–873. doi:10.1038/nrc2502. 
4. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression 
and hybridization array data repository. Nucleic Acids Research. 2002;30(1):207–210. 
doi:10.1093/nar/30.1.207. 
5. National Center for Biotechnology Information. Gene Expression Omnibus: 
GSE24455, Digital gene expression sequencing of 10 pairs samples between kidney 
normal tissue and cancer tissue. NCBI GEO. 
6. National Center for Biotechnology Information. Gene Expression Omnibus: 
GSE63183, Loss of 5-hydroxymethylcytosine is linked to gene body hypermethylation in 
kidney cancer. NCBI GEO. 
7. Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of 
Statistical Software. 2010;36(3):1–48. doi:10.18637/jss.v036.i03. 
8. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful 
approach to multiple testing. Journal of the Royal Statistical Society: Series B. 
1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
