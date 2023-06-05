# Pan Cancer Analysis Reveals Modecular Patterns Associated With Age

Age is a known driver of various diseases, including cancer. While there is a higher incidence of cancer in the older population with median age at diagnosis being > 60 years, cancer can still occur at any age. Given the occurrence across varied age groups,in addition to various known genetic drivers of cancer there could be age-dependent molecular differences at baseline in patients within each cancer type. Based on this hypothesis it would be interesting to identify age-related biological changes, particularly in terms of immune functions/immunocompetence. Therefore, the idea would be to interrogate age-dependent transcriptional regulation patterns within each tumor type by integrating RNA-Seq and other genomic/clinical data available from the TCGA consortium. Broadly set of analysis would encompass:

* Differential expression using RNA-Seq dataset across TCGA cancer types with varied occurrence across a broad range of ages/pathological stage.
* Functional annotation of transcriptional changes using Gene Set Enrichment Analysis to understand global biological pathway changes
* Are age-specific changes confounded with molecular subtypes in any given cancer? -Compare aggregate age specific gene signature scores across known molecular subtypes within a given cancer.
* Age-dependent immune cell type profiles – using cell type deconvolution tools such as xCell and cibersort to analyze differential age-related immune cell type enrichment.
Integrating transcriptional patterns observed with epigenetic changes for matched TCGA patients – can we also predict age based on some of the epigenetic markers/ explain age-dependent transcriptional changes?
* Use independent datasets (e.g. Mouse aging organs data/GTEX) to validate age-dependent markers in an organ/cancer type-specific basis.
* Tumor mutation burden is known to be tissue and age specific already – to validate this trend across subsets
* Age and ethnicity distribution (contingent on available diversified cohort)

This work culminated into the following publication: https://pubmed.ncbi.nlm.nih.gov/34879281/

Data generated can be obtained from: https://zenodo.org/record/5586331#.YXCGTtlKg-W
