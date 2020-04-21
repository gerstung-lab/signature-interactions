# Codes accompanying *Mutational signatures are jointly shaped by DNA damage and repair*

Codes and data for the manuscript "Mutational signatures are jointly shaped by DNA damage and repair" (accepted in Nature Communications 2020). `R-markdown` report is available at https://nvolkova.github.io/signature-interactions/ .

## Contents:

### worms
Codes and pre-calculated data to filter and classify mutations in *C. elegans* data and calculate the contributions of DNA repair and genotoxins to mutational spectra.

### cancer
Codes for analysing DNA damage-repair interactions in human cancers. Catalogues mutations in DNA repair genes, quantifies selection with dNdS analysis. This part requires a matrix with TCGA variants, and a metadata matrix with sample names, cancer types, and FPKM expression values for DNA repair genes.

### Supplementary_tables

Tables with signatures and interaction coefficients for *C. elegans* analysis, as well as with DNA repair genes and their mutations and selective pressure across TCGA samples.

### comparison

Compares genotoxin signatures in human cancers (as per COSMIC catalogue) and *C. elegans*.

### docs

Materials for `rmarkdown` report.

## R session information

Session information for the analyses:

```
R version 3.5.1 (2018-07-02)

Matrix products: default

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base   

other attached packages:
 [1] deconstructSigs_1.8.0                   org.Hs.eg.db_3.7.0                      xlsx_0.6.1                             
 [4] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GenomicFeatures_1.34.8                  AnnotationDbi_1.44.0                   
 [7] BSgenome.Hsapiens.UCSC.hg19_1.4.0       BSgenome_1.50.0                         rtracklayer_1.42.1                     
[10] dndscv_0.0.1.0                          bayesplot_1.6.0                         RColorBrewer_1.1-2                     
[13] greta_0.2.3                             reshape2_1.4.3                          ggplot2_3.1.0                          
[16] VariantAnnotation_1.28.10               Rsamtools_1.34.1                        Biostrings_2.50.2                      
[19] XVector_0.22.0                          SummarizedExperiment_1.12.0             DelayedArray_0.8.0                    
[22] BiocParallel_1.16.5                     matrixStats_0.54.0                      Biobase_2.42.0                         
[25] GenomicRanges_1.34.0                    GenomeInfoDb_1.18.1                     IRanges_2.16.0                         
[28] S4Vectors_0.20.1                        BiocGenerics_0.28.0     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0               lattice_0.20-35          prettyunits_1.0.2
 [4] assertthat_0.2.0         digest_0.6.18            R6_2.4.0
 [7] plyr_1.8.4               RSQLite_2.1.1            httr_1.4.0
[10] pillar_1.3.1             zlibbioc_1.28.0          rlang_0.3.1
[13] GenomicFeatures_1.34.3   progress_1.2.0           lazyeval_0.2.1
[16] blob_1.1.1               Matrix_1.2-14            stringr_1.3.1
[19] RCurl_1.95-4.11          bit_1.1-14               biomaRt_2.34.2
[22] munsell_0.5.0            compiler_3.5.1           pkgconfig_2.0.2
[25] tidyselect_0.2.5         tibble_2.0.1             GenomeInfoDbData_1.2.0
[28] codetools_0.2-15         XML_3.98-1.16            crayon_1.3.4
[31] dplyr_0.8.0.1            withr_2.1.2              GenomicAlignments_1.18.1
[34] MASS_7.3-50              bitops_1.0-6             grid_3.5.1
[37] gtable_0.2.0             DBI_1.0.0                magrittr_1.5
[40] scales_1.0.0             stringi_1.2.4            seqinr_3.4-5
[43] tools_3.5.1              ade4_1.7-13              bit64_0.9-7
[46] glue_1.3.0               purrr_0.3.0              hms_0.4.2
[49] AnnotationDbi_1.44.0     colorspace_1.4-0         memoise_1.1.0
[52] reticulate_1.11.1        listenv_0.7.0            base64enc_0.1-3
[55] tensorflow_1.10          future_1.12.0            coda_0.19-2
[58] tfruns_1.4               jsonlite_1.6             whisker_0.3-2
[61] globals_0.12.4
```
