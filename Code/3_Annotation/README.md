# Annotation
*updated 5/8/2021*

| **Filename**		| **Description** |
|:----------------------|:----------------|
| annotate_combp.R	| Annotates the comb-p results with overlapping and flanking genes from Gencode for both analyses |
| annotate_function.R	| Funtion to annoate the chromosome results files with overlapping and flanking genes from Gencode (can be called for either analysis |
| Download_code.sh	| Downloads and unzips the Gencode v31 hg38 annotation information |
| Extract_genes.R	| Extracts the relevant gene information from the Gencode file |

**Workflow:**
1) Download_code.sh
2) Extract_genes.R
3) annotate_function.R
4) annotate_combp.R