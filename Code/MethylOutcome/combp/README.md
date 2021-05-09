# Methylation Outcome Regional Analysis
*updated 5/8/2021*

| **Filename**			| **Description** |
|:------------------------------|:----------------|
| Comb_p_script.sh		| Sorts the input files using Bedtools and performs the comb-p regional analysis |
| combp_sites.R			| Merges the comb-p FDR single site results with the comb-p region results |
| Prep_data.R			| Re-formats the single site results for comb-p and then edits the comb-p results |

**Workflow:**
1) Prep_data.R until prompted to use Bedtools
2) Comb_p_script.sh until prompted to return to R and finish editing
3) Prep_data.R
4) combp_sites.R