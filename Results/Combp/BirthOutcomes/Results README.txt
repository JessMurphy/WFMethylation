#  NAME			DESCRIPTION

1. chrom		chromosome number
2. start		chromosome position for the beginning of the region
3. end  		chromosome position for the end of the region
4. min_p		minimum p-value of the original CpG's in the region?
5. n_probes 		number of CpG sites in the region
6. z_p			p-values for the region before the Sidak correction?
7. z_sidak_p 		reported p-values after the Sidak correction
8. genes		genes that directly overlap the region
9. genes_flanking	genes within 5kb of the region
10. gene_type		gene type (Ensembl or HAVANA) of the genes listed in the genes_flanking column

Note: A Sidak correction is the total bases covered by all input probes divided by the size of the given region
 
