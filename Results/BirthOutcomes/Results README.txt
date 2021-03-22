#  NAME			DESCRIPTION

1. methyl		chr.pos of the CpG
2. methyl.p		methylation p-value 
3. arm2.p		arm2 p-value 
4. arm3.p		arm3 p-value 
5. BMI.p		BMI p-value 
6. age.p		arm p-value 
7. expose.p		smoke exposure p-value
8. PC1.p		genetic PC1 p-value
9. methyl.beta		methylation beta coefficient 
10. arm2.beta		arm2 beta coefficient 
11. arm3.beta		arm3 beta coefficient 
12. BMI.beta		BMI beta coefficient 
13. age.beta		age beta coefficient 
14. expose.beta		smoke exposure beta coefficient 
15. PC1.beta		genetic PC1 beta coefficient
16. NAs			number of samples with NA methylation values at that CpG (excluding outliers)
17. outliers		number of samples considered an outlier at that CpG (> Q3 + 3*IQR or < Q1 - 3*IQR)
18. total.NAs		sum of NAs and outliers
19. zeros		number of samples with 0% methylation at that CpG
20. non.zeros		number of samples with nonzero methylation values at that CpG
21. ones		number of samples with 100% methylation at that CpG
22. non.ones		number of samples with non-one methylation values at that CpG
23. diff.miss		pvalue for the calculation of differential missingness
24. ME			metastable epiallele identifier (. = no, ME = yes)
25. genes		genes that directly overlap the CpG
26. genes_flanking	genes within 5kb of the CpG
27. gene_type		gene type (Ensembl or HAVANA) of the genes listed in the genes_flanking column
28. methyl.fdr		FDR adjusted p-value for methylation

*Note:	total.NAs + zeros + non.zeros = 83 (total samples)
	total.NAs + ones + non.ones = 83
