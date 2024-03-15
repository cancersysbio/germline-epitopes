## README
This directory includes source data supporting various controls. Each file is detailed below.

### TCGA Binding Thresholds
This is the source data for **Supplementary Figure 1i** demonstrating a negative association between GEB and HER2+ breast cancer considering various binding thresholds to determine GEB.\
**File:** tcga_her2_thresholds.txt

#### Columns
**sample:** sample name\
**wbs_0.5_2:** GEB considering "weak binders" as within 0.5-2% of naturally occuring random peptides\
**wbs_0.5_3:** GEB considering "weak binders" as within 0.5-3% of naturally occuring random peptides\
**wbs_0.5_4:** GEB considering "weak binders" as within 0.5-4% of naturally occuring random peptides\
**wbs_0.25_2:** GEB considering "weak binders" as within 0.25-2% of naturally occuring random peptides\
**wbs_0.25_3:** GEB considering "weak binders" as within 0.25-3% of naturally occuring random peptides\
**wbs_0.25_4:** GEB considering "weak binders" as within 0.25-4% of naturally occuring random peptides\
**pam50:** pam50 subtype\
**PC1-5:** genetic principal components 1-5 

### TCGA GEB using MHCflurry 
This is the source data for **Supplementary Figure 1j** demonstrating a negative association between GEB and HER2+ breast cancer calculating GEB using MHCflurry.\
**file:** tcga_her2_mhcflurry.txt

#### Columns
**sample:** sample name\
**wbs_0.5_2:** GEB considering "weak binders" as within 0.5-2% of naturally occuring random peptides by MHCflurry\
**wbs_0.5_3:** GEB considering "weak binders" as within 0.5-3% of naturally occuring random peptides by MHCflurry\
**wbs_0.5_5:** GEB considering "weak binders" as within 0.5-4% of naturally occuring random peptides by MHCflurry\
**wbs_0.25_2:** GEB considering "weak binders" as within 0.25-2% of naturally occuring random peptides by MHCflurry\
**wbs_0.25_3:** GEB considering "weak binders" as within 0.25-3% of naturally occuring random peptides by MHCflurry\
**wbs_0.25_5:** GEB considering "weak binders" as within 0.25-4% of naturally occuring random peptides by MHCflurry\

### Primary Null Associations
This is the source data for **Supplementary Figure 1k** demonstrating a negative association between GEB and somatic amplification goes away when calculating GEB with scrambled HLAs over 1,000 iterations.\
**file:** primary_null_associations.txt

#### Columns
**subtype:** sample name\
**gene:** gene(s)\
**coef:** coefficient\
**p:** p-value\
**se:** standard error\
**l95:** lower 95% confidence interval\
**u95:** upper 95% confidence interval\
**iteration:** iteration number

### Primary vs Metastasis Null Associations
This is the source data for **Supplementary Figure 3j** demonstrating a enrichment in of GEB metstatic tumors goes away when calculating GEB with scrambled HLAs over 1,000 iterations.\
**file:** primary_vs_metastatic_null_associations.txt

#### Columns
**subtype:** sample name\
**gene:** gene(s)\
**coef:** coefficient\
**p:** p-value\
**se:** standard error\
**l95:** lower 95% confidence interval\
**u95:** upper 95% confidence interval\
**iteration:** iteration number

### Variance explained by rare and common variants
The is the source data for **Supplementary Figure 1S** demonstrating majority of the variance explained is from common variants.\
**file:** rare_common_r2_estimates.txt

#### Columns
**subtype:** subtype\
**common:** variance explained (r2) by common variants
**rare:** variance explained (r2) by rare variants

