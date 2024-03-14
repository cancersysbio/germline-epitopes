## README
This directory includes source data for Supplementary Figures 1c&d (TCGA), Supplementary Figure 2c&e (ICGC), Supplementary Figure 2d&f (METABRIC), Supplementary Figure 2j&k (DCIS) and Supplementary Figure 3a&b (Hartwig). It also includes source data for Supplementary Figures 3e-h. The source data includes variant allele frequencies confirmed against gnomAD and HLA allele frequencies confirmed against the Allele Frequency Net Database.

### Variant allele frequencies 
Minor allele frequencies calculated in the cohort of interest and found in gnomAD\
**file:** *cohort*_maf_gnomad.txt

#### Columns
chr: chromosome\
start: position start\
gene: gene\
snp: snp id\
maf: minor allele frequency in cohort\
end: position end (same as start)\
gnomad: maf in gnomAD

### HLA allele frequencies in cohort of interest
HLA allele frequencies calculated in the cohort of interest\
**file:** *cohort*_hla_frequencies.txt

#### Columns
allele: HLA allele
freq: frequency in cohort

### HLA allelel frequencies in Allele Frequency Net Database
HLA allele frequencies downloaded from the Allele Frequency Net Database\
**file:** population_hla_a_frequencies.csv\
**file:** population_hla_b_frequencies.csv\
**file:** population_hla_c_frequencies.csv

#### Columns
allele: HLA allele\
frequency: frequency in cohort
