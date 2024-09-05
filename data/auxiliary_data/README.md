## README
This directory includes source data for auxillary analyses as detailed below. 

### HLA alleles
The following file provides HLA alleles for TCGA (discovery) cohort as provided by Shulka et al.\
Please download from Shulka et al. (Supplementary Table 11) \
**file:** Shukla_Wu_Getz_Polysolver_HLA_Types_2015.tsv

The following file provides HLA alleles for ICGC (replication) cohort.\
**file:** icgc_hlas.txt

The following file provides HLA alleles for METABRIC (replication) cohort:\
**file:** metabric_hlas.txt

The following file provides HLA alleles that binding GP2 and E75 used as source data in Figure 1c and Supplementary Figure 1b.\
**file:** gp2_e75_hlas.txt

The following file provides imputation accuracy estimates for HLA alleles imputed from SNP6 arrays compared to HLA inferred from WES by Shulka et al.\
**file:** tcga_cookhla_hla_alleles_imputation_accuracy.txt
#### Columns
hla: hla\
correct: number of samples with hla where SNP6 and WES agree\
total: number of samples with HLA called by WES\
accuracy: accuracy defined as correct/total

The following file provides proportion of HLA alleles estimated to have high "bindng promiscuity" according to Manczinger et al.\
**file:** hla_promiscuity_proportions.txt
#### Columns
sample: sample id (METABRIC)\
high: number of HLA alleles considered "highly promiscuous", i.e. a wide peptide repertoire breadth\
notna: number of HLA alleles with binding repertoire measurements provided by Manczinger et al.\
unique: number of unique HLA alleles\
prompprop: proportion of HLA alleles considered "highly promiscuous", calculated as high/notna

### Allele specific amplification 
The following files provide source data Supplementary Figure 2m&n.\
**file:** erbb2_snp_read_depths.txt\
**file:** tubd1_snp_read_depths.txt

#### Columns
sample: sample\
nref: number of reads supporting reference allele in normal\
nalt: number of reads supporting alternate allele in normal\
tref: number of reads supporting reference allele in tumor\
talt: number of reads supporting alternate allele in tumor\
normal_ratio: ratio of reads supporting alternate vs reference allele in normal\
tumor_ratio: ratio of reads supporting alternate vs reference allele in tumor\
index: index for plotting purposes

The following files provide source data Figures 2b&c and Supplementary Figure 2o&p.\
**file:** erbb2_snp_peptides.txt\
**file:** tubd1_snp_peptides.txt

#### Columns
var_nmer: variant peptide\
hla: hla\
chr: chromosome of snp\
pos: position of snp\
var_rank: rank of binding affinity of variant peptide\
wt_nmer: wildtype peptide\
wt_rank: rank of binding affinity of wildtype peptide\
rank_diff: difference in rank between variant and wildtype peptide\
sample: sample

### GWAS summary stats
The following provides breast cancer GWAS summary stats from Zhang et al. used as source data for Supplementary Figure 1r.\
**file:** gwas_snps_summary_stats.txt

#### Columns
snp: snp\
gene:gene\
beta.GWAS: GWAS coefficient\
P.value.Gwas: GWAS p-value

### MIBI myoepithelial integrity
The following file provides source data for Figure 4f.\
**file:** MIBI_ECAD_data.txt

#### Columns
sample: sample\
Myoep_cluster_fracECAD: fraction of ECAD in myoepithelial cluster as measurement of myoepithelial integrity

### METABRIC IMC 
The following file provides source data for Supplementary Figure 4l. Cell types frequencies from Ali et al.\
**file:** metabric_imc_proportions.txt

#### Columns
sample: sample\
celltype: cell type\
frequency: number of cells corresponding to cell type\
propotion: proportion of total cells corresponding to cell type

### References
S. A. Shukla et al., Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol 33, 1152–1158 (2015).\
H. Zhang et al., Genome-wide association study identifies 32 novel breast cancer susceptibility loci from overall and subtype-specific analyses. Nat Genet 52, 572–581 (2020).\
M. Manczinger et al., Negative trade-off between neoantigen repertoire breadth and the specificity of HLA-I molecules shapes antitumor immunity. Nat Cancer 2, 950–961 (2021).\
H. R. Ali et al., Imaging mass cytometry and multiplatform genomics define the phenogenomic landscape of breast cancer. Nat Cancer 1, 163–175 (2020).



