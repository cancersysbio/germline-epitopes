## README
This directory includes cohort summary tables to generate the majority of figures in the manuscript. The cohort summary tables are detailed below.

### TCGA
Summary table for TCGA (discovery) cohort.\
**file:** tcga_megatable.txt

#### Columns
sample: sample name\
GEB in genes of interest (ERBB2 through MIA): GEB in each gene of interest\
ER: ER status by mRNA\ 
HER2: HER2 status by mRNA\
pam50: pam50 subtype\
ic10: ic10 subtype by iC10\
PC1-6: principal components 1 through 6\
CNA of genes of interest (ERBB2_CNA through RPS6KB1_CNA): binary indication if sample has amplification in gene of interest defined as CN > 4\
somatic: number of somatic SNVs\
mRNA of immune genes (B2M through TGFB1): mRNA abundance of gene of interest\
cytotoxic_score: geometric mean of mRNA of GZMA and PRF1\
Macrophages.M2 through T.cells.CD8: immune transcriptional scores from Thorsson et al.\
HER2.newly.derived: HER2 status by IHC as defined by Thennavan et al.\
ERBB2_CNA_5copies: binary indicating amplification of ERBB2 considering gene amplified if CN > 5\
TME_subtype: TME subtype as defined by Bagaev et al.\
stage: tumor stage\
age: age at diagnosis\
Triple.Negative.Status: TNBC status as defined by Thennavan et al.\

#### References
A. Thennavan et al., Molecular analysis of TCGA breast cancer histologic types. Cell Genom 1, 100067 (2021).\
A. Bagaev et al., Conserved pan-cancer microenvironment subtypes predict response to immunotherapy. Cancer Cell 39, 845-865.e7 (2021).\
V. Thorsson et al., The Immune Landscape of Cancer. Immunity 48, 812-830.e14 (2018).\

### ICGC
Summary table for ICGC (replication) cohort.\
**file:** icgc_megatable.txt

#### Columns
sample_name: sample name\
Individual.ID: individual ID\
GEB in genes of interest (ERBB2 through MELK): GEB in each gene of interest\
PC1-6: principal components 1 through 6\
final.ER: ER status as provided by Nik-Zainal et al.\
final.HER2: HER2 status as provided by Nik-Zainel et al.\
CNA of genes of interest (ERBB2_CNA through RPS6KB1_CNA): binary indication if sample has amplification in gene of interest defined as CN > 4\
snvs: number of somatic SNVs\
ERBB2_CNA_5copies: binary indicating amplification of ERBB2 considering gene amplified if CN > 5\

#### References
S. Nik-Zainal et al., Landscape of somatic mutations in 560 breast cancer whole-genome sequences. Nature 534, 47–54 (2016).

### METABRIC
Summary table for METABRIC (replication) cohort.\
**file:** metabric_megatable.txt

#### Columns
sample: sample name\
GEB in genes of interest (ERBB2 through MYC): GEB in each gene of interest\
OS_MONTHS: overall survival in months\
AGE_AT_DIAGNOSIS: age at diagnosis\
INTCLUST: IC10 subtype as defined by Curtis et al.\
ER_IHC: ER status by IHC as defined by Curtis et al.\
HER2_SNP6: HER2 status by SNP6 array as defined by Curtis et al.\
CLAUDIN_SUBTYPE: PAM50 subtype as defined by Curtis et al.\
PC1-6: principal components 1 through 6\
CNA of genes of interest (RPS6KB1_CNA through MYC_CNA): binary indication if sample has amplification in gene of interest defined as CN > 4\
node: number of positive lymph nodes as defined by Rueda et al.\
size: tumor size as defined by Rueda et al.\
grade: tumor grade as defined by Rueda et al.\
stage: tumor stage as defined by Rueda et al.\
age: age at diagnosis as defined by Rueda et al.\
pga: percent genome altered

#### References
C. Curtis et al., The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature 486, 346–352 (2012).
O. M. Rueda et al., Dynamics of breast-cancer relapse reveal late-recurring ER-positive genomic subgroups. Nature 567, 399–404 (2019).

### DCIS
Summary table for DCIS cohort.\
**file:** dcis_megatable.txt

#### Columns
sample: sample name\
GEB in genes of interest (ERBB2 through MYC): GEB in each gene of interest\
EA_PC1, EA_PC2, AA_PC1, NA_PC1, NA_PC2, NA_PC3, CO_PC1, CO_PC2: ancestry principal components as calculated by SNPWEIGHTS\
Cohort: cohort name\
HER2_RNA: HER2 status defined by mRNA\
ER_RNA: ER status defined by mRNA\
IC11_RNA: IC10 subtype defined by iC10\
PAM50: pam50 subtype\
Mo_FU: follow-up time (in months)\
Diagnostic group: Diagnostic group indicating if individual experienced a DCIS (DCIS_with_DCIS_recurrence) or invasive breast cancer (DCIS_with_IBC_recurrence) recurrence\
CNA of genes of interest (ERBB2_CNA through RPS6KB1_CNA): binary indication if sample has amplification in gene of interest defined as CN > 4\
grade: lesion grade as defined by Strand et al.\

#### References
S. H. Strand et al., Molecular classification and biomarkers of clinical outcome in breast ductal carcinoma in situ: Analysis of TBCRC 038 and RAHBT cohorts. Cancer Cell 40, 1521-1536.e7 (2022).


