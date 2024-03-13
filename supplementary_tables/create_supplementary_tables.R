### CREATE SUPPLEMENTARY TABLES ###################################################################

### CALCULATE SUBTYPE EPITOPE BURDEN ##############################################################
calculate_subtype_epitope_burden <- function(dtf, subtype) {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32')
		)
	return(rowSums(sign(dtf[,genes[[subtype]]])))
}

### ASSIGN IC CNA SUBTYPE #########################################################################
assign_ic_cna_subtype <- function(dtf, subtype, cohort = 'TCGA') {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32')
		)
	if (cohort == 'TCGA') {
		return((dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER == 1)*1)
	} else if (cohort == 'ICGC') {
		return((dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$final.ER == 'positive')*1)
	} else if (cohort == 'Hartwig') {
		return((dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive'))*1)	
	} else if (cohort == 'DCIS') {
		return((dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER_RNA == '+')*1)
	} else if (cohort == 'METABRIC') {
		return((dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER_IHC == 'pos')*1)
	}
}

### TCGA ##########################################################################################
# read in tcga 
tcga <- read.delim(
	'tcga_megatable.txt',
	as.is = TRUE
	)

# create tcga st 
tcga_st <- data.frame(
	sample = tcga$sample,
	HER2_IHC = tcga$HER2.newly.derived,
	ER = tcga$ER,
	HER2_pam50 = (tcga$pam50 == 'Her2')*1,
	HER2_EB = sign(tcga$ERBB2),
	IC1_CNA = assign_ic_cna_subtype(tcga, 'IC1'),
	IC1_EB = calculate_subtype_epitope_burden(tcga, 'IC1'),
	IC2_CNA = assign_ic_cna_subtype(tcga, 'IC2'),
	IC2_EB = calculate_subtype_epitope_burden(tcga, 'IC2'),
	IC6_CNA = assign_ic_cna_subtype(tcga, 'IC6'),
	IC6_EB = calculate_subtype_epitope_burden(tcga, 'IC6'),
	IC9_CNA = assign_ic_cna_subtype(tcga, 'IC9'),
	IC9_EB = calculate_subtype_epitope_burden(tcga, 'IC9')
	)
# add pcs 
tcga_st <- cbind(tcga_st, tcga[,c('PC1','PC2','PC3','PC4','PC5','PC6','somatic')])

# write to table 
write.table(
	tcga_st,
	file = paste0(date, '_supplementary_table1_tcga.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### ICGC ##########################################################################################
# read in icgc 
icgc <- read.delim(
	'icgc_megatable.txt',
	as.is = TRUE
	)

# create icgc st 
icgc_st <- data.frame(
	sample = icgc$sample_name,
	HER2 = (icgc$final.HER2 == 'positive')*1,
	ER = (icgc$final.ER == 'positive')*1,
	HER2_ER = (icgc$final.HER2 == 'positive' & icgc$final.ER == 'negative')*1,
	HER2_EB = sign(icgc$ERBB2),
	IC1_CNA = assign_ic_cna_subtype(icgc, 'IC1', cohort = 'ICGC'),
	IC1_EB = calculate_subtype_epitope_burden(icgc, 'IC1'),
	IC2_CNA = assign_ic_cna_subtype(icgc, 'IC2', cohort = 'ICGC'),
	IC2_EB = calculate_subtype_epitope_burden(icgc, 'IC2'),
	IC6_CNA = assign_ic_cna_subtype(icgc, 'IC6', cohort = 'ICGC'),
	IC6_EB = calculate_subtype_epitope_burden(icgc, 'IC6'),
	IC9_CNA = assign_ic_cna_subtype(icgc, 'IC9', cohort = 'ICGC'),
	IC9_EB = calculate_subtype_epitope_burden(icgc, 'IC9')
	)
# add pcs 
icgc_st <- cbind(icgc_st, icgc[,c('PC1','PC2','PC3','PC4','PC5','PC6','snvs')])
colnames(icgc_st) <- gsub('snvs','somatic_snv', colnames(icgc_st))

# write to table 
write.table(
	icgc_st,
	file = paste0(date, '_supplementary_table1_icgc.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### DCIS ##########################################################################################
# read in dcis 
dcis <- read.delim(
	'dcis_megatable.txt',
	as.is = TRUE
	)

# create dcis st 
dcis_st <- data.frame(
	sample = dcis$sample,
	HER2_pam50 = (dcis$PAM50 == 'Her2')*1,
	HER2_EB = sign(dcis$ERBB2),
	IC1_CNA = assign_ic_cna_subtype(dcis, 'IC1', cohort = 'DCIS'),
	IC1_EB = calculate_subtype_epitope_burden(dcis, 'IC1'),
	IC2_CNA = assign_ic_cna_subtype(dcis, 'IC2', cohort = 'DCIS'),
	IC2_EB = calculate_subtype_epitope_burden(dcis, 'IC2')
	)
# add pcs 
dcis_st <- cbind(dcis_st, dcis[,c('EA_PC1','EA_PC2','AA_PC1','NA_PC1','NA_PC2','NA_PC3','CO_PC1','CO_PC2')])

# write to table 
write.table(
	dcis_st,
	file = paste0(date, '_supplementary_table1_dcis.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### METABRIC ######################################################################################
# read in metabric
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)

# create dcis st 
metabric_st <- data.frame(
	sample = metabric$sample,
	HER2_pam50 = (metabric$CLAUDIN_SUBTYPE == 'Her2')*1,
	HER2_EB = sign(metabric$ERBB2),
	IC1_CNA = assign_ic_cna_subtype(metabric, 'IC1', cohort = 'METABRIC'),
	IC1_EB = calculate_subtype_epitope_burden(metabric, 'IC1'),
	IC2_CNA = assign_ic_cna_subtype(metabric, 'IC2', cohort = 'METABRIC'),
	IC2_EB = calculate_subtype_epitope_burden(metabric, 'IC2'),
	IC9_CNA = assign_ic_cna_subtype(metabric, 'IC9', cohort = 'METABRIC'),
	IC9_EB = calculate_subtype_epitope_burden(metabric, 'IC9')
	)
metabric_st <- cbind(metabric_st, metabric[,c('PC1','PC2','PC3','PC4','PC5','PC6')])

# write to table 
write.table(
	metabric_st,
	file = paste0(date, '_supplementary_table1_metabric.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### CLINICAL ###
# consider only her2 patients
meta_her2 <- metabric[metabric$CLAUDIN_SUBTYPE == 'Her2',]
meta_her2$event <- (meta_her2$OS_STATUS == 'DECEASED')*1
meta_her2$bds <- (meta_her2$ERBB2 > 0)*1
# testing 5-year relapse 
meta_her2[meta_her2$OS_MONTHS > 60,'event'] <- 0
meta_her2[meta_her2$OS_MONTHS > 60,'OS_MONTHS'] <- 60
# reformat 
meta_her2_st <- meta_her2[,c('sample','bds','OS_MONTHS','event','PC1','PC2','AGE_AT_DIAGNOSIS','pga','INTCLUST','grade','stage')]
colnames(meta_her2_st) <- c('sample','ERBB2_EB','OS_MONTHS','OS_STATUS','PC1','PC2','AGE_AT_DIAGNOSIS','PGA','INTCLUST','GRADE','STAGE')

# write to table 
write.table(
	meta_her2_st,
	file = paste0(date, '_supplementary_table1_metabric_her2.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# ER+ high risk tumors
meta_prog <- metabric[which(rowSums(metabric[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & metabric$ER_IHC == 'pos'),]
meta_prog$burden <- sign(colSums(rbind(
	rowSums(sign(meta_prog[,c('MYC','SQLE','FBXO32')]))*meta_prog$MYC_CNA,
	rowSums(sign(meta_prog[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_prog$RPS6KB1_CNA,
	rowSums(sign(meta_prog[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_prog$RSF1_CNA,
	rowSums(sign(meta_prog[,c('ZNF703','FGFR1','LETM2')]))*meta_prog$ZNF703_CNA
	)))
# add pga because calculate burden across multiple amplicons
meta_prog$event <- (meta_prog$OS_STATUS == 'DECEASED')*1
meta_prog$burden <- (meta_prog$burden > 0)*1
# considering 5-year relapse
meta_prog[which(meta_prog$OS_MONTHS > 60),'event'] <- 0
meta_prog[which(meta_prog$OS_MONTHS > 60),'OS_MONTHS'] <- 60
# reformat 
meta_hr_st <- meta_prog[,c('sample','burden','OS_MONTHS','event','PC1','PC2','AGE_AT_DIAGNOSIS','pga','INTCLUST','grade','stage')]
colnames(meta_hr_st) <- c('sample','IC_SUM_EB','OS_MONTHS','OS_STATUS','PC1','PC2','AGE_AT_DIAGNOSIS','PGA','INTCLUST','GRADE','STAGE')

# write to table 
write.table(
	meta_hr_st,
	file = paste0(date, '_supplementary_table1_metabric_er.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### HARTWIG #######################################################################################
# read in hartwig 
hartwig <- read.delim(
	'hartwig_megatable.txt',
	as.is = TRUE
	)
# create hartwig st 
hartwig_st <- data.frame(
	sample = hartwig$sampleId,
	HER2_ER = (hartwig$hormone == 'ER-negative/HER2-positive')*1,
	HER2_EB = sign(hartwig$ERBB2),
	IC1_CNA = assign_ic_cna_subtype(hartwig, 'IC1', cohort = 'Hartwig'),
	IC1_EB = calculate_subtype_epitope_burden(hartwig, 'IC1'),
	IC2_CNA = assign_ic_cna_subtype(hartwig, 'IC2', cohort = 'Hartwig'),
	IC2_EB = calculate_subtype_epitope_burden(hartwig, 'IC2'),
	IC9_CNA = assign_ic_cna_subtype(hartwig, 'IC9', cohort = 'Hartwig'),
	IC9_EB = calculate_subtype_epitope_burden(hartwig, 'IC9')
	)
# add pcs 
hartwig_st <- cbind(hartwig_st, hartwig[,c('PC1','PC2','PC3','PC4','PC5','PC6','snvs')])
colnames(hartwig_st) <- gsub('snvs','somatic_snv', colnames(hartwig_st))

# write to table 
write.table(
	hartwig_st,
	file = paste0(date, '_supplementary_table1_hartwig.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
