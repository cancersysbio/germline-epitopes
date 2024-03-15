### CREATE FIGURE 1C ##############################################################################
# test presence of HLAs that can present IISAVVGIL with HER2 subtype
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

date <- Sys.Date()
### CREATE CONTINGENCY MULTIPLOT ##################################################################
create_gp2_barplot <- function(df, filename) {
	cont_table <- table(df[,c('hla','subtype')])
	bplot_data <- data.frame(
		ratio = c(
			cont_table[1,2]/cont_table[1,1],
			cont_table[2,2]/cont_table[2,1]
			),
		bds = c(0,1)
		)		
	fit <- glm(subtype ~ hla + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = df, family = 'binomial')
	or <- exp(coef(fit)[['hla']])
	pvalue_sci <- scientific.notation(summary(fit)$coefficients['hla',4], digits = 2,type = 'list');
	main1 <- as.expression(substitute(
	                                base, 
	                                list(base = paste0('OR = ', round(or, digits = 2)))
	                                ))
	main2 <- as.expression(substitute(
	                                base *' x '* 10^exponent, 
	                                list(base = paste0('P = ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
	                                ))
	        
	create.barplot(
		ratio ~ bds,
		ylimits = c(0,0.32),
		yat = seq(0,0.3,0.1),
		filename = filename,
		xaxis.lab = c('Low','High'),
		xlab.label = 'HLAs present GP2/E75',
		ylab.label = 'Ratio of HER2+/HER2-\n',
		data = bplot_data,
		add.text = TRUE,
		text.labels = c(main1, main2),
		text.x = c(2.08, 2.2),
		text.y = c(0.31,0.29),
		text.cex = 1.2,
		resolution = 300 
		)

}

### COUNT NUMBER OF BINDING ALLELES ###############################################################
count_number_binding_alleles <- function(dtf, alleles) {
	# bin patients by alleles that bind 
	hlasbinddf <- do.call(rbind, sapply(
		unique(dtf$sample),
		function(x) {
			tmp <- dtf[dtf$sample == x,]
			data.frame(
				sample = x,
				hla = sum(tmp$hla %in% alleles)
				)
			},
		simplify = FALSE
		))
	return(hlasbinddf)
}

### GET TCGA HLAs #################################################################################
get_tcga_hlas <- function(samples) {
	# read in hlas 
	hlas <- read.delim(
		'../../data/auxiliary_data/Shukla_Wu_Getz_Polysolver_HLA_Types_2015.tsv',
		header = FALSE,
		as.is = TRUE
		)
	hlas <- hlas[hlas$V1 %in% paste0('BRCA-', samples),]
	hlas$V1 <- gsub('BRCA-','',hlas$V1)
	hlas_rf <- data.frame(sample = rep(hlas$V1, 6), hla = unlist(hlas[,-1]))
	return(hlas_rf)
}

### MAIN ##################################################################################
# read in samples 
tcga <- read.delim(
	'../../data/cohort_megatables/tcga_megatable.txt',
	as.is = TRUE
	)
samples <- tcga$sample
# read in hlas 
gp2_e75_hlas <- read.delim(
	'../../data/auxiliary_data/gp2_e75_hlas.txt',
	as.is = TRUE,
	header = FALSE
	)
gp2_e75_hlas <- gp2_e75_hlas$V1

# read in hlas 
tcga_hlas <- get_tcga_hlas(samples)

# count number of binding alleles 
tcga_hlas_both <- count_number_binding_alleles(tcga_hlas, gp2_e75_hlas)

# annotate
tcga_hlas_both <- merge(tcga_hlas_both, tcga, by = 'sample')
# set subtype
tcga_hlas_both$hla <- (tcga_hlas_both$hla >= median(tcga_hlas_both$hla))*1
tcga_hlas_both$subtype <- (tcga_hlas_both$HER2.newly.derived == 'Positive')*1

# create barplot
create_gp2_barplot(tcga_hlas_both, filename = paste0(date, '_gp2_barplot.png'))


