### FIGURE 2CD ####################################################################################
# plot allele specific binding 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)
library(plyr)
library(vcfR)
library(bedr)

# Set the main path for repo
main_repo_path <- "~/stanford/projects/BC_GermlineEpitopes/code/germline-epitopes"

date <- Sys.Date()
### REFORMAT ANTIGENS ##############################################################################
reformat_antigens <- function(readsdf, antigens, gene, chr, pos, thresholds = c(0.2, 0.8)) {
	# set base direcotry and snp
	# split into high and low
	althigh <- readsdf[readsdf$tumor_ratio > thresholds[2],]
	altlow <- readsdf[readsdf$tumor_ratio < thresholds[1],]
	samples <- c(as.character(althigh$sample), as.character(altlow$sample))
	# keep only antigens mapping to position of interest
	antigens_snp <- antigens[which(antigens$pos == pos),c('sample','rank_diff')]
	# find any samples without peptides
	nopep_samples <- samples[!samples %in% antigens_snp$sample] 
	if (length(nopep_samples) > 0) {
		nopeptides <- data.frame(sample = nopep_samples, rank_diff = 0)
		antigens_snp <- rbind(
			antigens_snp,
			nopeptides
			)
	}
	antigens_snp$type <- NA
	antigens_snp[antigens_snp$sample %in% althigh$sample,'type'] <- 'high'
	antigens_snp[antigens_snp$sample %in% altlow$sample,'type'] <- 'low'
	colnames(antigens_snp) <- c('sample','rank_diff','type')

	return(antigens_snp)
}

### ERBB2 #########################################################################################
# identify samples with low and high alt allele counts
erbb2_readsdf <- read.delim(
  file.path(main_repo_path,'data','auxiliary_data','erbb2_snp_read_depths.txt'),
	as.is = TRUE
	)

# read in antigens generated from run_antigen_garnish.R
erbb2_antigens <- read.delim(
  file.path(main_repo_path,'data','auxiliary_data','erbb2_snp_peptides.txt'), 
  as.is = TRUE)

erbb2_antigens <- reformat_antigens(erbb2_readsdf, antigens = erbb2_antigens, 
	gene = 'ERBB2', chr = 17, pos = 39727784, thresholds = c(0.4, 0.6))

# test difference in median
es <- median(erbb2_antigens[erbb2_antigens$type == 'high','rank_diff'])-median(erbb2_antigens[erbb2_antigens$type == 'low','rank_diff'])
stats <- wilcox.test(
	erbb2_antigens[erbb2_antigens$type == 'high','rank_diff'],
	erbb2_antigens[erbb2_antigens$type == 'low','rank_dff'],
	conf.int = TRUE
	)
pvalue <- scientific.notation(stats$p.value, digits = 2, type = 'list');
main <- 'ERBB2|rs1058808'
filename <- paste0(date, '_ERBB2_rs1058808_antigen_boxplot.pdf')

# create barplot
create.boxplot(
	rank_diff ~ type,
	data = erbb2_antigens,
	add.stripplot = TRUE,
	xlab.label = 'Allele Amplified',
	xaxis.lab = c('Alt','Ref'),
	ylab.label = 'Differential binding\n(Alt binding strength)',
	ylimits = c(-1,4),
	yat = seq(-1,4,1),
	main = main,
	main.cex = 2,
	filename = filename,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         text = list(
                             lab = c(
                             	paste('ES =', round(es, digits = 4)),
                             	as.expression(substitute(
			                        base *' x '* 10^exponent, 
			                        list(base = paste0('P =', pvalue[[1]]), exponent = pvalue[[2]])
			                        ))
                             	)
                             ),
                         cex = 1.5
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1),
                 draw = FALSE
                 )
             ),
	height = 7,
	width = 6,
	resolution = 300
	)

### TUBD1 ##########################################################################################
# identify samples with low and high alt allele counts
tubd1_readsdf <- read.delim(
  file.path(main_repo_path,'data','auxiliary_data','tubd1_snp_read_depths.txt'),
	as.is = TRUE
	)

# read in antigens generated from run_antigen_garnish.R
tubd1_antigens <- read.delim(
  file.path(main_repo_path,'data','auxiliary_data','tubd1_snp_peptides.txt'),
  as.is = TRUE)

tubd1_antigens <- reformat_antigens(tubd1_readsdf, antigens = tubd1_antigens, 
	gene = 'TUBD1', chr = 17, pos = 59886176, thresholds = c(0.4, 0.6))

# test difference in median
es <- median(tubd1_antigens[tubd1_antigens$type == 'high','rank_diff'])-median(tubd1_antigens[tubd1_antigens$type == 'low','rank_diff'])
stats <- wilcox.test(
	tubd1_antigens[tubd1_antigens$type == 'high','rank_diff'],
	tubd1_antigens[tubd1_antigens$type == 'low','rank_diff'],
	conf.int = TRUE
	)
pvalue <- scientific.notation(stats$p.value, digits = 2, type = 'list');
main <- 'TUBD1|rs1292053'
filename <- paste0(date, '_TUBD1_rs1292053_antigen_boxplot.pdf')

# create barplot
create.boxplot(
  formula = rank_diff ~ type,
	data = tubd1_antigens,
	add.stripplot = TRUE,
	xlab.label = 'Allele Amplified',
	xaxis.lab = c('Alt','Ref'),
	ylab.label = 'Differential binding\n(Alt binding strength)',
	yat = seq(-0.5,2,0.5),
	ylimits = c(-1,2.2),
	main = main,
	main.cex = 2,
	filename = filename,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         text = list(
                             lab = c(
                             	paste('ES =', round(es, digits = 4)),
                             	as.expression(substitute(
			                        base *' x '* 10^exponent, 
			                        list(base = paste0('P =', pvalue[[1]]), exponent = pvalue[[2]])
			                        ))
                             	)
                             ),
                         cex = 1.5
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1),
                 draw = FALSE
                 )
             ),
	width = 6,
	height = 7,
	resolution = 300
	)