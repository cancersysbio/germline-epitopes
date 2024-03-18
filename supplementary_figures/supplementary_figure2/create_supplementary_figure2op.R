### SUPPLEMENTARY FIGURE 2OP ######################################################################
# create supplemenatary figure 2op
# test allele specific amplification 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(plyr)

main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

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
	file.path(main_repo_path, 'data', 'auxiliary_data', 'erbb2_snp_read_depths.txt'),
	as.is = TRUE
	)
# read in antigens generated from run_antigen_garnish.R
erbb2_antigens <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data','erbb2_snp_peptides.txt'), 
	as.is = TRUE
	)
erbb2_antigens <- reformat_antigens(erbb2_readsdf, antigens = erbb2_antigens, 
	gene = 'ERBB2', chr = 17, pos = 39727784, thresholds = c(0.4, 0.6))

# calculate median per sample
erbb2_antigens <- aggregate(erbb2_antigens$rank_diff, erbb2_antigens[,c('sample','type')], median)
colnames(erbb2_antigens) <- c('sample','type','rank_diff')

# test difference in median
es <- median(erbb2_antigens[erbb2_antigens$type == 'high','rank_diff'])-median(erbb2_antigens[erbb2_antigens$type == 'low','rank_diff'])
stats <- wilcox.test(
	erbb2_antigens[erbb2_antigens$type == 'high','rank_diff'],
	erbb2_antigens[erbb2_antigens$type == 'low','rank_dff'],
	alternative = 'greater'
	)
pvalue <- scientific.notation(stats$p.value, digits = 2, type = 'list');
main <- 'ERBB2|rs1058808'
filename <- paste0(date, '_supplementary_figure2o.png')

# create barplot
create.boxplot(
	rank_diff ~ type,
	data = erbb2_antigens,
	add.stripplot = TRUE,
	xlab.label = 'Allele Amplified',
	xaxis.lab = c('Alt','Ref'),
	ylab.label = 'Median Differential binding\n(Alt binding strength)',
	ylimits = c(-0.1,1),
	yat = seq(-0.1,1,0.2),
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
	resolution = 300
	)

### TUBD1 ##########################################################################################
# identify samples with low and high alt allele counts
tubd1_readsdf <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data','tubd1_snp_read_depths.txt'),
	as.is = TRUE
	)
# read in antigens generated from run_antigen_garnish.R
tubd1_antigens <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data','tubd1_snp_peptides.txt'), 
	as.is = TRUE
	)
tubd1_antigens <- reformat_antigens(tubd1_readsdf, antigens = tubd1_antigens, 
	gene = 'TUBD1', chr = 17, pos = 59886176, thresholds = c(0.4, 0.6))

# calculate median per sample
tubd1_antigens <- aggregate(tubd1_antigens$rank_diff, tubd1_antigens[,c('sample','type')], median)
colnames(tubd1_antigens) <- c('sample','type','rank_diff')

# test difference in median
es <- median(tubd1_antigens[tubd1_antigens$type == 'high','rank_diff'])-median(tubd1_antigens[tubd1_antigens$type == 'low','rank_diff'])
stats <- wilcox.test(
	tubd1_antigens[tubd1_antigens$type == 'high','rank_diff'],
	tubd1_antigens[tubd1_antigens$type == 'low','rank_diff'],
	alternative = 'less'
	)
pvalue <- scientific.notation(stats$p.value, digits = 2, type = 'list');
main <- 'TUBD1|rs1292053'
filename <- paste0(date, '_supplementary_figure2p.png')

# create barplot
create.boxplot(
	rank_diff ~ type,
	data = tubd1_antigens,
	add.stripplot = TRUE,
	xlab.label = 'Allele Amplified',
	xaxis.lab = c('Alt','Ref'),
	ylab.label = 'Median Differential binding\n(Alt binding strength)',
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
	resolution = 300
	)