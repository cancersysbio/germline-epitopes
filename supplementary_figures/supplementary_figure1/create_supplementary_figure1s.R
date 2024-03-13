### CREATE SUPPLEMENTARY FIGURE 1S ################################################################
# calculate r2 for common vs rare variants

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)
library(plyr)
library(vcfR)

date <- Sys.Date() 
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-d', '--dir', type = 'integer', help = 'base directory');

args <- parser$parse_args();
#### READ IN AND FORMAT TCGA SNP ##################################################################
read_in_and_format_tcga_snps <- function(gene, basedir) {
	# read in all samples that have epitope predictions
	samples <- read.delim(
                'tcga_megatable.txt',
                as.is = TRUE
                )
	# read in epiptope predictions
	gts <- list()
	for (id in samples$sample) {
		# find file flag
		file <- file.path(basedir, gene, id, paste0(id, '_', gene, '_snpeff_missense'))
		vcf <- read.vcfR(paste0(file, '.vcf'))
		if (nrow(vcf) > 0) {
			gt <- extract.gt(vcf)
			gt_rf <- data.frame(
				snp = rownames(gt),
				geno = gt[,1]
				)
			colnames(gt_rf) <- c('snp', id)
			gts[[id]] <- gt_rf
		} else {
			tmp <- data.frame(snp = '', geno = '0/0')
			colnames(tmp) <- c('snp', id)
			gts[[id]] <- tmp
		}
	}
	gts_df <- join_all(gts, by = 'snp', type = 'full')
	gts_df <- gts_df[gts_df$snp != '',]
	gts_df <- apply(gts_df, 2, as.character)
	gts_df[is.na(gts_df)] <- '0/0'
	rownames(gts_df) <- gts_df[,1]
	gts_df <- as.data.frame(t(gts_df[,-1]))
	gts_df$sample <- rownames(gts_df)

	# read in ancestry
	ancestry <- read.csv(
	        'yuan_tcga_breast_ancestry.csv',
	        as.is = TRUE,
	        header = FALSE
	        )
	colnames(ancestry) <- c('sample','cancer','self-reported','ancestry','symbol')
	ancestry <- ancestry[which(ancestry$symbol == 'EA'),]

	# only keep european ancestry 
	gts_df <- gts_df[rownames(gts_df) %in% ancestry$sample,]
	return(gts_df)
}

### FIND COMMON RARE VARIANTS #####################################################################
find_common_rare_variants <- function(snps, maf) {
	common <- maf[maf$maf > 0.01,'snp']
	rare <- maf[maf$maf < 0.01,'snp']
	# find number common vs rare 
	snp_count <- do.call(rbind, apply(
		snps,
		1,
		function(x) {
			x <- x[x != '0/0']
			data.frame(
				rare = sum(names(x) %in% rare),
				common = sum(names(x) %in% common)
				)
			}))
	snp_count$sample <- rownames(snp_count)
	return(snp_count)
}

### TEST MAF GEB ##################################################################################
test_maf_geb <- function(subtype, basedir) {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32')
		)
	# read in tcga
	tcga <- read.delim(
		'tcga_megatable', 
		as.is = TRUE
		)
	if (subtype == 'HER2') {
		snps <- read_in_and_format_tcga_snps('ERBB2', basedir = basedir)
		# read in maf
		maf <- read.delim(
			file.path('mafs', paste0('2022-06-02_tcga_ERBB2_maf.txt')),
			as.is = TRUE
			)
		# calculate number of rare vs common 
		snp_count <- find_common_rare_variants(snps, maf)
		# merge with geb 
		tcga$bds <- sign(tcga$ERBB2)
	} else {
		snp_count_tmp <- list()
		for (gene in genes[[subtype]][3:4]) {
			snps <- read_in_and_format_tcga_snps(gene, basedir = basedir)
			# read in maf
			maf <- read.delim(
				file.path('mafs', paste0('2022-06-02_tcga_', gene, '_maf.txt')),
				as.is = TRUE
				)
			# calculate number of rare vs common 
			tmp <- find_common_rare_variants(snps, maf)
			colnames(tmp) <- c(paste0(gene, '_rare'), paste0(gene, '_common'), 'sample')
			snp_count_tmp[[gene]] <- tmp
			}
		snp_count <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample", all.x = TRUE),
        		snp_count_tmp)
		snp_count$common <- rowSums(snp_count[,grep('_common', colnames(snp_count))])
		snp_count$rare <- rowSums(snp_count[,grep('_rare', colnames(snp_count))])
		# merge with geb 
		tcga$bds <- rowSums(sign(tcga[,genes[[subtype]]]))
	}
	# merge 
	dtf <- merge(snp_count, tcga[,c('sample','bds','PC1','PC2','PC3','PC4','PC5','PC6')], by = 'sample')
	dtf <- dtf[!is.na(dtf$bds),]
	# calculate correlation
	r_common <- cor(dtf$bds, dtf$common)^2
	r_rare <- cor(dtf$bds, dtf$rare)^2
	out <- data.frame(
		subtype = subtype,
		common = r_common,
		rare = r_rare
		)
	return(out)
}

### MAIN ##########################################################################################
# calculate correlation between GEB and number of rare or common variants
her2 <- test_maf_geb('HER2', basedir = args$dir)
ic1 <- test_maf_geb('IC1', basedir = args$dir)
ic2 <- test_maf_geb('IC2', basedir = args$dir)
ic9 <- test_maf_geb('IC9', basedir = args$dir)

plot_data <- rbind(her2, ic1, ic2, ic6, ic9)

write.table(
	plot_data,
	file = 'rare_common_r2_estimates.txt',
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# read in plot data
plot_data <- read.delim(
	'rare_common_r2_estimates.txt',
	as.is = TRUE
	)

plot_data <- gather(
	plot_data,
	value = 'r2',
	key = 'type',
	-subtype
	)

create.barplot(
	r2 ~ subtype,
	groups = plot_data$type,
	data = plot_data,
	stack = TRUE,
	col = default.colours(12)[c(10,8)],
	ylab.label = expression('R'^2),
	xlab.label = 'Subtype',
	ylimits = c(0,0.61),
	yat = seq(0,1,0.2),
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 2,
                             # reverse order to match stacked bar order
                             fill = default.colours(12)[c(10,8)]
                             ),
                         text = list(
                             # reverse order to match stacked bar order
                             lab = c('Common','Rare')
                             ),
                         padding.text = 3,
                         cex = 1
                         )
                     ),
                 x = 0.01,
                 y = 0.99
                 )
             ),
	filename = paste0(date, '_common_rare_r2_barplot.png'),
	resolution = 300
	)
