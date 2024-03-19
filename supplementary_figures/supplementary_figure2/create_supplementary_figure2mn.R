### CREATE SUPPLEMENTARY FIGURE 2MN ################################################################
# create read depth barplots
# create supplementary figure 2MN

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(plyr)

main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### CREATE READS BARPLOT ##########################################################################
create_reads_barplot <- function(readsdf, chr, pos, filename, thresholds = c(0.2, 0.8)) {
	# reformat plot data
	readsdf <- readsdf[order(as.numeric(readsdf$tumor_ratio)),]
	readsdf$index <- 1:nrow(readsdf)
	col <- rep('black', nrow(readsdf))
	col[which(readsdf$tumor_ratio < thresholds[1])] <- 'mediumpurple4'
	col[which(readsdf$tumor_ratio > thresholds[2])] <- 'darkolivegreen4'
	# create barplot
	bar1 <- create.barplot(
			tumor_ratio ~ index,
			data = readsdf,
			ylimits = c(0,1.1),
			main = 'Tumor',
			ylab.label = 'Fraction Alt Allele',
			xaxis.lab = rep('', nrow(readsdf)),
			abline.h = c(thresholds[1], 0.5, thresholds[2]),
			xlab.label = '',
			col = col,
			abline.lty = c(2,1,2),
			text.above.bars = list(labels = rowSums(readsdf[,c('tref','talt')]),
		                     padding = 0.05,
		                     bar.locations = 1:nrow(readsdf),
		                     rotation = 90
		                     ),
			legend = list(
	             inside = list(
	                 fun = draw.key,
	                 args = list(
	                     key = list(
	                         points = list(
	                             col = 'black',
	                             pch = 22,
	                             cex = 3,
	                             fill = c('mediumpurple4', 'darkolivegreen4')
	                             ),
	                         text = list(
	                             lab = c('Low Alt Freq','High Alt Freq')
	                             ),
	                         padding.text = 5,
	                         cex = 1.5
	                         )
	                     ),
	                     # Positioning legend on plot
	                     x = 0.01,
	                     y = 0.99
	                 )
	             ),
			resolution = 300
			)
	bar2 <- create.barplot(
			normal_ratio ~ index,
			data = readsdf,
			ylimits = c(0,1),
			main = 'Normal',
			ylab.label = 'Fraction Alt Allele',
			xaxis.lab = rep('', nrow(readsdf)),
			abline.h = 0.5,
			col = col,
			text.above.bars = list(labels = rowSums(readsdf[,c('nref','nalt')]),
		                     padding = 0.05,
		                     bar.locations = 1:nrow(readsdf),
		                     rotation = 90
		                     ),
			xlab.label = 'Sample',
			resolution = 300
			)
	create.multipanelplot(
			list(bar1, bar2),
			filename = filename,
			resolution = 300,
			width = 12
			)
	}

### MAIN ##########################################################################################
# create read frequency barplot for ERBB2 variant
erbb2_readsdf <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data', 'erbb2_snp_read_depths.txt'),
	as.is = TRUE
	)
create_reads_barplot(
	readsdf = erbb2_readsdf,
	chr = 17, 
	pos = 39727784,
	filename = paste0(date, '_supplementary_figure2m.png'), 
	thresholds = c(0.4, 0.6)
	)

# create read frequency barplot for TUBD1 variant
tubd1_readsdf <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data', 'tubd1_snp_read_depths.txt'),
	as.is = TRUE
	)
create_reads_barplot(
	readsdf = tubd1_readsdf,
	chr = 17, 
	pos = 59886176, 
	filename = paste0(date, '_supplementary_figure2n.png'),
	thresholds = c(0.4, 0.6)
	)

