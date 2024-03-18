### CREATE SUPPLEMENTARY FIGURE 1A #################################################################
# create barplot of HLA accuracy

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in allele accuracy
alleles <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data', 'tcga_cookhla_hla_alleles_imputation_accuracy.txt'),
	as.is = TRUE
	)
alleles$index <- 1:nrow(alleles)

# create barplot
create.barplot(
	accuracy ~ index,
	data = alleles,
	abline.h = 0.8,
	ylimits = c(0,1.05),
	yat = seq(0,1,0.2),
	xat = c(12.5,46,78.5),
	xaxis.lab = c('HLA-A','HLA-B','HLA-C'),
	add.rectangle = TRUE,
	xleft.rectangle = c(0,67.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(25.5,91),
	ytop.rectangle = 1.05,
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	xlab.label = 'HLA Alleles',
	ylab.label = 'Concordance',
	filename = paste0(date, '_supplementary_figure1a.png'),
	width = 8,
	height = 5,
	resolution = 300
	)
 
