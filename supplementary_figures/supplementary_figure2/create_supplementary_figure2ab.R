### CREATE SUPPLEMENTARY FIGURE 2AB #################################################################
# create supplementary figure 2ab
# create PCA plot of ICGC and METABRIC with reference

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### SUPPLEMENTARY FIGURE 2A #######################################################################
# read in pcs
iplot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ancestry_inference', 'icgc_pca_with_reference.txt'),
	as.is = TRUE
	)

# create barplot
ipc1_threshold <- max(iplot_data[iplot_data$ancestry == 'EUR','PC1'])+0.2
ipc2_threshold <- min(iplot_data[iplot_data$ancestry == 'EUR','PC2'])-0.2
create.scatterplot(
	PC2 ~ PC1,
	data = iplot_data,
	filename = paste0(date, '_supplementary_figure2a.png'),
	groups = iplot_data$ancestry,
	abline.h = ipc2_threshold,
	abline.v = ipc1_threshold,
	abline.lty = 2,
	col = default.colours(6),
	key = list(
             text = list(
                 lab = c('AFR','AME','E ASN','EUR','ICGC','OCE'),
                 cex = 1,
                 col = 'black'
                 ),
             points = list(
                 pch = 19,
                 col = default.colours(6),
                 cex = 1
                 ),
             x = 0.01,
             y = 0.01,
             corner = c(0,0),
             padding.text = 2
             ),
	resolution = 300
	)

### CREATE SUPPLEMENTARY FIGURE 2B ################################################################
# read in pcs
mplot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ancestry_inference', 'metabric_pca_with_reference.txt'),
	as.is = TRUE
	)

# create barplot
mpc1_threshold <- max(mplot_data[mplot_data$ancestry == 'EUR','PC1'])+0.2
mpc2_threshold <- min(mplot_data[mplot_data$ancestry == 'EUR','PC2'])-0.2
create.scatterplot(
	PC2 ~ PC1,
	data = mplot_data,
	filename = paste0(date, '_supplementary_figure2b.png'),
	groups = mplot_data$ancestry,
	abline.h = mpc2_threshold,
	abline.v = mpc1_threshold,
	abline.lty = 2,
	col = default.colours(6),
	key = list(
             text = list(
                 lab = c('AFR','AME','E ASN','EUR','METABRIC','OCE'),
                 cex = 1,
                 col = 'black'
                 ),
             points = list(
                 pch = 19,
                 col = default.colours(6),
                 cex = 1
                 ),
             x = 0.01,
             y = 0.01,
             corner = c(0,0),
             padding.text = 2
             ),
	resolution = 300
	)