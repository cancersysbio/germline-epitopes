### CREATE SUPPLEMENTARY FIGURE 3C #################################################################
# create supplementary figure 3c
# create PCA plot of Hartwig with reference

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### SUPPLEMENTARY FIGURE 3C #######################################################################
# read in pcs
plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ancestry_inference', 'hartwig_pca_with_reference.txt'),
	as.is = TRUE
	)

# create barplot
pc1_threshold <- max(plot_data[plot_data$ancestry == 'EUR','PC1'])+0.2
pc2_threshold <- min(plot_data[plot_data$ancestry == 'EUR','PC2'])-0.2
create.scatterplot(
	PC2 ~ PC1,
	data = plot_data,
	filename = paste0(date, '_supplementary_figure3c.png'),
	groups = plot_data$ancestry,
	abline.h = pc2_threshold,
	abline.v = pc1_threshold,
	abline.lty = 2,
	col = default.colours(6),
	key = list(
             text = list(
                 lab = c('AFR','AME','E ASN','EUR','Hartwig','OCE'),
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

