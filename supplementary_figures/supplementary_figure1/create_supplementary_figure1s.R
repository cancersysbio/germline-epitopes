### CREATE SUPPLEMENTARY FIGURE 1S ################################################################
# calculate r2 for common vs rare variants

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date() 
### MAIN ##########################################################################################
# read in plot data
plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'controls', 'rare_common_r2_estimates.txt'),
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
