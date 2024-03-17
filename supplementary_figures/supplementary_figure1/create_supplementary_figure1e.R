### SUPPLEMENTARY FIGURE 1E #########################################################################
# create supplementary figure 1e
# barplot of distribution of epitopes across regions

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN #######################################################################################
# read in summary data
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'), 
	as.is = TRUE
	)
         
# set genes to evaluate per subtype 
genes <- list(
	IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
	IC2 = c('RSF1','CCND1','PAK1','NARS2'),
	IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
	IC9 = c('MYC','SQLE','FBXO32'),
	HER2 = 'ERBB2',
	IC10 = c('FOXC1','MELK','MIA')
    )
# iterate over subtypes 
plot_data <- list()
for (subtype in names(genes)) {
	if (length(genes[[subtype]]) > 1) {
		bds <- sign(rowSums(sign(tcga[,genes[[subtype]]])))
	} else {
		bds <- sign(tcga[,genes[[subtype]]])
		}
	tmp <- as.data.frame(table(bds))
	tmp$subtype <- subtype
	plot_data[[subtype]] <- tmp
}
plot_data <- do.call(rbind, plot_data)
plot_data[plot_data$subtype == 'IC10','subtype'] <- 'TNBC'

# create barplot 
create.barplot(
	Freq ~ subtype,
	groups = bds,
	data = plot_data,
	filename = paste0(date, '_bds_per_regions_barplot.png'),
	col = default.colours(5, palette.type = 'spiral.afternoon')[c(1,3,5)],
	ylimits = c(0,650),
	yat = seq(0,650,100),
	xaxis.lab = c('ERBB2\n(HER2+)','17q23\n(IC1)','11q13\n(IC2)','8p12\n(IC6)','8q24\n(IC9)','Various\n(TNBC)'),
	ylab.label = 'Number of samples',
	xlab.label = 'Region (Subtype)',
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					points = list(
						col = 'black',
						pch = 22,
						cex = 2,
						fill = default.colours(5, palette.type = 'spiral.afternoon')[c(1,3,5)]
						),
					text = list(
						lab = c('<WT','WT','>WT')
						),
					padding.text = 3,
					cex = 1.5
					)
				),
			x = 0.05,
			y = 0.95
			)
		),
	width = 7,
	resolution = 300
	)
