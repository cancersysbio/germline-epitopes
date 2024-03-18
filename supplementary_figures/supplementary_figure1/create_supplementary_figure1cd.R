### CREATE SUPPLEMENTARY FIGURE 1C&D ##############################################################
# create supplementary figure 1c 
# scatterplot of maf in TCGA 
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()

### SUPPLEMENTARY FIGURE 1C #######################################################################
# read in maf in tcga 
allmafs <- read.delim(
	file.path(main_repo_path, 'data', 'mafs', 'tcga_maf_gnomad.txt'),
	as.is = TRUE
	)

create.scatterplot(
	gnomad ~ maf,
	data = allmafs,
	filename = paste0(date, '_supplementary_figure1c.png'),
	add.xyline = TRUE,
	ylab.label = 'Gnomad',
	xlab.label = 'TCGA',
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = allmafs$gnomad,
                         y = allmafs$maf,
                         label.items = c('pearson','pearson.p'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1)
                 )
             ),
	resolution = 300
	)

#### SUPPLEMENTARY FIGURE 1D ######################################################################
create_hla_scatterplot <- function(freq) {
        pop <- list()
        for (allele in c('A','B','C')) {
                pop[[allele]] <- read.csv(
                        file.path(main_repo_path, 'data', 'mafs', paste0("population_hla_", tolower(allele), "_frequencies.csv")),
                        as.is = TRUE
                        )
        }
        pop <- do.call(rbind, pop)
        # create plot data
        plot_data <- merge(freq, pop, by = 'allele', all.x = TRUE)
        stats <- cor.test(plot_data$freq, plot_data$frequency, method = 'spearman')
        pvalue_sci <- scientific.notation(stats$p.value, type = 'list')
        # create scatterplot 
        create.scatterplot(
                freq ~ frequency,
                data = plot_data,
                ylimits = c(0,0.30),
                xlimits = c(0,0.30),
                yat = seq(0,0.3,0.1),
                xat = seq(0,0.3,0.1),
                add.xyline = TRUE,
                filename = paste0(date, '_supplementary_figure1d.png'),
                ylab.label = 'Cohort Frequencies',
                xlab.label = 'Population Frequencies',
                legend = list(
                     inside = list(
                         fun = draw.key,
                         args = list(
                             key = get.corr.key(
                                 x = plot_data$freq,
                                 y = plot_data$frequency,
                                 label.items = c('pearson','pearson.p'),
                                 alpha.background = 0,
                                 key.cex = 1.2
                                 )
                             ),
                         x = 0.01,
                         y = 0.99,
                         corner = c(0,1)
                         )
                     ),
                resolution = 300
                )
}

# read in tcga hla frequencies 
hla_freq <- read.delim(
    file.path(main_repo_path, 'data', 'mafs', 'tcga_hla_frequencies.txt'), 
    as.is = TRUE)
# create scatterplot
create_hla_scatterplot(hla_freq)
