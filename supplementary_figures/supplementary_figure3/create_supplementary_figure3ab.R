### CREATE SUPPLEMENTARY FIGURE 3AB ###############################################################
# create supplementary figure 3ab 
# scatterplot of maf in Hartwig
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

date <- Sys.Date()

### CREATE MAF SCATTERPLOT ########################################################################
create_maf_scatterplot <- function(mafs, filename, xlab.label) {
       create.scatterplot(
            gnomad ~ maf,
            data = mafs,
            filename = filename,
            add.xyline = TRUE,
            ylab.label = 'Gnomad',
            xlab.label = xlab.label,
            ylimits = c(0,0.5),
            xlimits = c(0, 0.5),
            yat = seq(0, 0.4, 0.1),
            xat = seq(0, 0.4, 0.1),
            legend = list(
                     inside = list(
                         fun = draw.key,
                         args = list(
                             key = get.corr.key(
                                 x = mafs$gnomad,
                                 y = mafs$maf,
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
}

### CREATE HLA SCATTERPLOT ########################################################################
create_hla_scatterplot <- function(freq, filename) {
        pop <- list()
        for (allele in c('A','B','C')) {
                pop[[allele]] <- read.csv(
                        paste0("population_hla_", tolower(allele), "_frequencies.csv"),
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
                filename = filename,
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

### SUPPLEMENTARY FIGURE 3A #######################################################################
# read in maf in tcga 
hmafs <- read.delim(
	'hartwig_maf_gnomad.txt',
	as.is = TRUE
	)

# create maf scatterplot
create_maf_scatterplot(mafs = hmafs, filename = paste0(date, '_gnomad_hartwig_scatterplot.png'),
    xlab.label = 'Hartwig')


#### SUPPLEMENTARY FIGURE 3B ######################################################################
# read in tcga hla frequencies 
hartwig_freq <- read.delim('hartwig_hla_frequencies.txt', as.is = TRUE)
# create scatterplot
create_hla_scatterplot(hartwig_freq, filename = paste0(date, '_hartwig_hla_frequencies_scatterplot.png'))


