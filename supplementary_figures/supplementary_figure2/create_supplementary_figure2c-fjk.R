### CREATE SUPPLEMENTARY FIGURE 2C-F, J-K #########################################################
# create supplementary figure 1c 
# scatterplot of maf in ICGC, METABRIC, DCIS
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}


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

### SUPPLEMENTARY FIGURE 2C #######################################################################
# read in maf in tcga 
imafs <- read.delim(
	file.path(main_repo_path, 'data', 'mafs', 'icgc_maf_gnomad.txt'),
	as.is = TRUE
	)

# create maf scatterplot
create_maf_scatterplot(mafs = imafs, filename = paste0(date, '_supplementary_figure2c.png'),
    xlab.label = 'ICGC')

### SUPPLEMENTARY FIGURE 2D #######################################################################
# read in maf in tcga 
mmafs <- read.delim(
    file.path(main_repo_path, 'data', 'mafs', 'metabric_maf_gnomad.txt'),
    as.is = TRUE
    )

# create maf scatterplot
create_maf_scatterplot(mafs = mmafs, filename = paste0(date, '_supplementary_figure2d.png'),
    xlab.label = 'METABRIC')

### SUPPLEMENTARY FIGURE 2J #######################################################################
# read in maf in tcga 
dmafs <- read.delim(
    file.path(main_repo_path, 'data', 'mafs', 'dcis_maf_gnomad.txt'),
    as.is = TRUE
    )

# create maf scatterplot
create_maf_scatterplot(mafs = dmafs, filename = paste0(date, '_supplementary_figure2j.png'),
    xlab.label = 'DCIS')

#### SUPPLEMENTARY FIGURE 2E ######################################################################
# read in tcga icgc frequencies 
icgc_freq <- read.delim(
    file.path(main_repo_path, 'data', 'mafs','icgc_hla_frequencies.txt'), 
    as.is = TRUE
    )
# create scatterplot
create_hla_scatterplot(icgc_freq, filename = paste0(date, '_supplementary_figure2e.png'))


#### SUPPLEMENTARY FIGURE 2F ######################################################################
# read in metabric hla frequencies 
metabric_freq <- read.delim(
    file.path(main_repo_path, 'data', 'mafs','metabric_hla_frequencies.txt'), 
    as.is = TRUE
    )
# create scatterplot
create_hla_scatterplot(metabric_freq, filename = paste0(date, '_supplementary_figure2f.png'))

#### SUPPLEMENTARY FIGURE 2K ######################################################################
# read in dcis hla frequencies 
dcis_freq <- read.delim(
    file.path(main_repo_path, 'data', 'mafs','dcis_hla_frequencies.txt'), 
    as.is = TRUE
    )
# create scatterplot
create_hla_scatterplot(dcis_freq, filename = paste0(date, '_supplementary_figure2k.png'))
