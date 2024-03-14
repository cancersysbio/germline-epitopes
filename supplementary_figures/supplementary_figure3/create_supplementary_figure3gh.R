### CREATE SUPPLEMENTARY FIGURE 3GH ###############################################################
# create supplementary figure 3gh 
# scatterplot of hla frequencies in primary vs metastatic tumors
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

date <- Sys.Date()
### SUPPLEMENTARY FIGURE 3G #######################################################################
# read in tcga and hartwig frequencies 
tcga <- read.delim('tcga_hla_frequencies.txt', as.is = TRUE)
hartwig <- read.delim('hartwig_hla_frequencies.txt', as.is = TRUE)

th_plot_data <- merge(tcga, hartwig, by = 'allele', all = TRUE)
colnames(th_plot_data) <- c('allele','tcga','hartwig')

th_plot_data[is.na(th_plot_data)] <- 0

# create scatterplot 
create.scatterplot(
        hartwig ~ tcga,
        data = th_plot_data,
        ylimits = c(0,0.30),
        xlimits = c(0,0.30),
        yat = seq(0,0.3,0.1),
        xat = seq(0,0.3,0.1),
        add.xyline = TRUE,
        filename = paste0(date, '_tcga_vs_hartwig_HLA_frequency_scatterplot.png'),
        ylab.label = 'Hartwig',
        xlab.label = 'TCGA',
        legend = list(
                inside = list(
                        fun = draw.key,
                        args = list(
                             key = get.corr.key(
                                 x = th_plot_data$tcga,
                                 y = th_plot_data$hartwig,
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

#### SUPPLEMENTARY FIGURE 3H ######################################################################
# read in icgc 
icgc <- read.delim('icgc_hla_frequencies.txt', as.is = TRUE)

# merge 
ih_plot_data <- merge(icgc, hartwig, by = 'allele', all = TRUE)
colnames(ih_plot_data) <- c('allele','icgc','hartwig')

ih_plot_data[is.na(ih_plot_data)] <- 0

# create scatterplot 
create.scatterplot(
        hartwig ~ icgc,
        data = ih_plot_data,
        ylimits = c(0,0.30),
        xlimits = c(0,0.30),
        yat = seq(0,0.3,0.1),
        xat = seq(0,0.3,0.1),
        add.xyline = TRUE,
        filename = paste0(date, '_icgc_vs_hartwig_HLA_frequency_scatterplot.png'),
        ylab.label = 'Hartwig',
        xlab.label = 'ICGC',
        legend = list(
                inside = list(
                        fun = draw.key,
                        args = list(
                             key = get.corr.key(
                                 x = ih_plot_data$icgc,
                                 y = ih_plot_data$hartwig,
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

