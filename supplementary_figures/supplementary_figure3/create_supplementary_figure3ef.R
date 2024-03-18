### CREATE SUPPLEMENTARY FIGURE 3EF ###############################################################
# create supplementary figure 3ef 
# scatterplot of allele frequencies in primary vs metastatic tumors
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### SUPPLEMENTARY FIGURE 3E #######################################################################
# read in tcga and hartwig frequencies 
tcga <- read.delim(
    file.path(main_repo_path, 'data', 'mafs', 'tcga_maf_gnomad.txt'), 
    as.is = TRUE
    )
hartwig <- read.delim(
    file.path(main_repo_path, 'data', 'mafs','hartwig_maf_gnomad.txt'), 
    as.is = TRUE
    )

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

# lift tcga to hg19
tcga$chr <- paste0('chr', tcga$chr) 
bed <- tcga[order(tcga$chr, tcga$start),]
bed_gr <- makeGRangesFromDataFrame(bed)
hg19 <- liftOver(bed_gr, ch)
hg19 <- as.data.frame(hg19)
mapping <- cbind(
    bed[,c('chr','start')],
    hg19[,c('seqnames','start')]
    )
mapping <- mapping[,-3]
colnames(mapping) <- c('chr','start','hg19')
tcga <- merge(tcga, mapping, by = c('chr','start'), all.x = TRUE)
# create snp id with hg19 coordinates
tcga$snp <- paste(gsub('chr','', tcga$chr), tcga$hg19, sep = '_')

th_plot_data <- merge(tcga[,c('snp','gene','maf')], hartwig[,c('snp','gene','maf')], by = c('snp','gene'), all = TRUE)
colnames(th_plot_data) <- c('gene','snp','tcga','hartwig')

th_plot_data[is.na(th_plot_data)] <- 0

# create scatterplot 
create.scatterplot(
        hartwig ~ tcga,
        data = th_plot_data,
        ylimits = c(0,0.5),
        xlimits = c(0,0.5),
        yat = seq(0,0.4,0.1),
        xat = seq(0,0.4,0.1),
        add.xyline = TRUE,
        ylab.label = 'Hartwig',
        xlab.label = 'TCGA',
        filename = paste0(date, '_supplementary_figure3e.png'),
        legend = list(
             inside = list(
                     fun = draw.key,
                     args = list(
                         key = get.corr.key(
                             x = th_plot_data$tcga,
                             y = th_plot_data$hartwig,
                             label.items = c('pearson','pearson.p'),
                             alpha.background = 0,
                             key.cex = 1.5
                             )
                         ),
                     x = 0.04,
                     y = 0.95,
                     corner = c(0,1)
                     )
                 ),
        resolution = 300
        )


#### SUPPLEMENTARY FIGURE 3F ######################################################################
# read in icgc 
icgc <- read.delim(
    file.path(main_repo_path, 'data', 'mafs','icgc_maf_gnomad.txt'), 
    as.is = TRUE
    )

# merge 
ih_plot_data <- merge(icgc[,c('snp','gene','maf')], hartwig[,c('snp','gene','maf')], by = c('gene','snp'), all = TRUE)
colnames(ih_plot_data) <- c('gene','snp','icgc','hartwig')

ih_plot_data[is.na(ih_plot_data)] <- 0

# plot_data 
create.scatterplot(
    hartwig ~ icgc, 
    data = ih_plot_data,
    ylab.label = 'Hartwig',
    xlab.label = 'ICGC',
    ylimits = c(0,0.5),
    xlimits = c(0,0.5),
    yat = seq(0,0.4,0.1),
    xat = seq(0,0.4,0.1),
    filename = paste0(date, '_supplementary_figure3f.png'),
    add.xyline = TRUE,
    legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = ih_plot_data$icgc,
                         y = ih_plot_data$hartwig,
                         label.items = c('pearson','pearson.p'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.04,
                 y = 0.95,
                 corner = c(0,1)
                 )
             ),
    resolution = 300
    )

