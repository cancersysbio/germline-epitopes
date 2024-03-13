### CREATE SUPPLEMENTARY FIGURES L-O, U ###########################################################
# create supplementary figures l-o and u
# create boxplots of mRNA abundance

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

date <- Sys.Date()
### CREATE IC BOXPLOT ##############################################################################
create_ic_boxplot <- function(ic, icgenes, ylimits = NULL, yat = TRUE) {
     # create boxplots for each subtype 
        icdata <- as.data.frame(t(rna_tumor[which(rna_tumor$id %in% icgenes),-1]))
        colnames(icdata) <- rna_tumor$id[which(rna_tumor$id %in% icgenes)]
        icdata <- icdata[rownames(icdata) %in% BRCA_annotation[BRCA_annotation$ic11 == tolower(ic),'ID'],]
        icdata <- as.data.frame(apply(icdata, 2, as.numeric))
        # calculate median 
        medians <- apply(icdata, 2, median, na.rm = TRUE)
        medians <- medians[order(-medians)]
        # reformat
        icdata <- gather(icdata, value = 'rna', key = 'gene')
        # order genes by median
        icdata$index <- NA
        for (j in 1:length(medians)) {
                icdata[icdata$gene == names(medians)[j],'index'] <- j
        }
        icdata$index <- factor(icdata$index, levels = 1:length(medians))
        if (ic != 'ic9') {
                rectright <- 4.5
        } else {
                rectright <- 3.5
        }

        # create boxplot 
        create.boxplot(
                rna ~ index,
                data = icdata,
                main = toupper(ic),
                ylimits = ylimits,
                yat = yat,
                xaxis.rot = 45,
                xaxis.lab = names(medians),
                add.stripplot = TRUE,
                ylab.label = 'mRNA abundance',
                xlab.label = 'Genes',
                add.rectangle = TRUE,
                xleft.rectangle = 0,
                ybottom.rectangle = 0,
                xright.rectangle = rectright,
                ytop.rectangle = ylimits[2],
                xaxis.fontface = 4,
                col.rectangle = 'grey',
                alpha.rectangle = 0.5,
                filename = paste0(date, '_TCGA_', toupper(ic), '_genes_rna_boxplot.png'),
                resolution = 300
                )   
}

### MAIN ##########################################################################################
# read in rna 
# file downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas
rna <- read.delim(
        'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
        as.is = TRUE
        )
rna_genes <- sapply(rna$gene_id, function(x) strsplit(x, '\\|')[[1]][1])
tissues <- sapply(colnames(rna), function(x) strsplit(x, '\\.')[[1]][4])
# remove normal tissue 
normal <- which(tissues %in% c('11A','11B','11C'))
rna_tumor <- rna[,-normal]
colnames(rna_tumor) <- gsub('.', '-', colnames(rna_tumor), fixed = TRUE)
rna_tumor$id <- rna_genes

# read in TCGA annotations
load('genefuAnnotations.RData')
BRCA_annotation$donor <- substr(BRCA_annotation$ID, 1, 12)
# keep only primary
BRCA_annotation <- BRCA_annotation[order(BRCA_annotation$ID),]
BRCA_annotation <- BRCA_annotation[!duplicated(BRCA_annotation$donor),]

# create list of genes and IC subtypes 
ic1 <- c('RPS6KB1','TUBD1','DHX40','BCAS3','BRIP1','PRR11','TBX2','PPM1E','HSF5','CA4','C17orf64','TBC1D3P2')
ic6 <- c('ZNF703','FGFR1','LETM2','EIF4EBP1','STAR')
ic9 <- c('MYC','SQLE','FBXO32','ADCY8')
ic2 <- c('CCND1','RSF1','PAK1',"NARS2",'FGF3','OMP')

ne <- c('KRT34','KRT71','KRT74','KRT82')

ic1ne <- c('CSH1','CSH2','CSHL1','GH1')
ic2ne <- c('FGF4','CABP2','DEFB108B','OR2AT4')
ic9ne <- c('CYP11B1','GML', 'CYP11B2')
ic5ne <- c('HNF1B')


create_ic_boxplot(ic = 'ic1', icgenes = ic1, ylimits = c(0,10000), yat = seq(0,10000,2000))
create_ic_boxplot(ic = 'ic6', icgenes = ic6, ylimits = c(0,60000), yat = seq(0,60000,20000))
create_ic_boxplot(ic = 'ic9', icgenes = ic9, ylimits = c(0,15000), yat = seq(0,15000,5000))
create_ic_boxplot(ic = 'ic2', icgenes = ic2, ylimits = c(0,150000), yat = seq(0,150000,50000))


all_genes <- c(
        'DHX40','TUBD1','RPS6KB1','BCAS3', 
        ic1ne,
        'CCND1','RSF1','NARS2','PAK1', 
        ic2ne,
        'MYC','SQLE','FBXO32', 
        ic9ne,
        'ERBB2',
        ic5ne,
        ne
        )

plot_data <- as.data.frame(t(rna_tumor[rna_tumor$id %in% all_genes,-ncol(rna_tumor)]))
colnames(plot_data) <- sapply(as.character(unlist(plot_data[1,])), function(x) strsplit(x, '\\|')[[1]][1])
plot_data <- plot_data[-1,]
plot_data <- gather(plot_data, key = 'gene', value = 'rna')
plot_data$index <- NA
for (i in 1:length(all_genes)) {
        plot_data[plot_data$gene == all_genes[i],'index'] <- i
}
plot_data$index <- factor(plot_data$index)
plot_data$rna <- as.numeric(plot_data$rna)

# create boxplot 
create.boxplot(
                log10(rna) ~ index,
                data = plot_data,
                #ylimits = ylimits,
                #yat = yat,
                xaxis.rot = 45,
                xaxis.cex = 1.2,
                xaxis.lab = all_genes,
                add.text = TRUE,
                text.labels = c('IC1','IC2','IC9','HER2','Keratins'),
                text.x = c(4.5,12.5,19.5,23.5,26.5),
                text.y = rep(6, 5),
                add.stripplot = TRUE,
                ylab.label = expression('log'[10]*' mRNA abundance'),
                xlab.label = 'Genes',
                add.rectangle = TRUE,
                xleft.rectangle = c(0,8.5,16.5,22.5,24.5),
                ybottom.rectangle = -100,
                xright.rectangle = c(4.5,12.5,19.5,23.5,31.5),
                ytop.rectangle = 1000,
                col.rectangle = 'grey',
                alpha.rectangle = 0.5,
                xaxis.fontface = 4,
                filename = paste0(date, '_all_genes_rna_boxplot.png'),
                width = 10,
                height = 8,
                resolution = 300
                ) 

### CREATE SUPPLEMENTARY FIGURE 1W ################################################################
# create plot_data
plot_data <- rna_tumor[rna_tumor$id %in% c('FOXC1','MIA','MELK'),]
rownames(plot_data) <- plot_data$id
plot_data <- as.data.frame(t(plot_data[,-1]))
plot_data$sample <- substr(rownames(plot_data), 1, 12)
plot_data <- merge(plot_data, tcga[,c('sample','Triple.Negative.Status')], by = 'sample')
plot_data$FOXC1 <- as.numeric(as.character(plot_data$FOXC1))
plot_data$MIA <- as.numeric(as.character(plot_data$MIA))
plot_data$MELK <- as.numeric(as.character(plot_data$MELK))

# reformat plot data
plot_data <- gather(
        plot_data[,c('sample','Triple.Negative.Status','FOXC1','MIA','MELK')],
        value = 'mrna',
        key = 'gene',
        -Triple.Negative.Status,
        -sample
        )

plot_data$logmrna <- log10(plot_data$mrna)
create.boxplot(
        logmrna ~ Triple.Negative.Status | gene,
        data = plot_data[which(plot_data$gene %in% c('FOXC1','MELK','MIA') & plot_data$Triple.Negative.Status != 'N/A'),],
        add.stripplot = TRUE,
        xlab.label = 'Triple Negative Status',
        ylab.label = 'mRNA Abundance',
        yat = seq(0,6,2),
        ylimits = c(0,6),
        yaxis.lab = c(expression(0), expression(10^2), expression(10^4), expression(10^6)),
        yaxis.cex = 1.3,
        filename = paste0(date, '_genes_rna_tnbc_boxplot.pdf'),
        resolution = 300
        )