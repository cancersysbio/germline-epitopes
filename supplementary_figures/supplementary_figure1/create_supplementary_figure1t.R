### CREATE SUPPLEMENTARY FIGURE 1T ################################################################
# create plot of benign vs pathogenic variants
# create supplementary figure 1T 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)
library(vcfR)

date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-d', '--dir', type = 'integer', help = 'base directory');

args <- parser$parse_args();
### FILTER MISSENSE VARIANTS ######################################################################
filter_missense_variants <- function(vcf, id, gene) {
        # extract info
        info <- INFO2df(vcf)
        missense <- sapply(info$ANN, function(i) strsplit(i, '\\|')[[1]][2])
        benign <- sapply(info$ANN, function(i) strsplit(i, '\\|')[[1]][3])
        out <- data.frame(
                snp = paste(vcf@fix[,'CHROM'], vcf@fix[,'POS'], sep = '_'),
                sample = id,
                gene = gene,
                missense = missense,
                benign = benign
                )
        rownames(out) <- 1:nrow(out)
        return(out)
        }

### FIND SNPS #####################################################################################
find_snps <- function(samples, dir) {
        # read in epiptope predictions
        gts <- list()
        for (gene in c('ERBB2','MYC','SQLE','FBXO32','RSF1','CCND1','PAK1','NARS2',
                'RPS6KB1','TUBD1','DHX40','BCAS3')) {
                for (id in samples) {
                        # find file flag
                        file <- file.path(dir, gene, id, paste0(id, '_', gene, '_snpeff_missense'))
                        vcf <- read.vcfR(paste0(file, '.vcf'))
                        if (nrow(vcf) > 0) {
                                gts[[paste0(id, gene)]] <- filter_missense_variants(vcf, id = id, gene = gene)
                        }
                }
        }
        gts <- do.call(rbind, gts)
        return(gts)
}

find_clinvar_annotations <- function(gts, dir) {
        genes <- c('ERBB2',
                'RPS6KB1','TUBD1','DHX40','BCAS3',
                'RSF1','CCND1','PAK1','NARS2',
                'MYC','SQLE','FBXO32'
                )
        clinvar <- do.call(rbind, sapply(
               genes,
                function(x) {
                        tmp <- read.delim(file.path(dir, 'clinvar', paste0(x, '_clinvar_result.txt')), as.is = TRUE)
                        tmp <- tmp[,c('GRCh38Chromosome','GRCh38Location','Clinical.significance..Last.reviewed.')]
                        tmp$snp <- paste0(
                                'chr', tmp$GRCh38Chromosome,
                                '_', sapply(tmp$GRCh38Location, function(y) strsplit(y, ' -')[[1]][1])
                                )
                        tmp <- tmp[tmp$snp %in% gts$snp,]
                        if (nrow(tmp) > 0) {
                                tmp$gene <- x
                                return(tmp)
                                }
                        },
                        simplify = FALSE
                ))
        # merge clinvar annotations with snps
        gts <- merge(gts, clinvar[,c('snp','Clinical.significance..Last.reviewed.')], by = 'snp', all.x = TRUE)
        gts$Clinical.significance..Last.reviewed. <- sapply(
                gts$Clinical.significance..Last.reviewed.,
                function(x) {
                        strsplit(x, '\\(')[[1]][1]
                }
                )
        # find pathogenic variants
        gts$pathogenic <- apply(
                gts, 
                1,
                function(x) {
                        ifelse(x['benign'] == 'HIGH' | x['Clinical.significance..Last.reviewed.'] %in% c('Likely pathogenic'),'yes', 'no') 
                        })

        gts_snp <- unique(gts[,c('snp','gene','pathogenic')])
        gts_snp$gene <- factor(gts_snp$gene, levels = genes)
        plot_data <- as.data.frame(table(gts_snp[,c('pathogenic','gene')]))
        return(plot_data)
}


### MAIN ##########################################################################################
# read in all samples that have epitope predictions
samples <- read.delim(
                'tcga_megatable.txt',
                as.is = TRUE
                )
 genes <- c('ERBB2','RPS6KB1','TUBD1','DHX40','BCAS3',
        'RSF1','CCND1','PAK1','NARS2','MYC','SQLE','FBXO32'
        )

# read in plot data which was generated from
#gts <- find_snps(samples = samples$sample, dir = args$dir)
#plot_data <- find_clinvar_annotations(gts = gts, dir = args$dir)
plot_data <- read.delim(
        'pathogenic_benign_variants.txt',
        as.is = TRUE
        )

plot_data$index <- NA
for (i in 1:length(genes)) {
        plot_data[plot_data$gene == genes[i],'index'] <- letters[i]
}

create.barplot(
        Freq ~ index,
        groups = plot_data$pathogenic,
        data = plot_data,
        xaxis.lab = genes,
        xaxis.cex = 1,
        xaxis.rot = 45,
        width = 9,
        ylab.label = 'Number Unique Variants',
        xlab.label = 'Gene',
        ylimits = c(0,40),
        yat = seq(0,40,10),
        col = c('darkgrey','black'),
        add.rectangle = TRUE,
        xleft.rectangle = c(5.5,9.5),
        ybottom.rectangle = 0,
        xright.rectangle = c(1.5,13.5),
        ytop.rectangle = 40,
        col.rectangle = 'grey',
        alpha.rectangle = 0.5,
        add.text = TRUE,
        text.labels = c('HER2+','IC1','IC2','IC9'),
        text.x = c(1,3.5,7.5,11),
        text.y = 38,
        xaxis.fontface = 4,
        filename = paste0(date, '_pathogenic_benign_variants_barplot.pdf'),
        legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('darkgrey', 'black')
                             ),
                         text = list(
                             lab = c('Benign','Pathogenic')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     ),
                     # Positioning legend on plot
                     x = 0.8,
                     y = 0.85
                 )
             ),
        resolution = 300
        )




