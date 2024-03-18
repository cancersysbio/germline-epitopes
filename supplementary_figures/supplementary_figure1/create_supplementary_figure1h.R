#### SUPPLEMENTARY FIGURE 1H ######################################################################
# create supplementary figure 1h

### PREAMBLE ######################################################################################
library(GenomicRanges)
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST SUBTYPE ASSOCIATION ######################################################################
run_subtype_associations <- function(dtf, subtype, gene = NULL, hla_correct = FALSE) {
        genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1','NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32'),
                IC10 = c('FOXC1','MIA','MELK','KNTC2')
                )

        if (!is.null(gene)) {
                dtf$bds <- sign(dtf[,gene])
                if (subtype == 'Her2') {
                        dtf$subtype <- (dtf$pam50 == 'Her2')*1
                } else {
                        dtf$subtype <- (dtf[,paste0(gene, '_CNA')] == 1 & dtf$ER == 1)*1
                        }
        } else if (subtype == 'IC10') {
                dtf$subtype <- (dtf$pam50 == 'Basal')*1
                dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
        } else {
                dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
                dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER == 1)*1
        }

        
        # run association 
        if (!hla_correct) {
                fit <- glm(
                        subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic,
                        data = dtf,
                        family = 'binomial'
                        )
        } else {
                fit <- glm(
                        subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic + HLA_CN,
                        data = dtf,
                        family = 'binomial'
                        )
                }
        ci <- confint(fit)

        out <- data.frame(
                subtype = subtype,
                gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
                coef = coef(fit)[['bds']],
                p = summary(fit)$coefficients['bds',4],
                se = summary(fit)$coefficients['bds',2],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                number_subtype = sum(dtf$subtype)
                )
        return(out)
}

### FIND HLA CN ###################################################################################
find_hla_cn <- function(cna) {
        # find CN of 6p21
        cna_gr <- makeGRangesFromDataFrame(
                        cna,
                        keep.extra.columns=TRUE,
                        seqnames.field = 'Chromosome',
                        start.field = 'Start',
                        end.field = 'End'
                        )
        gene_gr <- GRanges(
                        seqnames = 6,
                        # defining HLA region as start of HLA-A - 1Mbp to end of HLA-B + 1Mbp
                        ranges = IRanges(start = 28942532, end = 32357179)
                        )
        overlap <- subsetByOverlaps(
                        cna_gr,
                        gene_gr
                        )

        hla_cna <- aggregate(mcols(overlap)[,'Modal_Total_CN'], list(mcols(overlap)[,'Sample']), median)
        colnames(hla_cna) <- c('Sample','HLA_CN')
        hla_cna$sample <- substr(hla_cna$Sample, 1, 12)
        # remove met sample 
        hla_cna <- hla_cna[!hla_cna$Sample %in% c('TCGA-AC-A6IX-06','TCGA-E2-A15A-06','TCGA-E2-A15E-06','TCGA-E2-A15K-06'),]
        return(hla_cna)
}


### MAIN ##########################################################################################
# read in tcga data 
tcga <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
        as.is = TRUE
        )
# read in cna
cnafile <- file.path(main_repo_path, 'data', 'auxiliary_data', 'tcga_cna_segments_chr6.txt')
cna <- read.delim(
        cnafile,
        as.is = TRUE
        )

# remove met sample 
hla_cna <- find_hla_cn(cna)

# all hla cna to tcga 
tcga_hla <- merge(tcga, hla_cna[,c('sample','HLA_CN')], by = 'sample', all.x = TRUE)

# test subtype association 
plot_data_subtype_hla <- rbind(
        run_subtype_associations(tcga_hla, gene = 'ERBB2', subtype = 'Her2', hla_correct = TRUE),
        run_subtype_associations(tcga_hla, subtype = 'IC1', hla_correct = TRUE),
        run_subtype_associations(tcga_hla, subtype = 'IC2', hla_correct = TRUE),
        run_subtype_associations(tcga_hla, subtype = 'IC9', hla_correct = TRUE)
        )
plot_data_subtype_hla$experiment <- 'hla'

# read in cn calls as 2+ploidy
cn_2ploidy <- read.delim(
        file.path(main_repo_path, 'data', 'auxiliary_data', 'tcga_cn_2ploidy.txt'),
        as.is = TRUE
        )

# all hla cna to tcga 
tcga_2ploidy <- merge(tcga[,!grepl('CNA', colnames(tcga))], cn_2ploidy, by = 'sample', all.x = TRUE)

# test subtype association 
plot_data_subtype_2ploidy <- rbind(
        run_subtype_associations(tcga_2ploidy, gene = 'ERBB2', subtype = 'Her2'),
        run_subtype_associations(tcga_2ploidy, subtype = 'IC1'),
        run_subtype_associations(tcga_2ploidy, subtype = 'IC2'),
        run_subtype_associations(tcga_2ploidy, subtype = 'IC9')
        )
plot_data_subtype_2ploidy$experiment <- 'ploidy'

### SUBTYPE PLOT ####
plot_data_subtype <- rbind(plot_data_subtype_hla, plot_data_subtype_2ploidy)
plot_data_subtype <- plot_data_subtype[order(plot_data_subtype$subtype, plot_data_subtype$experiment),]
plot_data_subtype$index <- 1:nrow(plot_data_subtype)

cov <- list(
        rect = list(
                col = 'transparent',
                fill = rep(c('mediumpurple1', 'mediumpurple4'), 5)
                )
        );

cov.grob <- covariates.grob(
        covariates = cov,
        ord = c(1:nrow(plot_data_subtype)),
        side = 'right',
        size = 1
        );

cov.legend <- list(
        legend = list(
                colours =  c('mediumpurple1', 'mediumpurple4'),
                labels = c('Ploidy','HLA CN'),
                title = 'Controling for',
                border = 'transparent'
                )
        );

cov.legend.grob <- legend.grob(
        legends = cov.legend
        );

create.scatterplot(
        index ~ coef,
        data = plot_data_subtype,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        filename = paste0(date, '_supplementary_figure1h.png'),
        xat = log(c(0.2, 0.5, 1, 2, 5)),
        xaxis.lab = c('0.2','0.5','1.0','2.0','5.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = seq(1.5, 8, 2),
        ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(0,0), x = 0.01, y = 0.05)
                ),
        add.rectangle = TRUE,
        xleft.rectangle = -2.5,
        ybottom.rectangle = c(0, seq(4.5, nrow(plot_data_subtype), 4)),
        xright.rectangle = 2.5,
        ytop.rectangle = c(seq(2.5, nrow(plot_data_subtype), 4)),
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        add.text = TRUE,
        text.labels = plot_data_subtype$gene,
        text.y = 1:nrow(plot_data_subtype),
        text.x = rep(1.65,nrow(plot_data_subtype)),
        text.cex = 0.8,
        x.error.right = plot_data_subtype$u95-plot_data_subtype$coef,
        x.error.left = plot_data_subtype$coef-plot_data_subtype$l95,
        width = 8.5,
        height = 3,
        top.padding = 2,
        left.padding = 5,
        resolution = 300
        )

