### CREATE SUPPLEMENTARY FIGURE 1K ################################################################
# compare true associations to associations run with scrambled HLA alleles
# supplementary figure 1k

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(plyr)

date <- Sys.Date()
### TEST ASSOCIATION ##############################################################################
run_subtype_associations <- function(dtf, subtype, gene = NULL) {
        genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40', 'BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1', 'NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32')
                )

        if (!is.null(gene)) {
                dtf$bds <- sign(dtf[,gene])
                if (subtype == 'Her2') {
                        dtf$subtype <- (dtf$pam50 == 'Her2')*1
                } else {
                        dtf$subtype <- (dtf[,paste0(gene, '_CNA')] == 1 & dtf$ER == 1)*1
                        }
        } else {
                dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
                dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER == 1)*1
        }
        dtf <- dtf[!is.na(dtf$somatic),]

        # run association 
        fit <- glm(
                subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic,
                data = dtf,
                family = 'binomial'
                )
        ci <- confint(fit)

        out <- data.frame(
                subtype = subtype,
                gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
                coef = coef(fit)[['bds']],
                p = summary(fit)$coefficients['bds',4],
                se = summary(fit)$coefficients['bds',2],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                nonzero = length(which(dtf$bds != 0)),
                nonzero_subtype = length(which(dtf$bds != 0 & dtf$subtype == 1))
                )
        return(out)
}

### RUN NULL ANALYSES #############################################################################
# code to run to generate null results
run_null_analyses <- function(tcgadir) {
        # set genes
        genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1', 'NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32'),
                IC5 = 'ERBB2'
                )
        # calculate null stats 
        null <- list()
        for (subtype in c('IC5','IC1','IC2','IC9')) {
                for (i in 1:1000) {
                        tmp <- list()
                        for (gene in genes[[subtype]]) {
                                tmp[[gene]] <- read.delim(
                                        file.path(tcgadir, gene, 'null', paste0('avg_binders_int', i, '.txt')),
                                        as.is = TRUE
                                        )
                                colnames(tmp[[gene]]) <- c('sample', gene)

                        }
                        tmp <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample", all.x = TRUE),
                                tmp)
                        tmp <- merge(
                                tmp, 
                                tcga[,c(grep('CNA|PC|somatic|^ER$', colnames(tcga), value = TRUE), 'pam50', 'sample')],
                                by = 'sample'
                                )
                        if (subtype == 'IC5') {
                                stats <- run_subtype_associations(tmp, gene = 'ERBB2', subtype = 'Her2')
                        } else {
                                stats <- run_subtype_associations(tmp, subtype = subtype)
                        }
                        stats$iteration <- i
                        null[[paste0(subtype, i)]] <- stats
                }
        }
        null <- do.call(rbind, null)
        return(null)
}

### MAIN ##########################################################################################
# read in summary table 
tcga <- read.delim(
        'tcga_megatable.txt',
        as.is = TRUE
        )
# test subtype association 
plot_data_subtype <- rbind(
        run_subtype_associations(tcga, gene = 'ERBB2', subtype = 'Her2'),
        run_subtype_associations(tcga, subtype = 'IC1'),
        run_subtype_associations(tcga, subtype = 'IC2'),
        run_subtype_associations(tcga, subtype = 'IC9')
        )

# calculate null analyses
# file generated with run_null_analyses(tcgadir)
# read in null results
null <- read.delim(
        'primary_null_associations.txt',
        as.is = TRUE
        )

# calculate null median and 95CI 
plot_data_null <- list()
for (subtype in c('Her2','IC1','IC2','IC9')) {
        median <- median(null[null$subtype == subtype,'coef'])
        l95 <- median(null[null$subtype == subtype,'l95'])
        u95 <- median(null[null$subtype == subtype,'u95'])
        plot_data_null[[subtype]] <- data.frame(
                subtype = subtype,
                coef = median,
                l95 = l95,
                u95 = u95,
                type = 'null'
                )
}
plot_data_null <- do.call(rbind, plot_data_null)
plot_data_subtype$type <- 'real'
plot_data_subtype <- rbind(plot_data_subtype[,c('subtype','coef','l95','u95','type')], plot_data_null)
plot_data_subtype <- plot_data_subtype[order(plot_data_subtype$subtype, plot_data_subtype$type),]
plot_data_subtype$index <- 1:nrow(plot_data_subtype)

cov <- list(
        rect = list(
                col = 'transparent',
                fill = rep(c('lightpink', 'lightpink3'), 5)
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
                colours =  c('lightpink', 'lightpink3'),
                labels = c('Scrambled','Real'),
                title = 'Analysis',
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
        filename = paste0(date, '_null_compared_real_scatterplot.pdf'),
        xat = c(-2,-1, 0, 1, 2),
        xaxis.lab = c(0.15,0.37, 1, 2.70, 7.40),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = seq(1.5, nrow(plot_data_subtype), 2),
        ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        add.rectangle = TRUE,
        xleft.rectangle = -2.5,
        ybottom.rectangle = c(0, seq(4.5, nrow(plot_data_subtype), 4)),
        xright.rectangle = 2.5,
        ytop.rectangle = c(seq(2.5, nrow(plot_data_subtype), 4)),
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(0,0), x = 0.01, y = 0.05)
                ),
        x.error.right = plot_data_subtype$u95-plot_data_subtype$coef,
        x.error.left = plot_data_subtype$coef-plot_data_subtype$l95,
        width = 8,
        height = 3,
        top.padding = 2,
        resolution = 300
        )

