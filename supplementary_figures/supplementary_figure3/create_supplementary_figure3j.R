### CREATE SUPPLEMENTARY FIGURE 3J ################################################################
# create supplementary figure 3J
# compare primary vs metastatic null

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(plyr)

date <- Sys.Date()

### COMPARE PRIM VS MET ###########################################################################
compare_prim_vs_met <- function(hartwig, tcga, subtype, gene = NULL) {
        genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1','NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32'),
                IC5 = 'ERBB2'
                )

        if (!is.null(gene)) {
                tcga$bds <- sign(tcga[,gene])
                hartwig$bds <- sign(hartwig[,gene])
                if (gene == 'ERBB2') {
                        tcga_st <- tcga[tcga$pam50 == 'Her2',]
                        hartwig_st <- hartwig[which(hartwig$hormone %in% c('ER-negative/HER2-positive')),]
                } else {
                        tcga_st <- tcga[which(tcga[,paste0(gene, '_CNA')] == 1 & tcga$ER == 1),]
                        hartwig_st <- hartwig[which(hartwig[,paste0(gene, '_CNA')] == 1 & hartwig$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive')),]
                        }
        } else {
                tcga$bds <- rowSums(sign(tcga[,genes[[subtype]]]))
                hartwig$bds <- rowSums(sign(hartwig[,genes[[subtype]]]))

                tcga_st <- tcga[which(tcga[,paste0(genes[[subtype]][1], '_CNA')] == 1 & tcga$ER == 1),]
                hartwig_st <- hartwig[which(hartwig[,paste0(genes[[subtype]][1], '_CNA')] == 1 & hartwig$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive')),]
        }

        # create model data 
        model_data_tcga <- data.frame(
                disease = c(rep(0, nrow(tcga_st)), rep(1, nrow(hartwig_st))),
                bds = c(tcga_st$bds, hartwig_st$bds)
                )

        fit_tcga <- glm(disease ~ bds, data = model_data_tcga, family = 'binomial')
        ci_tcga <- confint(fit_tcga)

        out <- data.frame(
                subtype = subtype,
                gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
                coef = coef(fit_tcga)[['bds']], 
                p = summary(fit_tcga)$coefficients['bds',4],
                se = summary(fit_tcga)$coefficients['bds',2],
                l95 = ci_tcga['bds',1], 
                u95 = ci_tcga['bds',2]
                )
        return(out)
}

### RUN NULL ANALYSES #############################################################################
calculate_avg_binders <- function(allbds, thres) {
        res_avg <- aggregate(allbds[,thres], allbds[,c('sample','HLA')], mean)
        res_avg <- aggregate(res_avg$x, list(res_avg$sample), mean)
        colnames(res_avg) <- c('sample','bds')
        return(res_avg)
        }

create_null_genes_bds_table <- function(cohort, genes, dir) {
        dtf <- list()
        for (gene in genes) {
                if (cohort == 'TCGA') {
                        file <- file.path(dir, gene, 'null', paste0('avg_binders_int', i, '.txt'))
                } else {
                        file <- file.path(dir, gene, 'null', paste0('avg_binders_int', i, '.txt'))
                        }
                tmp <- read.delim(
                        file,
                        as.is = TRUE
                        )
                if (ncol(tmp) > 2) {
                        tmp <- calculate_avg_binders(tmp, thres = 'wbs_0.5_3')
                }
                colnames(tmp) <- c('sample', gene)
                dtf[[gene]] <- tmp
        }
        dtf <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample", all.x = TRUE),
                        dtf)
        return(dtf)
}

# code to run to generate null results
run_null_analyses <- function(tcgadir, hartwigdir) {
        # set genes
        genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1', 'NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32'),
                IC5 = 'ERBB2'
                )
        # calculate null stats 
        # calculate null stats 
        null <- list()
        for (subtype in c('IC5','IC1','IC2','IC9')) {
                for (i in 1:1000) {
                        tcgatmp <-  create_null_genes_bds_table('TCGA', genes[[subtype]], dir = tcgadir)
                        tcgatmp <- merge(
                                tcgatmp, 
                                tcga[,c(grep('CNA|PC|^ER$', colnames(tcga), value = TRUE), 'pam50', 'sample')],
                                by = 'sample'
                                )

                        harttmp <-  create_null_genes_bds_table('Hartwig', genes[[subtype]], dir = hartwigdir)
                        harttmp <- merge(
                                harttmp, 
                                hartwig[,c(grep('CNA|PC', colnames(hartwig), value = TRUE), 'hormone', 'sample')],
                                by = 'sample'
                                )
                        if (subtype == 'IC5') {
                                stats <- compare_prim_vs_met(tcga = tcgatmp, hartwig = harttmp, gene = 'ERBB2', subtype = 'IC5')
                        } else {
                                stats <- compare_prim_vs_met(tcga = tcgatmp, hartwig = harttmp, subtype = subtype)
                                }
                        stats$iteration <- i
                        null[[paste0(subtype, i)]] <- stats
                }
        }
        null <- do.call(rbind, null)
        return(null)
}

### PRIMARY VS METASTATIC #########################################################################
# read in summary table 
tcga <- read.delim(
        'tcga_megatable.txt',
        as.is = TRUE
        )
hartwig <- read.delim(
        'hartwig_megatable.txt',
        as.is = TRUE
        )

# test primary vs met association 
plot_data_pvm <- rbind(
        compare_prim_vs_met(tcga = tcga, hartwig = hartwig, gene = 'ERBB2', subtype = 'IC5'),
        compare_prim_vs_met(tcga = tcga, hartwig = hartwig, subtype = 'IC1'),
        compare_prim_vs_met(tcga = tcga, hartwig = hartwig, subtype = 'IC2'),
        compare_prim_vs_met(tcga = tcga, hartwig = hartwig, subtype = 'IC9')
        )

# calculate null stats 
# file generated with run_null_analyses(tcgadir, hartwigdir)
null <- read.delim(
        'primary_vs_metastatic_null_associations.txt',
        as.is = TRUE
        )

# calculate null median and 95CI 
plot_data_null <- list()
for (subtype in c('IC5','IC1','IC2','IC9')) {
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
plot_data_pvm$type <- 'real'
plot_data_pvm <- rbind(plot_data_pvm[,c('subtype','coef','l95','u95','type')], plot_data_null)
plot_data_pvm <- plot_data_pvm[order(plot_data_pvm$subtype, plot_data_pvm$type),]
plot_data_pvm$index <- 1:nrow(plot_data_pvm)

cov <- list(
        rect = list(
                col = 'transparent',
                fill = rep(c('lightpink', 'lightpink3'), 5)
                )
        );

cov.grob <- covariates.grob(
        covariates = cov,
        ord = c(1:nrow(plot_data_pvm)),
        side = 'right',
        size = 1
        );

cov.legend <- list(
        legend = list(
                colours =  rev(c('lightpink', 'lightpink3')),
                labels = rev(c('Scrambled','Real')),
                title = 'Analysis',
                border = 'transparent'
                )
        );

cov.legend.grob <- legend.grob(
        legends = cov.legend
        );

create.scatterplot(
        index ~ coef,
        data = plot_data_pvm,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        filename = paste0(date, '_null_compared_real_prim_vs_met_scatterplot.pdf'),
        xat = c(-2,-1, 0, 1, 2),
        xaxis.lab = c(0.15,0.37, 1, 2.70, 7.40),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('Her2','IC5', unique(plot_data_pvm$subtype)),
        yat = seq(1.5, nrow(plot_data_pvm), 2),
        ylimits = c(0.5, nrow(plot_data_pvm)+0.75),
        main.cex = 2,
        abline.v = 0,
        add.rectangle = TRUE,
        xleft.rectangle = -2.5,
        ybottom.rectangle = c(0, seq(4.5, nrow(plot_data_pvm), 4)),
        xright.rectangle = 2.5,
        ytop.rectangle = c(seq(2.5, nrow(plot_data_pvm), 4)),
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(0,0), x = 0.01, y = 0.05)
                ),
        x.error.right = plot_data_pvm$u95-plot_data_pvm$coef,
        x.error.left = plot_data_pvm$coef-plot_data_pvm$l95,
        width = 8,
        height = 3.75,
        top.padding = 2,
        resolution = 300
        )


