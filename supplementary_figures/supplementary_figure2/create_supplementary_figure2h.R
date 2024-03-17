### SUPPLEMENTARY FIGURE 2H #######################################################################
# test all definitions of her2 in TCGA, ICGC, GEL and METABRIC
# create supplementary figure 2H

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST ASSOCIATIONS #############################################################################
test_tcga_icgc_association <- function(tcga, subtype) {
        fit <- glm(
                as.formula(paste(subtype, '~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic')),
                data = tcga,
                family = 'binomial'
                )
        ci <- confint(fit)
        res <- data.frame(
                subtype = subtype,
                coef = coef(fit)[['bds']],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                num = sum(tcga[,subtype] == 1)
                )
        return(res)
}
test_metabric_association <- function(tcga, subtype) {
        fit <- glm(
                as.formula(paste(subtype, '~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6')),
                data = tcga,
                family = 'binomial'
                )
        ci <- confint(fit)
        res <- data.frame(
                subtype = subtype,
                coef = coef(fit)[['bds']],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                num = sum(tcga[,subtype] == 1)
                )
        return(res)
}

### TCGA ##########################################################################################
# read in summary data
tcga <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
        as.is = TRUE
        )
tcga$bds <- sign(tcga$ERBB2)

tcga$pam50_subtype <- (tcga$pam50 == 'Her2')*1
tcga$cna_subtype <- (tcga$ERBB2_CNA_5copies == 1)*1
tcga$ihc_subtype <- (tcga$HER2.newly.derived == 'Positive')*1
tcga$ihc_erneg_subtype <- (tcga$HER2.newly.derived == 'Positive' & tcga$ER == 0)*1
tcga$ihc_erpos_subtype <- (tcga$HER2.newly.derived == 'Positive' & tcga$ER == 1)*1
tcga$cna_erneg_subtype <- (tcga$ERBB2_CNA_5copies == 1 & tcga$ER == 0)*1


tplot_data <- rbind(
        test_tcga_icgc_association(tcga, 'pam50_subtype'),
        test_tcga_icgc_association(tcga, 'cna_subtype'),
        test_tcga_icgc_association(tcga, 'ihc_subtype'),
        test_tcga_icgc_association(tcga, 'ihc_erneg_subtype'),
        test_tcga_icgc_association(tcga, 'ihc_erpos_subtype'),
        test_tcga_icgc_association(tcga, 'cna_erneg_subtype')
        )
tplot_data$cohort <- 'TCGA'

### METABRIC #####################################################################################
# read in summary data
metabric <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables','metabric_megatable.txt'),
        as.is = TRUE
        )
metabric$bds <- sign(metabric$ERBB2)
metabric$pam50_subtype <- (metabric$CLAUDIN_SUBTYPE== 'Her2')*1

# test subtype association 
mplot_data <- test_metabric_association(metabric, 'pam50_subtype')
mplot_data$cohort <- 'METABRIC'

### ICGC ##########################################################################################
# read in summary data
icgc <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables','icgc_megatable.txt'),
	as.is = TRUE
	)
icgc$bds <- sign(icgc$ERBB2)
icgc$ihc_subtype <- (icgc$final.HER2 == 'positive')*1
icgc$ihc_erneg_subtype <- (icgc$final.HER2 == 'positive' & icgc$final.ER == 'negative')*1
icgc$ihc_erpos_subtype <- (icgc$final.HER2 == 'positive' & icgc$final.ER == 'positive')*1
icgc$cna_subtype <- (icgc$ERBB2_CNA_5copies == 1)*1
icgc$cna_erneg_subtype <- (icgc$ERBB2_CNA_5copies == 1 & icgc$final.ER == 'negative')*1
colnames(icgc) <- gsub('snvs','somatic', colnames(icgc))

iplot_data <- rbind(
        test_tcga_icgc_association(icgc, 'ihc_subtype'),
        test_tcga_icgc_association(icgc, 'ihc_erneg_subtype'),
        test_tcga_icgc_association(icgc, 'ihc_erpos_subtype'),
        test_tcga_icgc_association(icgc, 'cna_subtype'),
        test_tcga_icgc_association(icgc, 'cna_erneg_subtype')
        )
iplot_data$cohort <- 'ICGC'

### GEL ###########################################################################################
# read in summary data
gplot_data <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables','gel_her2_definition_associations.txt'), 
        as.is = TRUE
        )

### PLOT ###########################################################################################
# create plot data
plot_data <- rbind(tplot_data, iplot_data, gplot_data, mplot_data)
plot_data <- plot_data[order(plot_data$subtype, plot_data$cohort),]
# reorder slightly
plot_data$index <- c(13:15,4:6,10:12,1:3, 7:9,16,17)
plot_data <- plot_data[order(plot_data$index),]

create.scatterplot(
        index ~ coef,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        filename = paste0(date, '_GEL_ICGC_TCGA_METABRIC_her2_definitions_scatterplot.png'),
        xat = c(-2,-1, 0, 1, 2),
        xaxis.lab = c(0.15,0.37, 1, 2.70, 7.40),
        xlab.label = 'Odds Ratio',
        ylab.label = 'HER2 Definition',
        yaxis.lab = c('Clinical\nHER2+/ER+','CNA (>5)\nHER2+','Clinical\nHER2+', 'Clinical\nHER2+/ER-', 'CNA (>5)\nHER2+/ER-','PAM50 HER2+'),
        yat = c(2,5,8,11,14,16.5),
        ylimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        add.text = TRUE,
        text.labels = paste0(plot_data$cohort, '(n=', plot_data$num, ')'),
        text.y = 1:nrow(plot_data),
        text.x = rep(1.65,nrow(plot_data)),
        text.cex = 0.8,
        x.error.right = plot_data$u95-plot_data$coef,
        x.error.left = plot_data$coef-plot_data$l95,
        abline.h = c(3.5,6.5,9.5,12.5,15.5),
        abline.lty = 2,
        width = 8.5,
        height = 5,
        top.padding = 2,
        left.padding = 5,
        resolution = 300
        )

