### CREATE SUPPLEMENTARY FIGURE 3I ################################################################
# per cohort primary vs metastatic comparisons
# creates supplementary figure 3i

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

### COMPARE PRIMARY VS METASTASIS #################################################################
compare_prim_vs_met <- function(hartwig, tcga, icgc, subtype, gene = NULL) {
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
		icgc$bds <- sign(icgc[,gene])
		if (gene == 'ERBB2') {
			tcga_st <- tcga[tcga$pam50 == 'Her2',]
			icgc_st <- icgc[which(icgc$final.HER2 == 'positive' & icgc$final.ER == 'negative'),]
			hartwig_st <- hartwig[which(hartwig$hormone %in% c('ER-negative/HER2-positive')),]
		} else {
			tcga_st <- tcga[which(tcga[,paste0(gene, '_CNA')] == 1 & tcga$ER == 1),]
			icgc_st <- icgc[which(icgc[,paste0(gene, '_CNA')] == 1 & icgc$final.ER == 'positive'),]
			hartwig_st <- hartwig[which(hartwig[,paste0(gene, '_CNA')] == 1 & hartwig$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive')),]
			}
	} else {
		tcga$bds <- rowSums(sign(tcga[,genes[[subtype]]]))
		icgc$bds <- rowSums(sign(icgc[,genes[[subtype]]]))
		hartwig$bds <- rowSums(sign(hartwig[,genes[[subtype]]]))

		tcga_st <- tcga[which(tcga[,paste0(genes[[subtype]][1], '_CNA')] == 1 & tcga$ER == 1),]
		icgc_st <- icgc[which(icgc[,paste0(genes[[subtype]][1], '_CNA')] == 1 & icgc$final.ER == 'positive'),]
		hartwig_st <- hartwig[which(hartwig[,paste0(genes[[subtype]][1], '_CNA')] == 1 & hartwig$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive')),]
	}

	# create model data 
	model_data_tcga <- data.frame(
		disease = c(rep(0, nrow(tcga_st)), rep(1, nrow(hartwig_st))),
		bds = c(tcga_st$bds, hartwig_st$bds)
		)

	fit_tcga <- glm(disease ~ bds, data = model_data_tcga, family = 'binomial')
	ci_tcga <- confint(fit_tcga)

	model_data_icgc <- data.frame(
		disease = c(rep(0, nrow(icgc_st)), rep(1, nrow(hartwig_st))),
		bds = c(icgc_st$bds, hartwig_st$bds)
		)

	fit_icgc <- glm(disease ~ bds, data = model_data_icgc, family = 'binomial')
	ci_icgc <- confint(fit_icgc)

	out <- data.frame(
		subtype = subtype,
		cohort = c('TCGA', 'ICGC'),
		gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
		coef = c(coef(fit_tcga)[['bds']], coef(fit_icgc)[['bds']]),
		p = c(
			summary(fit_tcga)$coefficients['bds',4],
			summary(fit_icgc)$coefficients['bds',4]
			),
		se = c(
			summary(fit_tcga)$coefficients['bds',2],
			summary(fit_icgc)$coefficients['bds',2]
			),
		l95 = c(ci_tcga['bds',1], ci_icgc['bds',1]),
		u95 = c(ci_tcga['bds',2], ci_icgc['bds',2])
		)
	return(out)
}
### SUBTYPE #######################################################################################
# read in summary data
hartwig <- read.delim(
	'hartwig_megatable.txt', 
	as.is = TRUE
	)
# read in tcga summary data 
tcga <- read.delim(
	'tcga_megatable.txt', 
	as.is = TRUE
	)
icgc <- read.delim(
	'icgc_megatable.txt', 
	as.is = TRUE
	)

# test primary vs met association 
plot_data_pvm <- rbind(
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, gene = 'ERBB2', subtype = 'IC5'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC1'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC2'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC9')
	)

plot_data_pvm <- plot_data_pvm[order(plot_data_pvm$subtype, plot_data_pvm$cohort),]
plot_data_pvm$index <- 1:nrow(plot_data_pvm)


cohort_col <- c("darkorchid4", 'deeppink3')
names(cohort_col) <- c('TCGA','ICGC')
cohort_cov <- as.character(plot_data_pvm$cohort)
for (i in names(cohort_col)) {
        cohort_cov[cohort_cov == i] <- cohort_col[i]
}


cov <- list(
        rect = list(
                col = 'transparent',
                fill = cohort_cov
                )
        );

cov.grob <- covariates.grob(
        covariates = cov,
        ord = c(1:length(cohort_cov)),
        side = 'right',
        size = 1
        );

cov.legend <- list(
        legend = list(
                colours = cohort_col[c('TCGA','ICGC')],
                labels = c('TCGA','ICGC'),
                title = 'Cohort',
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
        xlimits = c(-2.2,2.2),
        filename = paste0(date, '_primary_vs_met_scatterplot.pdf'),
        xat = c(log(c(0.2,0.5)),0, log(c(2,6))),
        xaxis.lab = c('0.2','0.5', '1.0', '2.0', '6.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('IC5','HER2+', unique(plot_data_pvm$subtype)),
        yat = seq(1.5, nrow(plot_data_pvm), 2),
        ylimits = c(0.5, nrow(plot_data_pvm)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(0,0), x = 0.01, y = 0.025)
                ),
        add.rectangle = TRUE,
        xleft.rectangle = -3,
        ybottom.rectangle = c(0, seq(4.5, nrow(plot_data_pvm), 4)),
        xright.rectangle = 3,
        ytop.rectangle = seq(2.5, nrow(plot_data_pvm), 4),
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        x.error.right = plot_data_pvm$u95-plot_data_pvm$coef,
        x.error.left = plot_data_pvm$coef-plot_data_pvm$l95,
        width = 8,
        height = 3,
        top.padding = 2,
        resolution = 300
        )
