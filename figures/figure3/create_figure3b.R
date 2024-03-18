### CREATE FIGURE3B ###############################################################################
# create figure 3B

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}
date <- Sys.Date()

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
  file.path(main_repo_path,'data','cohort_megatables','hartwig_megatable.txt'), 
	as.is = TRUE
	)
# read in tcga summary data 
tcga <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','tcga_megatable.txt'), 
	as.is = TRUE
	)
icgc <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','icgc_megatable.txt'), 
	as.is = TRUE
	)

# test primary vs met association 
plot_data_pvm <- rbind(
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, gene = 'ERBB2', subtype = 'IC5'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC1'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC2'),
	compare_prim_vs_met(tcga = tcga, icgc = icgc, hartwig = hartwig, subtype = 'IC9')
	)

meta_plot_data <- list()
for (subtype in unique(plot_data_pvm$subtype)) {
	meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                        data = plot_data_pvm[plot_data_pvm$subtype == subtype,])
        res <- rma(yi, vi, data=meta_data)
        meta_plot_data[[subtype]] <- data.frame(
                        subtype = subtype,
                        coef = res$beta[[1]],
                        l95 = res$ci.lb,
                        u95 = res$ci.ub,
                        p = res$pval
                        )
}
meta_plot_data <- do.call(rbind, meta_plot_data)
meta_plot_data$index <- 1:nrow(meta_plot_data)

create.scatterplot(
        index ~ coef,
        data = meta_plot_data,
        horizontal = TRUE,
        xlimits = c(-2,2),
        filename = paste0(date, '_figure3b.png'),
        xat = c(log(c(0.2,0.5)),0, log(c(2,6))),
        xaxis.lab = c('0.2','0.5', '1.0', '2.0', '6.0'),
        xlab.label = '',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('IC5','HER2+', unique(meta_plot_data$subtype)),
        yat = 1:nrow(meta_plot_data),
        ylimits = c(0.5, nrow(meta_plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        x.error.right = meta_plot_data$u95-meta_plot_data$coef,
        x.error.left = meta_plot_data$coef-meta_plot_data$l95,
        width = 8,
        height = 3,
        top.padding = 2,
        resolution = 300
        )