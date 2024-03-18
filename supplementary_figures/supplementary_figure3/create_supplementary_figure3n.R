### CREATE SUPPLEMENTARY FIGURE3N #################################################################
# create supplementary figure 3N
# clinical feature forest plot

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST CLINICAL ASSOCIATIONS ####################################################################
test_clinical_associations <- function(dtf, features, covariate = NULL) {
	# set covarites 
	covariates <- 'bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6'
	if (!is.null(covariate)) {
		covariates <- paste0(covariates, '+', covariate)
	}
	# test associations
	res <- list()
	for (i in features) {
		fit <- glm(as.formula(paste(i, ' ~ ', covariates)), 
			data = dtf, family = 'binomial')
		ci <- confint(fit)
		res[[i]] <- data.frame(
			subtype = 'HER2+',
			feature = i,
			coef = coef(fit)[['bds']],
			se = summary(fit)$coefficients['bds',2],
			p = summary(fit)$coefficients['bds',4],
			l95 = ci['bds',1],
			u95 = ci['bds',2]
			)
	}
	res <- do.call(rbind, res)
	return(res)
}

### RUN META ANALYSIS #############################################################################
run_meta_analysis <- function(dtf, subtype, feature) {
	meta_data <- escalc(measure="OR", 
		yi = coef, 
		sei = se, 
		data = dtf
		)
	res <- rma(yi, vi, data=meta_data)
	res_dtf <- data.frame(
		subtype = subtype,
		feature = feature,
		coef = res$beta[[1]],
		l95 = res$ci.lb,
		u95 = res$ci.ub,
		p = res$pval
		)
	return(res_dtf)
	}

run_meta_analysis_wrapper <- function(dtf, subtype, features) {
	plot_data <- list()
	for (i in features) {
		plot_data[[i]] <- run_meta_analysis(
			dtf = dtf[dtf$feature == i,],
			feature = i,
			subtype = subtype
			)
		}
	plot_data <- do.call(rbind, plot_data)
	return(plot_data)
	}

### MAIN ##########################################################################################
# read in summary data
metabric <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'metabric_megatable.txt'),
	as.is = TRUE
	)

tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','tcga_megatable.txt'),
	as.is = TRUE
	)

icgc <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','icgc_megatable.txt'),
	as.is = TRUE
	)

### HER2+ METABRIC ################################################################################
# consider only her2 patients
meta_her2 <- metabric[metabric$CLAUDIN_SUBTYPE == 'Her2',]
meta_her2$event <- (meta_her2$OS_STATUS == 'DECEASED')*1
meta_her2$bds <- (meta_her2$ERBB2 > 0)*1
meta_her2$node_bin <- (meta_her2$node > 1)*1
meta_her2$stage <- (meta_her2$stage > 2)*1
meta_her2$age <- (meta_her2$age > 45)*1

res_her2 <- test_clinical_associations(dtf = meta_her2, 
	features = c('stage','age','node_bin'))
### HER2+ TCGA ####################################################################################
tcga_her2 <- tcga[which(tcga$pam50 == 'Her2'),]
tcga_her2$bds <- (tcga_her2$ERBB2 > 0)*1
tcga_her2$stage <- (tcga_her2$stage %in% c('III','IV'))*1
tcga_her2[which(tcga_her2$stage == '[Not vailable]'),'stage_num'] <- NA
tcga_her2$age <- (tcga_her2$age > 45)*1
tcga_her2$ic10 <- as.factor(tcga_her2$ic10)

tres_her2 <- test_clinical_associations(dtf = tcga_her2, 
	features = c('stage','age'),
	covariate = 'ic10')

### ER+ HIGH RISK METABRIC ########################################################################
# generate km plot for ER+ high risk tumors
meta_prog <- metabric[which(rowSums(metabric[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & metabric$ER_IHC == 'pos'),]
meta_prog$bds <- sign(colSums(rbind(
	rowSums(sign(meta_prog[,c('MYC','SQLE','FBXO32')]))*meta_prog$MYC_CNA,
	rowSums(sign(meta_prog[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_prog$RPS6KB1_CNA,
	rowSums(sign(meta_prog[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_prog$RSF1_CNA,
	rowSums(sign(meta_prog[,c('ZNF703','FGFR1','LETM2')]))*meta_prog$ZNF703_CNA
	)))
meta_prog$bds <- (meta_prog$bds > 0)*1
meta_prog$stage <- (meta_prog$stage > 2)*1
meta_prog$node_bin <- (meta_prog$node > 1)*1
meta_prog$age <- (meta_prog$age > 45)*1

res_er <- test_clinical_associations(dtf = meta_prog, 
	features = c('stage','age','node_bin'),
	covariate = 'INTCLUST')

### ER+ HIGH TCGA #################################################################################
tcga_er <- tcga[which(rowSums(tcga[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & tcga$ER == 1),]
tcga_er$bds <- sign(colSums(rbind(
	rowSums(sign(tcga_er[,c('MYC','SQLE','FBXO32')]))*tcga_er$MYC_CNA,
	rowSums(sign(tcga_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*tcga_er$RPS6KB1_CNA,
	rowSums(sign(tcga_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*tcga_er$RSF1_CNA,
	rowSums(sign(tcga_er[,c('ZNF703','FGFR1','LETM2')]))*tcga_er$ZNF703_CNA
	)))
tcga_er$bds <- (tcga_er$bds > 0)*1
tcga_er$stage <- (tcga_er$stage %in% c('III','IV'))*1
tcga_er[which(tcga_er$stage == '[Not vailable]'),'stage'] <- NA
tcga_er$age <- (tcga_er$age > 45)*1
tcga_er$ic10 <- as.factor(tcga_er$ic10)

tres_er <- test_clinical_associations(dtf = tcga_er, 
	features = c('stage','age'),
	covariate = 'ic10')
### ER+ HIGH ICGC #################################################################################
icgc_er <- icgc[which(rowSums(icgc[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & icgc$final.ER == 'positive'),]
icgc_er$bds <- sign(colSums(rbind(
	rowSums(sign(icgc_er[,c('MYC','SQLE','FBXO32')]))*icgc_er$MYC_CNA,
	rowSums(sign(icgc_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*icgc_er$RPS6KB1_CNA,
	rowSums(sign(icgc_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*icgc_er$RSF1_CNA,
	rowSums(sign(icgc_er[,c('ZNF703','FGFR1','LETM2')]))*icgc_er$ZNF703_CNA
	)))
icgc_er$bds <- (icgc_er$bds > 0)*1
icgc_er$stage <- (as.numeric(gsub('T','',icgc_er$T_stage)) > 2)*1
icgc_er$node_bin <- (icgc_er$N_stage %in% c('N2','N3'))*1
icgc_er[which(icgc_er$N_stage %in% c('Nx','NX')),'node_bin'] <- NA
icgc_er$age <- (icgc_er$age > 45)*1

ires_er <- test_clinical_associations(dtf = icgc_er, 
	features = c('stage','age','node_bin'))
### RUN META_ANALYSIS #############################################################################
# run meta analysis
plot_data_her2 <- run_meta_analysis_wrapper(
	dtf = rbind(res_her2, tres_her2),
	subtype = 'HER2',
	features = c('stage','age','node_bin')
	)
plot_data_er <- run_meta_analysis_wrapper(
	dtf = rbind(res_er, tres_er, ires_er),
	subtype = 'ER',
	features = c('stage','age','node_bin')
	)

plot_data <- rbind(
	plot_data_er,
	plot_data_her2
	)
plot_data$index <- 1:nrow(plot_data)

create.scatterplot(
        coef ~ index,
        data = plot_data,
        ylimits = c(-2.5,2.5),
        #ylimits = c(-0.1,0.1),
        filename = paste0(date, '_supplementary_figure3n.png'),
        ylab.label = 'Coefficient',
        xlab.label = '',
        xaxis.lab = rep(c('Stage','Age','Lymph Node'),2),
        xat = 1:nrow(plot_data),
        xaxis.rot = 90,
        xaxis.cex = 1.2,
        xlimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.h = 0,
        add.text = TRUE,
        text.labels = c('HER2+','ER+'),
        text.x = c(5,2),
        text.y = c(2.3,2.3),
        text.cex = 1.2,
        y.error.up = plot_data$u95-plot_data$coef,
        y.error.down = plot_data$coef-plot_data$l95,
        add.rectangle = TRUE,
        xleft.rectangle = 0,
        ybottom.rectangle = -2.5,
        xright.rectangle = 3.5,
        ytop.rectangle = 2.5,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        width = 4.5,
        height = 5.5,
        top.padding = 2,
        resolution = 300
        )