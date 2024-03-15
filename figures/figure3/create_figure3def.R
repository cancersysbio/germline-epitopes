### CREATE FIGURE3 ################################################################################
# create figure 3

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survminer)
library(survival)
library(caret)
library(pec)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}
date <- Sys.Date()

### TEST RISK PREDICTIONS COMPARED TO CLINICAL ####################################################
bootstrap_cindex <- function(dtf, model1, model2, iterations = 100) {
	# bootstrap the predictions 
	res <- list()
	for (i in 1:iterations) {
		bdata <- dtf[sample(1:nrow(dtf), nrow(dtf), replace = TRUE),]
		# fit model with only IHC
		fit1 <- coxph(as.formula(paste('Surv(OS_MONTHS, event) ~', model1)), data = bdata, x = TRUE)
		predict1 <- pec::cindex(fit1, formula = as.formula(paste('Surv(OS_MONTHS, event) ~', model1)), 
			data = bdata)
		# fit model with IHC and IC
		fit2 <- coxph(as.formula(paste('Surv(OS_MONTHS, event) ~', model2)), data = bdata, x = TRUE)
		predict2 <- pec::cindex(fit2, formula = as.formula(paste('Surv(OS_MONTHS, event) ~', model2)), 
			data = bdata)
		# add results to dataframe
		res[[i]] <- data.frame(
			index = i,
			model1 = predict1$AppCindex[[1]],
			model2 = predict2$AppCindex[[1]]
			)
	}
	res <- do.call(rbind, res)
	return(res)
}

calculate_ci <- function(x) {
	x <- x[order(x)]
	lrank <- round(length(x)*0.025)
	urank <- round(length(x)*0.975)
	return(x[c(lrank, urank)])
}

### MAIN #######################################################################################
# read in summary data
metabric <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','metabric_megatable.txt'),
	as.is = TRUE
	)

### FIGURE 3B - HER2+ #############################################################################
# consider only her2 patients
meta_her2 <- metabric[metabric$CLAUDIN_SUBTYPE == 'Her2',]
meta_her2$event <- (meta_her2$OS_STATUS == 'DECEASED')*1
meta_her2$bds <- (meta_her2$ERBB2 > 0)*1
# testing 5-year relapse 
meta_her2[meta_her2$OS_MONTHS > 60,'event'] <- 0
meta_her2[meta_her2$OS_MONTHS > 60,'OS_MONTHS'] <- 60
fit1 <- survfit(Surv(OS_MONTHS,event)~bds,meta_her2)
family_type <- 'sans'
pdf(paste0(date, '_metabric_her2_survival_KM.pdf'))
ggsurvplot(fit1, data = meta_her2, pval = F, conf.int = F, risk.table = TRUE,
	   palette = c("#8fbdd9ff", "#fa8c57ff"), legend.labs = c("Low","High"), legend = c(0.8, 0.25), legend.title="",
           xlab = 'Month',
           ggtheme = theme_classic() + theme(axis.text.x = element_text(family=family_type,size = 14),
                                             axis.title =  element_text(family=family_type,size = 16,face = "bold"),
                                             axis.text.y = element_text(family=family_type,size = 14),
                                             legend.background = element_rect(fill = F),
                                             title =  element_text(family=family_type),
                                             text = element_text(family=family_type),
                                             legend.text =element_text(family=family_type, size = 14)
                                             ))
dev.off()

# calculate stats to report correcting for PCs, age
fit <- coxph(
	Surv(OS_MONTHS,event)~bds + PC1 + PC2 + AGE_AT_DIAGNOSIS + stage + grade + pga + INTCLUST,
	data = meta_her2
	)


### FIGURE 3C - ER+ HIGH RISK #####################################################################
# generate km plot for ER+ high risk tumors
meta_prog <- metabric[which(rowSums(metabric[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & metabric$ER_IHC == 'pos'),]
meta_prog$burden <- sign(colSums(rbind(
	rowSums(sign(meta_prog[,c('MYC','SQLE','FBXO32')]))*meta_prog$MYC_CNA,
	rowSums(sign(meta_prog[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_prog$RPS6KB1_CNA,
	rowSums(sign(meta_prog[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_prog$RSF1_CNA,
	rowSums(sign(meta_prog[,c('ZNF703','FGFR1','LETM2')]))*meta_prog$ZNF703_CNA
	)))


# add pga because calculate burden across multiple amplicons
meta_prog$event <- (meta_prog$OS_STATUS == 'DECEASED')*1
meta_prog$burden <- (meta_prog$burden > 0)*1
# considering 5-year relapse
meta_prog[which(meta_prog$OS_MONTHS > 60),'event'] <- 0
meta_prog[which(meta_prog$OS_MONTHS > 60),'OS_MONTHS'] <- 60
fitall <- survfit(Surv(OS_MONTHS,event)~burden,meta_prog)
family_type <- 'sans'
pdf(paste0(date, '_metabric_high_risk_survival_KM.pdf'))
ggsurvplot(fitall, data = meta_prog, pval = F, conf.int = F, risk.table = TRUE,
	   palette = c("#8fbdd9ff", "#fa8c57ff"), legend.labs = c("Low","High"), legend = c(0.8, 0.25), legend.title="",
           xlab = 'Month',
           ggtheme = theme_classic() + theme(axis.text.x = element_text(family=family_type,size = 14),
                                             axis.title =  element_text(family=family_type,size = 16,face = "bold"),
                                             axis.text.y = element_text(family=family_type,size = 14),
                                             legend.background = element_rect(fill = F),
                                             title =  element_text(family=family_type),
                                             text = element_text(family=family_type),
                                             legend.text =element_text(family=family_type, size = 14)
                                             ))
dev.off()

# calculate stats to report correcting for PCs, age, intclust and pga
fitcox <- coxph(
	Surv(OS_MONTHS,event)~burden + PC1 + PC2 + AGE_AT_DIAGNOSIS + INTCLUST + stage + grade + pga, 
	data = meta_prog
	)

### FIGURE 3D #####################################################################################
# reformat metabric data
meta_prog <- metabric
meta_prog$event <- (meta_prog$OS_STATUS == 'DECEASED')*1
meta_prog[which(meta_prog$OS_MONTHS > 60),'event'] <- 0
meta_prog[which(meta_prog$OS_MONTHS > 60),'OS_MONTHS'] <- 60
meta_prog$ER <- NA
meta_prog[which(meta_prog$ER_IHC == 'pos'),'ER'] <- 1
meta_prog[which(meta_prog$ER_IHC == 'neg'),'ER'] <- 0
meta_prog$HER2 <- NA
meta_prog[which(meta_prog$HER2_SNP6 == 'GAIN'),'HER2'] <- 1
meta_prog[which(meta_prog$HER2_SNP6 == 'NEUT'),'HER2'] <- 0
meta_prog[which(meta_prog$HER2_SNP6 == 'LOSS'),'HER2'] <- -1
# truncate number lymph nodes at 10, same as Rueda et al. 
meta_prog[which(meta_prog$node > 10),'node'] <- 10
# remove any NAs
meta_prog <- meta_prog[-unique(which(is.na(meta_prog[,c('node','size','grade','stage','HER2')]), arr.ind = TRUE)[,'row']),]

# run model in ER+ high risk samples 
meta_prog_er <- meta_prog[which(rowSums(meta_prog[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & meta_prog$ER_IHC == 'pos'),]
meta_prog_er$burden <- sign(colSums(rbind(
	rowSums(sign(meta_prog_er[,c('MYC','SQLE','FBXO32')]))*meta_prog_er$MYC_CNA,
	rowSums(sign(meta_prog_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_prog_er$RPS6KB1_CNA,
	rowSums(sign(meta_prog_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_prog_er$RSF1_CNA,
	rowSums(sign(meta_prog_er[,c('ZNF703','FGFR1','LETM2')]))*meta_prog_er$ZNF703_CNA
	)))
res_er <- bootstrap_cindex(meta_prog_er, 
	model1 = 'INTCLUST',
	model2 = 'INTCLUST + burden',
	iteration = 1000
	)

# run model in ER+ high risk samples 
meta_prog_her2 <- meta_prog[which(meta_prog$CLAUDIN_SUBTYPE == 'Her2' & !is.na(meta_prog$ER)),]
meta_prog_her2$burden <- sign(meta_prog_her2$ERBB2)

res_her2 <- bootstrap_cindex(meta_prog_her2, 
	model1 = 'INTCLUST',
	model2 = 'INTCLUST + burden',
	iteration = 1000
	)

plot_data <- data.frame(
	group = c('ER_model1','ER_model2','HER2_model1','HER2_model2'),
	cindex = c(median(res_er$model1), median(res_er$model2), median(res_her2$model1), median(res_her2$model2)),
	l95 = c(calculate_ci(res_er$model1)[1], calculate_ci(res_er$model2)[1], 
		calculate_ci(res_her2$model1)[1], calculate_ci(res_her2$model2)[1]),
	u95 = c(calculate_ci(res_er$model1)[2], calculate_ci(res_er$model2)[2], 
		calculate_ci(res_her2$model1)[2], calculate_ci(res_her2$model2)[2]),
	index = 1:4
	)

fc_er  <- median(res_er$model2)/median(res_er$model1)
fc_her2  <- median(res_her2$model2)/median(res_her2$model1)

p_er <- 1-(sum(res_er$model2 > res_er$model1)/nrow(res_er))
p_her2 <- 1-(sum(res_her2$model2 > res_her2$model1)/nrow(res_her2))

# create plot
create.scatterplot(
        index ~ cindex,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(0.45,0.85),
        filename = paste0(date, '_cindex_model_scatterplot.pdf'),
        xlab.label = 'C-Index',
        ylab.label = 'Model',
        yaxis.lab = rep(c('IntClust','IntClust\nEpitopes'), 2),
        yat = 1:nrow(plot_data),
        ylimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        add.rectangle = TRUE,
        xleft.rectangle = 0.45,
        ybottom.rectangle = 0,
        xright.rectangle = 0.85,
        ytop.rectangle = 2.5,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        x.error.right = plot_data$u95-plot_data$cindex,
        x.error.left = plot_data$cindex-plot_data$l95,
        width = 8,
        height = 3,
        top.padding = 2,
        add.text = TRUE,
        text.labels = c('HER2+','ER+',
			paste0('FC=', round(fc_her2, digits = 2)), paste0('FC=', round(fc_er, digits = 2)), 
			paste0('P=', round(p_her2, digits = 2)), paste0('P=', round(p_er, digits = 2))
			),
        text.y = c(4.3,2.2,3.8,1.7,3.4,1.3),
        text.fontface = c(rep('bold', 2), rep('plain',4)),
        text.x = rep(0.48,6),
        text.col = 'black',
        text.cex = 1.2,
        resolution = 300
        )