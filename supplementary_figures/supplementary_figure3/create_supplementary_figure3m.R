### CREATE SUPPLEMENTARY FIGURE3M #################################################################
# create supplementary figure 3M


### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survminer)
library(survival)
library(caret)
library(pec)

main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST RISK PREDICTIONS COMPARED TO CLINICAL ####################################################
bootstrap_cindex <- function(dtf, model1, model2, iterations = 100) {
	# bootstrap the predictions 
	res <- list()
	for (i in 765:iterations) {
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

### MAIN ##########################################################################################
# read in summary data
metabric <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'metabric_megatable.txt'),
	as.is = TRUE
	)

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
res_er_all <- bootstrap_cindex(meta_prog_er, 
	model1 = 'node + size + grade + age + HER2 + INTCLUST',
	model2 = 'node + size + grade + age + HER2 + INTCLUST + burden',
	iteration = 1000
	)

# run model in ER+ high risk samples 
meta_prog_her2 <- meta_prog[which(meta_prog$CLAUDIN_SUBTYPE == 'Her2' & !is.na(meta_prog$ER)),]
meta_prog_her2$burden <- sign(meta_prog_her2$ERBB2)

res_her2_all <- bootstrap_cindex(meta_prog_her2, 
	model1 = 'node + size + grade + age + ER + INTCLUST',
	model2 = 'node + size + grade + age + ER + INTCLUST + burden',
	iteration = 1000
	)

fc_er  <- median(res_er_all$model2)/median(res_er_all$model1)
fc_her2  <- median(res_her2_all$model2)/median(res_her2_all$model1)

p_er <- 1-(sum(res_er_all$model2 > res_er_all$model1)/nrow(res_er_all))
p_her2 <- 1-(sum(res_her2_all$model2 > res_her2_all$model1)/nrow(res_her2_all))


# create boxplot 
plot_data <- data.frame(
	group = c(
		rep('ER_model1', nrow(res_er_all)),
		rep('ER_model2', nrow(res_er_all)),
		rep('HER2_model1', nrow(res_her2_all)),
		rep('HER2_model2', nrow(res_her2_all))
		),
	cindex = c(
		unlist(res_er_all$model1),
		unlist(res_er_all$model2),
		unlist(res_her2_all$model1),
		unlist(res_her2_all$model2)
		)
	)

create.boxplot(
	cindex ~ group,
	data = plot_data,
	add.stripplot = TRUE,
	filename = paste0(date, '_supplementary_figure3m.png'),
	ylimits = c(0.54, 0.95),
	yat = seq(0.6,0.9,0.1),
	ylab.label = 'C-Index',
	xaxis.lab = rep(c('IC\nClinical',"IC\nClinical\nEB"),2),
	xlab.label = 'Model',
	add.rectangle = TRUE,
	xleft.rectangle = 2.5,
	ybottom.rectangle = 0,
	xright.rectangle = 7,
	ytop.rectangle = 1,
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	add.text = TRUE,
	text.labels = c('ER+','HER2+',
		paste0('FC=', round(fc_er, digits = 2)), paste0('FC=', round(fc_her2, digits = 2)), 
		paste0('P=', round(p_er, digits = 2)), paste0('P=', round(p_her2, digits = 2))
		),
	text.y = c(rep(0.93, 2), rep(0.57, 2), rep(0.55, 2)),
	text.fontface = c(rep('bold', 2), rep('plain',4)),
	text.x = rep(c(1.5,3.5),3),
	#text.x = rep(c(1.5,3.5),3),
	text.col = 'black',
	 text.cex = 1.2,
	resolution = 300
	)
