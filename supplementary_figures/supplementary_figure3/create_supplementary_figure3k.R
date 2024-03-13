### CREATE SUPPLEMENTARY FIGURE 3K ################################################################
# create supplementary figure 3k
# test prognostic association in ER+ subtypes individually

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survminer)
library(survival)
library(caret)
library(pec)

### MAIN #######################################################################################
# read in summary data
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)

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

### Subtypes ######################################################################################
subtypes <- list()
subtypes[['IC1']] <- meta_prog[which(meta_prog$RPS6KB1_CNA > 0 & meta_prog$ER_IHC == 'pos'),]
subtypes[['IC2']] <- meta_prog[which(meta_prog$RSF1_CNA > 0 & meta_prog$ER_IHC == 'pos'),]
subtypes[['IC6']] <- meta_prog[which(meta_prog$ZNF703_CNA > 0 & meta_prog$ER_IHC == 'pos'),]
subtypes[['IC9']] <- meta_prog[which(meta_prog$MYC_CNA > 0 & meta_prog$ER_IHC == 'pos'),]

# calculate stats to report correcting for PCs, age, intclust and pga
res <- list()
for (i in c('IC1','IC2','IC6','IC9')) {
	fitcox <- coxph(
		Surv(OS_MONTHS,event)~burden + PC1 + PC2 + AGE_AT_DIAGNOSIS + INTCLUST, 
		data = subtypes[[i]]
		)
	ci <- confint(fitcox)
	res[[i]] <- data.frame(
		subtype = i,
		hr = coef(fitcox)['burden'],
		p = summary(fitcox)$coefficients['burden',5],
		l95 = ci['burden',1], 
		u95 = ci['burden',2]
		)
}
res <- do.call(rbind, res)


create.scatterplot(
        index ~ hr,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-0.5,1.7),
        xat = log(c(0.15, 0.5, 1, 2.5, 8)),
        xaxis.lab = c('0.15','0.50','1.00','2.50','8.00'),
        filename = paste0(date, '_IC_survival_scatterplot.pdf'),
        xlab.label = 'Hazard Ratio',
        ylab.label = 'Subtype',
        ylimits = c(0.5, nrow(plot_data)+0.5),
        yaxis.lab = paste0('IC', c(1,2,6,9)),
        yat = 1:nrow(plot_data),
        abline.v = 0,
        width = 10,
        height = 3,
        top.padding = 2,
        x.error.right = plot_data$u95-plot_data$hr,
        x.error.left = plot_data$hr-plot_data$l95,
        resolution = 300
        )
