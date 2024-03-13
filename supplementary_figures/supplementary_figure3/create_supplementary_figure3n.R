### CREATE SUPPLEMENTARY FIGURE3N #################################################################
# create supplementary figure 3N
# clinical feature forest plot

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

date <- Sys.Date()
### MAIN ##########################################################################################
# read in summary data
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)

### HER2+ #########################################################################################
# consider only her2 patients
meta_her2 <- metabric[metabric$CLAUDIN_SUBTYPE == 'Her2',]
meta_her2$event <- (meta_her2$OS_STATUS == 'DECEASED')*1
meta_her2$bds <- (meta_her2$ERBB2 > 0)*1

res_her2 <- list()
for (i in c('size','age','node')) {
	fit_clinical <- glm(as.formula(paste('bds ~ ', i, ' + PC1 + PC2')), 
		data = meta_her2, family = 'binomial')
	ci <- confint(fit_clinical)
	res_her2[[i]] <- data.frame(
		subtype = 'HER2+',
		feature = i,
		coef = coef(fit_clinical)[[i]],
		p = summary(fit_clinical)$coefficients[i,4],
		l95 = ci[i,1],
		u95 = ci[i,2]
		)
}
res_her2 <- do.call(rbind, res_her2)

### ER+ HIGH RISK #################################################################################
# generate km plot for ER+ high risk tumors
meta_prog <- metabric[which(rowSums(metabric[,c('MYC_CNA','RSF1_CNA','RPS6KB1_CNA','ZNF703_CNA')]) > 0 & metabric$ER_IHC == 'pos'),]
meta_prog$burden <- sign(colSums(rbind(
	rowSums(sign(meta_prog[,c('MYC','SQLE','FBXO32')]))*meta_prog$MYC_CNA,
	rowSums(sign(meta_prog[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_prog$RPS6KB1_CNA,
	rowSums(sign(meta_prog[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_prog$RSF1_CNA,
	rowSums(sign(meta_prog[,c('ZNF703','FGFR1','LETM2')]))*meta_prog$ZNF703_CNA
	)))
meta_prog$burden <- (meta_prog$burden > 0)*1

res_er <- list()
for (i in c('size','age','node')) {
	fit_clinical <- glm(as.formula(paste('burden ~ ', i, ' + PC1 + PC2 + INTCLUST')), 
		data = meta_prog, family = 'binomial')
	ci <- confint(fit_clinical)
	res_er[[i]] <- data.frame(
		subtype = 'ER+',
		feature = i,
		coef = coef(fit_clinical)[[i]],
		p = summary(fit_clinical)$coefficients[i,4],
		l95 = ci[i,1],
		u95 = ci[i,2]
		)
}
res_er <- do.call(rbind, res_er)

### CLINICAL FEATURES ##############################################################################
plot_data <- rbind(res_er, res_her2)
plot_data$index <- 1:nrow(plot_data)

create.scatterplot(
        coef ~ index,
        data = plot_data,
        ylimits = c(-0.1,0.1),
        filename = paste0(date, '_metabric_clinical_features_scatterplot.pdf'),
        ylab.label = 'Coefficient',
        xlab.label = '',
        xaxis.lab = rep(c('Tumor Size','Age','Lymph Node'),2),
        xat = 1:nrow(plot_data),
        xaxis.rot = 90,
        xaxis.cex = 1.2,
        xlimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.h = 0,
        add.text = TRUE,
        text.labels = c('HER2+','ER+'),
        text.x = c(4.5,1.5),
        text.y = c(0.085, 0.085),
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