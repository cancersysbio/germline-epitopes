### CREATE SUPPLEMENTARY FIGURE 3L ################################################################
# correct survival outcomes for HLA breadth
# create supplementary figure 3L

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survminer)
library(survival)

### MAIN #######################################################################################
# read in summary data
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)

# read in hlas
hlas <- read.delim(
	'hla_promiscuity_proportions.txt',
	as.is = TRUE
	)

### HER2+ ######################################################################################
# consider only her2 patients
meta_her2 <- metabric[metabric$CLAUDIN_SUBTYPE == 'Her2',]
meta_her2$event <- (meta_her2$OS_STATUS == 'DECEASED')*1
meta_her2$bds <- (meta_her2$ERBB2 > 0)*1
# testing 5-year relapse 
meta_her2[meta_her2$OS_MONTHS > 60,'event'] <- 0
meta_her2[meta_her2$OS_MONTHS > 60,'OS_MONTHS'] <- 60

meta_her2 <- merge(meta_her2, hlas, by = 'sample', all.x = TRUE)
# calculate stats to report correcting for PCs, age
her2_fit1 <- coxph(
	Surv(OS_MONTHS,event)~bds + PC1 + PC2 + AGE_AT_DIAGNOSIS + pga + INTCLUST + stage + grade + promprop,
	data = meta_her2
	)
her2_ci1 <- confint(her2_fit1)
her2_fit2 <- coxph(
	Surv(OS_MONTHS,event)~bds + PC1 + PC2 + AGE_AT_DIAGNOSIS + pga + INTCLUST + stage + grade,
	data = meta_her2
	)
her2_ci2 <- confint(her2_fit2)

### ER+ HIGH RISK #################################################################################
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

meta_prog <- merge(meta_prog, hlas, by = 'sample', all.x = TRUE)
# calculate stats to report correcting for PCs, age
er_fit1 <- coxph(
	Surv(OS_MONTHS,event)~burden + PC1 + PC2 + AGE_AT_DIAGNOSIS + pga + INTCLUST + stage + grade + promprop,
	data = meta_prog
	)
er_ci1 <- confint(er_fit1)
er_fit2 <- coxph(
	Surv(OS_MONTHS,event)~burden + PC1 + PC2 + AGE_AT_DIAGNOSIS + pga + INTCLUST + stage + grade,
	data = meta_prog
	)
er_ci2 <- confint(er_fit2)

# create plot data 
plot_data <- data.frame(
	subtype = rep(c('HER2','ER'), each = 2), 
	comparison = rep(c('adjusted','unadjusted'), 2),
	hr = c(coef(her2_fit1)['bds'], coef(her2_fit2)['bds'],
		coef(er_fit1)['burden'], coef(er_fit2)['burden']),
	l95 = c(her2_ci1['bds',1], her2_ci2['bds',1],
		er_ci1['burden',1], er_ci2['burden',1]),
	u95 = c(her2_ci1['bds',2], her2_ci2['bds',2],
		er_ci1['burden',2], er_ci2['burden',2])
	)
plot_data$index <- rev(1:nrow(plot_data))

create.scatterplot(
        index ~ hr,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-2,2),
        #xat = log(c(0.15, 0.5, 1, 2.5, 8)),
        #xaxis.lab = c('0.15','0.50','1.00','2.50','8.00'),
        filename = paste0(date, '_metabric_hla_adjusted_survival_scatterplot.pdf'),
        xlab.label = 'Hazard Ratio',
        ylab.label = 'GEB',
        ylimits = c(0.5, nrow(plot_data)+0.5),
        yaxis.lab = rep(c('Unadjusted','Adjusted'),2),
        yat = 1:nrow(plot_data),
        add.rectangle = TRUE,
        xleft.rectangle = -5,
        ybottom.rectangle = 0,
        xright.rectangle = 5,
        ytop.rectangle = 2.5,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        abline.v = 0,
        width = 10,
        height = 3,
        add.text = TRUE,
        text.labels = c('PAM50 HER2-enriched','ER+ (IC1/IC2/IC6/IC9)'),
        text.y = c(4,2),
        text.x = rep(-1.5,2),
        top.padding = 2,
        x.error.right = plot_data$u95-plot_data$hr,
        x.error.left = plot_data$hr-plot_data$l95,
        resolution = 300
        )