### CREATE SUPPLEMENTARY FIGURE 2G ################################################################
# create supplementary figure 2G

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

date <- Sys.Date()
### TEST SUBTYPE ASSOCIATION ######################################################################
run_icgc_subtype_associations <- function(dtf, subtype, gene = NULL) {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32'),
		IC10 = c('FOXC1','MIA','MELK'),
		IC5 = 'ERBB2'
		)

	if (!is.null(gene)) {
		dtf$bds <- sign(dtf[,gene])
		if (subtype == 'Her2') {
			dtf$subtype <- (dtf$final.HER2 == 'positive' & dtf$final.ER == 'negative')*1
		} else {
			dtf$subtype <- (dtf[,paste0(gene, '_CNA')] == 1 & dtf$final.ER == 'positive')*1
			}
	} else if (subtype == 'IC10') {
		dtf$subtype <- (dtf$final.ER == 'negative' & dtf$final.HER2 == 'negative')*1
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
		dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$final.ER == 'positive')*1
	}

	# run association 
	fit <- glm(
		subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + snvs,
		data = dtf,
		family = 'binomial'
		)
	ci <- confint(fit)

	out <- data.frame(
		subtype = subtype,
		gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
		coef = coef(fit)[['bds']],
		p = summary(fit)$coefficients['bds',4],
		se = summary(fit)$coefficients['bds',2],
		l95 = ci['bds',1],
		u95 = ci['bds',2],
		number_subtype = sum(dtf$subtype, na.rm = TRUE)
		)
	return(out)
}

run_metabric_subtype_associations <- function(dtf, subtype, gene = NULL) {
	genes <- list(
		IC1 = c('RPS6KB1','DHX40','BCAS3','TUBD1'), 
		IC2 = c('RSF1','CCND1','PAK1','NARS2'), 
		IC6 = c('ZNF703','FGFR1','LETM2'), # EIF4EBP1 no variants
		IC9 = c('MYC','SQLE','FBXO32'),
		IC10 = c('FOXC1','MIA','MELK'),
		IC5 = 'ERBB2'
		)

	if (!is.null(gene)) {
		dtf$bds <- sign(dtf[,gene])
		if (subtype == 'Her2') {
			dtf$subtype <- (dtf$CLAUDIN_SUBTYPE == 'Her2')*1
		} else {
			dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER_IHC == 'pos')
			}
	} else if (subtype == 'IC10') {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
		dtf$subtype <- (dtf$CLAUDIN_SUBTYPE == 'Basal')*1
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
		dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER_IHC == 'pos')

	}

	# run association 
	fit <- glm(
		subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
		data = dtf,
		family = 'binomial'
		)
	ci <- confint(fit)

	out <- data.frame(
		subtype = subtype,
		gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
		coef = coef(fit)[['bds']],
		p = summary(fit)$coefficients['bds',4],
		se = summary(fit)$coefficients['bds',2],
		l95 = ci['bds',1],
		u95 = ci['bds',2],
		number_subtype = sum(dtf$subtype, na.rm = TRUE)
		)
	return(out)
}

### SUBTYPE #######################################################################################
# read in summary data
icgc <- read.delim(
	'icgc_megatable.txt', 
	as.is = TRUE
	)
# test subtype association 
icgc_data_subtype <- rbind(
	run_icgc_subtype_associations(icgc, gene = 'ERBB2', subtype = 'Her2'),
	run_icgc_subtype_associations(icgc, subtype = 'IC1'),
	run_icgc_subtype_associations(icgc, subtype = 'IC2'),
	run_icgc_subtype_associations(icgc, subtype = 'IC9')
	)

# read in summary data
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)
# test subtype association 
metabric_data_subtype <- rbind(
	run_metabric_subtype_associations(metabric, gene = 'ERBB2', subtype = 'Her2'),
	run_metabric_subtype_associations(metabric, subtype = 'IC1'),
	run_metabric_subtype_associations(metabric, subtype = 'IC2'),
	run_metabric_subtype_associations(metabric, subtype = 'IC9')
	)

# add GEL results 
gel_data_subtype <- read.delim('gel_subtype_associations.txt', as.is = TRUE)
gel_data_subtype <- gel_data_subtype[gel_data_subtype$subtype != 'IC10',]

### SUBTYPE PLOT ####
icgc_data_subtype$cohort <- 'ICGC'
metabric_data_subtype$cohort <- 'METABRIC'
gel_data_subtype$cohort <- 'GEL'
plot_data_subtype <- rbind(icgc_data_subtype, metabric_data_subtype, gel_data_subtype)
plot_data_subtype <- plot_data_subtype[order(plot_data_subtype$subtype, plot_data_subtype$cohort),]
plot_data_subtype$index <- rev(c(3,1,2,6,4,5,9,7,8,12,10,11))
plot_data_subtype <- plot_data_subtype[order(plot_data_subtype$index),]
# create forest plot of all cohorts individually
cohort_col <- c("darkorchid4", 'deeppink3',"chartreuse4",'dodgerblue')
names(cohort_col) <- c('GEL','ICGC','DCIS','METABRIC')
cohort_cov <- as.character(plot_data_subtype$cohort)
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
                colours = cohort_col[c('GEL','ICGC','METABRIC')], 
                labels = c('GEL','ICGC','METABRIC'),
                title = 'Cohort',
                border = 'transparent'
                )
        );

cov.legend.grob <- legend.grob(
        legends = cov.legend
        );

create.scatterplot(
        index ~ coef,
        data = plot_data_subtype,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        xat = c(log(c(0.1,0.5)),0, log(c(2,8))),
        xaxis.lab = c('0.1','0.5', '1.0', '2.0', '8.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = '',
        yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = seq(2, nrow(plot_data_subtype), 3),
        ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        add.text = TRUE,
        text.labels = paste0('n=', plot_data_subtype$number_subtype),
        text.x = 1,
        text.y = 1:nrow(plot_data_subtype),
        filename = paste0(date, '_icgc_metabric_gel_scatterplot.pdf'),
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(1,0), x = 0.99, y = 0.01)
                ),
        add.rectangle = TRUE,
        xleft.rectangle = -2.75,
        ybottom.rectangle = c(0, seq(6.5, nrow(plot_data_subtype), 6)),
        xright.rectangle = 2.75,
        ytop.rectangle = c(seq(3.5, nrow(plot_data_subtype), 6)),
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        width = 8,
        height = 5,
        x.error.right = plot_data_subtype$u95-plot_data_subtype$coef,
        x.error.left = plot_data_subtype$coef-plot_data_subtype$l95,
        top.padding = 2,
        resolution = 300
        )