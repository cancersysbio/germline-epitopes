### CREATE FIGURE 1H ###############################################################################
# create figure 1H

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

meta_plot_data <- list()
for (subtype in unique(plot_data_subtype$subtype)) {
	meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                        data = plot_data_subtype[plot_data_subtype$subtype == subtype,])
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
meta_plot_data <- meta_plot_data[order(meta_plot_data$index),]

# create plot
create.scatterplot(
        index ~ coef,
        data = meta_plot_data,
        horizontal = TRUE,
        xlimits = c(-1.4, 1.4),
        filename = paste0(date, '_meta_scatterplot.pdf'),
         xat = c(log(c(0.2,0.5)),0, log(c(2,6))),
        xaxis.lab = c('0.2','0.5', '1.0', '2.0', '6.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('IC10', 'TNBC', gsub('Her2','HER2+', unique(meta_plot_data$subtype))),
        yat = 1:nrow(meta_plot_data),
        ylimits = c(0.5, nrow(meta_plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        x.error.right = meta_plot_data$u95-meta_plot_data$coef,
        x.error.left = meta_plot_data$coef-meta_plot_data$l95,
        width = 8.5,
        height = 3,
        top.padding = 2,
        resolution = 300
        )