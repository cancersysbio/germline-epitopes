### CREATE SUPPLEMENTARY FIGURE2I #################################################################
# create supplementary figure 2i
# replicate TNBC associations in ICGC and GEL

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

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

### SUBTYPE #######################################################################################
# read in summary data
icgc <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'icgc_megatable.txt'), 
	as.is = TRUE
	)
# test subtype association 
icgc_data_subtype <- run_icgc_subtype_associations(icgc, subtype = 'IC10')
icgc_data_subtype$cohort <- 'ICGC'

# add GEL results 
gel_data_subtype <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','gel_subtype_associations.txt')
	, as.is = TRUE
	)
gel_data_subtype <- gel_data_subtype[gel_data_subtype$subtype == 'IC10',]
gel_data_subtype$cohort <- 'GEL'

### SUBTYPE PLOT ####
plot_data_tnbc <- rbind(icgc_data_subtype, gel_data_subtype)

meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                        data = plot_data_tnbc)
res <- rma(yi, vi, data=meta_data)
meta_plot_data <- data.frame(
                        subtype = 'IC10',
                        coef = res$beta[[1]],
                        l95 = res$ci.lb,
                        u95 = res$ci.ub,
                        p = res$pval
                        )
meta_plot_data$cohort <- 'meta-analysis'

plot_data_tnbc <- rbind(plot_data_tnbc[,c('subtype','coef','l95','u95','p','cohort')], meta_plot_data)
plot_data_tnbc$index <- 1:nrow(plot_data_tnbc)


create.scatterplot(
        index ~ coef,
        data = plot_data_tnbc,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        xat = c(log(c(0.1,0.5)),0, log(c(2,8))),
        xaxis.lab = c('0.1','0.5', '1.0', '2.0', '8.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = '',
        yaxis.lab = c('ICGC','GEL','Meta-analysis'),
   	yat = 1:nrow(plot_data_tnbc),
        ylimits = c(0.5, nrow(plot_data_tnbc)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        filename = paste0(date, '_supplementary_figure2i.png'),
        width = 8,
        height = 3,
        add.rectangle = TRUE,
        xleft.rectangle = -2.75,
        ybottom.rectangle = 2.5,
        xright.rectangle = 2.75,
        ytop.rectangle = 4,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        x.error.right = plot_data_tnbc$u95-plot_data_tnbc$coef,
        x.error.left = plot_data_tnbc$coef-plot_data_tnbc$l95,
        top.padding = 2,
        resolution = 300
        )
