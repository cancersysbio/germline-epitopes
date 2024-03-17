### CREATE SUPPLEMENTARY FIGURE2L ##################################################################
# create supplementary figure 2l
# GEB associations in DCIS

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST SUBTYPE ASSOCIATION ######################################################################
run_dcis_subtype_associations <- function(dtf, subtype, gene = NULL) {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32'),
		IC5 = 'ERBB2'
		)

	if (!is.null(gene)) {
		dtf$bds <- sign(dtf[,gene])
		if (subtype == 'Her2') {
			dtf$subtype <- (dtf$PAM50 == 'Her2')*1
		} else {
			dtf$subtype <- (dtf[,'IC11_RNA'] == gsub('IC','',subtype))*1
			}
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
		dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER_RNA == '+')*1
	}

	# run association 
	fit <- glm(
		subtype ~ bds + EA_PC1 + EA_PC2 + AA_PC1 + NA_PC1 + NA_PC2 + NA_PC3 + CO_PC1 + CO_PC2,
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

### MAIN ##########################################################################################
# read in summary data
dcis <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'dcis_megatable.txt'), 
	as.is = TRUE
	)
# test subtype association 
dcis_data_subtype <- rbind(
	run_dcis_subtype_associations(dcis, gene = 'ERBB2', subtype = 'Her2'),
	run_dcis_subtype_associations(dcis, subtype = 'IC1'),
	run_dcis_subtype_associations(dcis, subtype = 'IC2')
	)
dcis_data_subtype$index <- 1:nrow(dcis_data_subtype)

create.scatterplot(
        index ~ coef,
        data = dcis_data_subtype,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        xat = c(log(c(0.1,0.5)),0, log(c(2,8))),
        xaxis.lab = c('0.1','0.5', '1.0', '2.0', '8.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = '',
        yaxis.lab = gsub('Her2','HER2+', unique(dcis_data_subtype$subtype)),
   		yat = 1:nrow(dcis_data_subtype),
        ylimits = c(0.5, nrow(dcis_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        filename = paste0(date, '_dcis_scatterplot.png'),
        width = 8,
        height = 3,
        x.error.right = dcis_data_subtype$u95-dcis_data_subtype$coef,
        x.error.left = dcis_data_subtype$coef-dcis_data_subtype$l95,
        top.padding = 2,
        resolution = 300
        )