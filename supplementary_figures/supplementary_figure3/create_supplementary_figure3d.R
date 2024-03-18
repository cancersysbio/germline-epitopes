### CREATE SUPPLEMENTARY FIGURE 3D #################################################################
# create forest plot of subtype associations in hartwig
# create supplementary figure 3d

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
run_subtype_associations <- function(dtf, subtype, gene = NULL) {
	genes <- list(
		IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
		IC2 = c('RSF1','CCND1','PAK1','NARS2'),
		IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
		IC9 = c('MYC','SQLE','FBXO32'),
		IC5 = 'ERBB2'
		)

	if (!is.null(gene)) {
		dtf$bds <- sign(dtf[,gene])
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
	}

	if (!is.null(gene)) {
		dtf$bds <- sign(dtf[,gene])
		if (gene == 'ERBB2') {
                        dtf$subtype <- (dtf$hormone == 'ER-negative/HER2-positive')*1
		} else {
			dtf$subtype <- (dtf[,paste0(gene, '_CNA')] == 1 & dtf$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive'))*1
			}
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))

		dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$hormone %in% c('ER-positive/HER2-negative','ER-positive/HER2-positive'))*1
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
		nonzero = sum(dtf$bds != 0),
		num_subtype = sum(dtf$subtype)
		)
	return(out)
}

### MAIN ##########################################################################################
# read in summary data
hartwig <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','hartwig_megatable.txt'), 
	as.is = TRUE
	)

# test subtype association 
plot_data_subtype <- rbind(
	run_subtype_associations(hartwig, gene = 'ERBB2', subtype = 'Her2'),
	run_subtype_associations(hartwig, subtype = 'IC1'),
	run_subtype_associations(hartwig, subtype = 'IC2'),
	run_subtype_associations(hartwig, subtype = 'IC9')
	)
plot_data_subtype$index <- 1:nrow(plot_data_subtype)

create.scatterplot(
        index ~ coef,
        data = plot_data_subtype,
        horizontal = TRUE,
        xlimits = c(-2,2),
        filename = paste0(date, '_supplementary_figure3d.png'),
        xat = c(log(c(0.2,0.5)),0, log(c(2,6))),
        xaxis.lab = c('0.2','0.5', '1.0', '2.0', '6.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = 1:nrow(plot_data_subtype),
        ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        x.error.right = plot_data_subtype$u95-plot_data_subtype$coef,
        x.error.left = plot_data_subtype$coef-plot_data_subtype$l95,
        width = 8,
        height = 3,
        top.padding = 2,
        resolution = 300
        )
