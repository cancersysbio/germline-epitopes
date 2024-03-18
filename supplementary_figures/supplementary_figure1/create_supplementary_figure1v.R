### CREATE SUPPLEMENTARY FIGURE 1V ################################################################
# test keratin and unexpressed genes 
# plot supplementary figure 1V

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST SUBTYPE ASSOCIATION ######################################################################
run_subtype_associations <- function(dtf, subtype, genes_to_test, gene = NULL) {
	# set cna genes
	cna_genes <- c('RPS6KB1','RSF1','ZNF703','MYC')
	names(cna_genes) <- c('IC1','IC2','IC6','IC9')

	if (subtype == 'IC5') {
		if (length(genes_to_test[['IC5']]) > 1) {
			dtf$bds <- rowSums(sign(dtf[,genes_to_test[['IC5']]]))
		} else {
			dtf$bds <- sign(dtf[,genes_to_test[['IC5']]])
			}
		dtf$subtype <- (dtf$pam50 == 'Her2')*1
	} else {
		dtf$bds <- rowSums(sign(dtf[,genes_to_test[[subtype]]]))
		dtf$subtype <- (dtf[,paste0(cna_genes[subtype], '_CNA')] == 1 & dtf$ER == 1)*1
	}

	# run association 
	fit <- glm(
		subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic,
		data = dtf,
		family = 'binomial'
		)
	ci <- confint(fit)

	out <- data.frame(
		subtype = subtype,
		gene = ifelse(!is.null(gene), gene, paste(genes_to_test[[subtype]], collapse = '|')),
		coef = coef(fit)[['bds']],
		p = summary(fit)$coefficients['bds',4],
		se = summary(fit)$coefficients['bds',2],
		l95 = ci['bds',1],
		u95 = ci['bds',2],
		nonzero = length(which(dtf$bds != 0)),
		nonzero_subtype = length(which(dtf$bds != 0 & dtf$subtype == 1))
		)
	return(out)
}

### TEST KERATINS ################################################################################
test_keratins <- function(dtf, keratin, subtype) {
	dtf$bds <- sign(dtf[,keratin])
	dtf$subtype <- (dtf$pam50 == subtype)
	# run association 
	fit <- glm(
		subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic,
		data = dtf,
		family = 'binomial'
		)
	ci <- confint(fit)

	out <- data.frame(
		subtype = subtype,
		gene = keratin,
		coef = coef(fit)[['bds']],
		p = summary(fit)$coefficients['bds',4],
		se = summary(fit)$coefficients['bds',2],
		l95 = ci['bds',1],
		u95 = ci['bds',2],
		nonzero = length(which(dtf$bds != 0)),
		nonzero_subtype = length(which(dtf$bds != 0 & dtf$subtype == 1))
		)
	return(out)
}

### MAIN ##########################################################################################
# read in negative controls 
negcontrols <- read.delim(
	file.path(main_repo_path, 'data', 'controls', 'negative_control_geb_table.txt'),
	as.is = TRUE
	)

# set genes to test
neg_genes <- list(
		IC1 = c('CSH1','CSH2','CSHL1','GH1'),
		IC2 = c('FGF4','CABP2','DEFB108B','OR2AT4'),
		IC6 = c('KCNU1','GOT1L1','ADRB3','UNC5D'),
		IC9 = c('CYP11B1','GML', 'CYP11B2'),
		IC5 = c('CCL1','HNF1B')
		)

# test subtype association 
neg_plot_data <- rbind(
	run_subtype_associations(negcontrols, subtype = 'IC5', genes_to_test = neg_genes),
	run_subtype_associations(negcontrols, subtype = 'IC1', genes_to_test = neg_genes),
	run_subtype_associations(negcontrols, subtype = 'IC2', genes_to_test = neg_genes),
	run_subtype_associations(negcontrols, subtype = 'IC9', genes_to_test = neg_genes)
	)

# read in tcga megatable
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','tcga_megatable.txt'),
	as.is = TRUE
	)

# set genes to test
sub_genes <- list(
                IC1 = c('RPS6KB1','TUBD1','DHX40','BCAS3'),
                IC2 = c('RSF1','CCND1','PAK1','NARS2'),
                IC6 = c('ZNF703','FGFR1','LETM2','EIF4EBP1'),
                IC9 = c('MYC','SQLE','FBXO32'),
                IC5 = 'ERBB2'
                )

# test subtype association 
sub_plot_data <- rbind(
	run_subtype_associations(tcga, subtype = 'IC5', genes_to_test = sub_genes),
	run_subtype_associations(tcga, subtype = 'IC1', genes_to_test = sub_genes),
	run_subtype_associations(tcga, subtype = 'IC2', genes_to_test = sub_genes),
	run_subtype_associations(tcga, subtype = 'IC9', genes_to_test = sub_genes)
	)

# test keratins
keratins <- do.call(rbind, sapply(
	c("KRT71","KRT74","KRT82","KRT34"),
	function(keratin) {
		tmp <- do.call(rbind, sapply(
			c('Her2','LumA','LumB','Basal'),
			test_keratins,
			dtf = negcontrols,
			keratin = keratin,
			simplify = FALSE
			))
		return(tmp)
		},
	simplify = FALSE
	))

# generate plot data
plot_data_control <- data.frame(
	coef = c(keratins$coef, neg_plot_data$coef, sub_plot_data$coef),
	type = c(rep('a', nrow(keratins)), rep('b', nrow(neg_plot_data)), rep('c', nrow(sub_plot_data))),
	subtype = c(keratins$subtype, as.character(neg_plot_data$subtype), as.character(sub_plot_data$subtype))
	)

es_k <- median(plot_data_control[plot_data_control$type == 'c','coef'])-median(plot_data_control[plot_data_control$type == 'a','coef'])
p_k <- wilcox.test(
        plot_data_control[plot_data_control$type == 'c','coef'],
        plot_data_control[plot_data_control$type == 'a','coef']
        )$p.value
p_k_sci <- scientific.notation(p_k, type = 'list')

es_u <- median(plot_data_control[plot_data_control$type == 'c','coef'])-median(plot_data_control[plot_data_control$type == 'b','coef'])
p_u <- wilcox.test(
        plot_data_control[plot_data_control$type == 'c','coef'],
        plot_data_control[plot_data_control$type == 'b','coef']
        )$p.value
p_u_sci <- scientific.notation(p_u, type = 'list')

col <- rep('darkgrey', nrow(plot_data_control))
col[which(plot_data_control$subtype == 'IC1')] <- "darkorange2"
col[which(plot_data_control$subtype == 'IC2')] <- "chartreuse4"
col[which(plot_data_control$subtype == 'IC6')] <- "gold"
col[which(plot_data_control$subtype == 'IC9')] <- "#F8B4E3"
col[which(plot_data_control$subtype %in% c('IC5','Her2'))] <- "#8B0000"
col[which(plot_data_control$type == 'a')] <- "darkgrey"

create.boxplot(
        coef ~ type,
        data = plot_data_control,
        add.stripplot = TRUE,
        filename = paste0(date, '_supplementary_figure1v.png'),
        resolution = 300,
        ylimits = c(-1.1,1.1),
        yat = seq(-1, 1, 0.5),
        points.col = col,
        points.cex = 1,
        ylab.label = 'Coefficient',
        xlab.label = 'Association',
        xaxis.lab = c('Non-expressed\nkeratins','Non-expressed\nsubtype-specific\nproteins','Expressed\nsubtype-specific\nproteins'),
        xaxis.cex = 1,
        legend = list(
	        inside = list(
	                 fun = draw.key,
	                 args = list(
	                     key = list(
	                         text = list(
	                             lab = c(
	                                paste('ES:', round(es_k, digits = 2)),
	                                as.expression(substitute(
	                                        base *' x '* 10^exponent, 
	                                        list(base = paste('P:', p_k_sci[[1]]), exponent = p_k_sci[[2]])
	                                        ))
	                                )
	                             ),
	                         cex = 1.5
	                         )
	                     ),
	                 x = 0.01,
	                 y = 0.01,
	                 corner = c(0,0),
	                 draw = FALSE
	                 ),
	        inside = list(
	                 fun = draw.key,
	                 args = list(
	                     key = list(
	                         text = list(
	                             lab = c(
	                                paste('ES:', round(es_u, digits = 2)),
	                                as.expression(substitute(
	                                        base *' x '* 10^exponent, 
	                                        list(base = paste('P:', p_u_sci[[1]]), exponent = p_u_sci[[2]])
	                                        ))
	                                )
	                             ),
	                         cex = 1.5
	                         )
	                     ),
	                 x = 0.99,
	                 y = 0.99,
	                 corner = c(1,1),
	                 draw = FALSE
	                 )
	             )
        )


