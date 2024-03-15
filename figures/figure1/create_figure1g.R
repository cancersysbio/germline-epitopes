### CREATE FIGURE1G ###############################################################################
# create figure 1g

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

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
                IC10 = c('FOXC1','MIA','MELK','KNTC2')
                )

        if (!is.null(gene)) {
                dtf$bds <- sign(dtf[,gene])
                if (subtype == 'Her2') {
                        dtf$subtype <- (dtf$pam50 == 'Her2')*1
                } else {
                        dtf$subtype <- (dtf[,paste0(gene, '_CNA')] == 1 & dtf$ER == 1)*1
                        }
        } else if (subtype == 'IC10') {
                dtf$subtype <- (dtf$pam50 == 'Basal')*1
                dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
        } else {
                dtf$bds <- rowSums(sign(dtf[,genes[[subtype]]]))
                dtf$subtype <- (dtf[,paste0(genes[[subtype]][1], '_CNA')] == 1 & dtf$ER == 1)*1
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
                gene = ifelse(!is.null(gene), gene, paste(genes[[subtype]], collapse = '|')),
                coef = coef(fit)[['bds']],
                p = summary(fit)$coefficients['bds',4],
                se = summary(fit)$coefficients['bds',2],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                number_subtype = sum(dtf$subtype)
                )
        return(out)
}

### MAIN #######################################################################################
# read in summary data
tcga <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','tcga_megatable.txt'), 
	as.is = TRUE
	)
# test subtype association 
plot_data_subtype <- rbind(
	run_subtype_associations(tcga, gene = 'ERBB2', subtype = 'Her2'),
	run_subtype_associations(tcga, subtype = 'IC1'),
	run_subtype_associations(tcga, subtype = 'IC2'),
	run_subtype_associations(tcga, subtype = 'IC9')
	)

### SUBTYPE PLOT ####
plot_data_subtype <- plot_data_subtype[plot_data_subtype$nonzero_subtype > 6,]
plot_data_subtype$index <- 1:nrow(plot_data_subtype)
plot_data_subtype <- plot_data_subtype[order(plot_data_subtype$index),]

create.scatterplot(
        index ~ coef,
        data = plot_data_subtype,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        filename = paste0(date, '_TCGA_only_scatterplot.pdf'),
        xat = log(c(0.2, 0.5, 1, 2, 5)),
        xaxis.lab = c('0.2','0.5','1.0','2.0','5.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Subtype',
        yaxis.lab = gsub('IC10', 'TNBC', gsub('Her2','HER2+', unique(plot_data_subtype$subtype))),
        yat = 1:nrow(plot_data_subtype),
        ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        add.text = TRUE,
        text.labels = plot_data_subtype$gene,
        text.y = 1:nrow(plot_data_subtype),
        text.x = rep(1.65,nrow(plot_data_subtype)),
        text.cex = 0.8,
        x.error.right = plot_data_subtype$u95-plot_data_subtype$coef,
        x.error.left = plot_data_subtype$coef-plot_data_subtype$l95,
        width = 8.5,
        height = 3,
        top.padding = 2,
        left.padding = 5,
        resolution = 300
        )