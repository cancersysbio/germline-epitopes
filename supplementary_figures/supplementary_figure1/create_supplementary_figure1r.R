### CREATE SUPPLEMENTARY FIGURE 1R ################################################################
# create summary of gwas summary stats for snps of interest

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in snps 
snps <- read.delim(
	file.path(main_repo_path, 'data', '/auxiliary_data', 'gwas_snps_summary_stats.txt'),
	as.is = TRUE
	)
snps <- snps[order(snps$beta.Gwas),]
snps$index <- 1:nrow(snps)

# add genes as covariates 
gene.covs <- as.character(snps$gene)
gene.covs[which(gene.covs == 'RSF1')] <- default.colours(7)[1]
gene.covs[which(gene.covs == 'ERBB2')] <- default.colours(7)[2]
gene.covs[which(gene.covs == 'ZNF703')] <- default.colours(7)[3]
gene.covs[which(gene.covs == 'MYC')] <- default.colours(7)[4]
gene.covs[which(gene.covs == 'BCAS3')] <- default.colours(7)[5]
gene.covs[which(gene.covs == 'TUBD1')] <- default.colours(7)[6]
gene.covs[which(gene.covs == 'SQLE')] <- default.colours(7)[7]

covariates.object <- list(
         rect = list(
             col = 'white',
             fill = gene.covs,
             lwd = 1.5
             )
        )
covariate.object.grob <- covariates.grob(
         covariates = covariates.object,
         ord = c(1:nrow(snps)),
         side = 'top'
         );

covariates.legends <- list(
         legend = list(
             colours = default.colours(7),
             labels = c('RSF1','ERBB2','ZNF703','MYC','BCAS3','TUBD1','SQLE'),
             title = 'Genes',
             border = 'white'
             )
         )
covariate.legend.grob <- legend.grob(
         legends = covariates.legends,
         title.just = 'left'
         );

p1 <- create.barplot(
	beta.Gwas ~ index,
	data = snps,
	xaxis.lab = rep('', nrow(snps)),
	xlab.label = '',
	ylab.label = 'Beta',
	ylimits = c(-1.5, 1.5),
	key = NULL,
	legend = list(
             top = list(fun = covariate.object.grob),
             right = list(fun = covariate.legend.grob)
             ),
	resolution = 300
	)
p2 <- create.barplot(
	-log10(P.value.Gwas) ~ index,
	data = snps,
	xaxis.lab = paste0('chr', gsub('_',':', snps$snp)),
	ylimits = c(0,1.5),
	yat = c(0, 0.3, 1.3),
	abline.h = -log10(0.05),
	abline.lty = 2,
	yaxis.lab = c('1','0.5','0.05'),
	ylab.label = 'P-value',
	xaxis.rot = 45,
	xlab.label = 'SNP',
	resolution = 300
	)
create.multipanelplot(
	list(p1, p2),
	plot.objects.heights = c(0.35,0.65),
	file = paste0(date, '_supplementary_figure1r.png'),
	resolution = 300,
	width = 9,
	height = 7
	)