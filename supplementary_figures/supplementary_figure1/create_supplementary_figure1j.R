### CREATE SUPPLEMENTARY FIGURE 1J ################################################################
# create forest plot of mhcflurry results
# create supplementary figure1j 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in summary tcga file
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
	as.is = TRUE
	)

# read in mhcflurry results
res_avg <- read.delim(
	file.path(main_repo_path, 'data', 'controls', 'tcga_her2_mhcflurry.txt'),
	as.is = TRUE
	)

# merge with summary tcga
tcga <- merge(res_avg, tcga, by = 'sample')
tcga$subtype <- (tcga$pam50 == 'Her2')*1

flurryres <- list()
for (i in colnames(res_avg)[-1]) {
	fit <- glm(
		as.formula(paste('subtype ~ ', i, '+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6')), 
		data = tcga, family = 'binomial'
		)
	ci <- confint(fit)
	flurryres[[i]] <- data.frame(
		thres = i,
		coef = coef(fit)[[i]],
		l95 = ci[i,1],
		u95 = ci[i,2],
		p = summary(fit)$coefficients[i,4]
		)
}
flurryres <- do.call(rbind, flurryres)
flurryres <- flurryres[flurryres$thres %in% c('wbs_0.25_5','wbs_0.25_3','wbs_0.25_2','wbs_0.5_5','wbs_0.5_3','wbs_0.5_2'),]
flurryres$index <- 1:nrow(flurryres)

create.scatterplot(
        index ~ coef,
        data = flurryres,
        horizontal = TRUE,
        xlimits = c(-2,2),
        xat = log(c(0.15, 0.5, 1, 2.5, 8)),
        xaxis.lab = c('0.15','0.50','1.00','2.50','8.00'),
        filename = paste0(date, '_supplementary_figure1j.png'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Rank Thresholds',
        ylimits = c(0.5, nrow(flurryres)+0.5),
		yaxis.lab = c('0.25-2','0.25-3','0.25-5','0.5-2','0.5-3','0.5-5'),
        yat = 1:nrow(flurryres),
        abline.v = 0,
        width = 10,
        height = 4,
        top.padding = 2,
        x.error.right = flurryres$u95-flurryres$coef,
        x.error.left = flurryres$coef-flurryres$l95,
        resolution = 300
        )