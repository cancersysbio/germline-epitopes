### CREATE SUPPLEMENTARY FIGURE 1I ################################################################
# create forest plot of associations with varying binder thresholds
# supplementary figure 1i

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# create thresholds plot 
her2_thres <- read.delim(
	file.path(main_repo_path, 'data', 'controls', 'tcga_her2_thresholds.txt'), 
	as.is = TRUE
	)
res <- list()
for (thres in c('wbs_0.5_2','wbs_0.5_3','wbs_0.5_4','wbs_0.25_2','wbs_0.25_3','wbs_0.25_4')) {
	tmp <- her2_thres
	tmp$bds <- sign(tmp[,thres])
	tmp$subtype <- (tmp$pam50 == 'Her2')*1
	fit <- glm(
		subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
		data = tmp,
		family = 'binomial'
		)
	ci <- confint(fit)
	res[[thres]] <- data.frame(
		threshold = thres,
		coef = coef(fit)[['bds']],
		l95 = ci['bds',1],
		u95 = ci['bds',2],
		p = summary(fit)$coefficients['bds',4]
		)
}
res <- do.call(rbind, res)
res$index <- c(4,5,6,1,2,3)

create.scatterplot(
        index ~ coef,
        data = res,
        horizontal = TRUE,
        xlimits = c(-1.5,1.5),
        xat = log(c(0.15, 0.5, 1, 2.5, 8)),
        xaxis.lab = c('0.15','0.50','1.00','2.50','8.00'),
        filename = paste0(date, '_her2_thresholds_scatterplot.png'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Rank Thresholds',
        ylimits = c(0.5, nrow(res)+0.5),
        yaxis.lab = c('0.25-2','0.25-3','0.25-4','0.5-2','0.5-3','0.5-4'),
        yat = 1:nrow(res),
        abline.v = 0,
        width = 10,
        height = 4,
        top.padding = 2,
        x.error.right = res$u95-res$coef,
        x.error.left = res$coef-res$l95,
        resolution = 300
        )