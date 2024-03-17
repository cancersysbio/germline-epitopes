### SUPPLEMENTARY FIGURE 1G #######################################################################
# test all definitions of her2 in TCGA
# create supplementary figure 1G

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### TEST ASSOCIATIONS #############################################################################
test_association <- function(tcga, subtype) {
        fit <- glm(
                as.formula(paste(subtype, '~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic')),
                data = tcga,
                family = 'binomial'
                )
        ci <- confint(fit)
        res <- data.frame(
                subtype = subtype,
                coef = coef(fit)[['bds']],
                p = summary(fit)$coefficients['bds',4],
                l95 = ci['bds',1],
                u95 = ci['bds',2],
                num = sum(tcga[,subtype] == 1)
                )
        return(res)
}

### TCGA ##########################################################################################
# read in summary data
tcga <- read.delim(
        file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
	as.is = TRUE
	)
tcga$bds <- sign(tcga$ERBB2)

tcga$pam50_subtype <- (tcga$pam50 == 'Her2')*1
tcga$cna_subtype <- (tcga$ERBB2_CNA_5copies == 1)*1
tcga$ihc_subtype <- (tcga$HER2.newly.derived == 'Positive')*1
tcga$ihc_erneg_subtype <- (tcga$HER2.newly.derived == 'Positive' & tcga$ER == 0)*1
tcga$ihc_erpos_subtype <- (tcga$HER2.newly.derived == 'Positive' & tcga$ER == 1)*1
tcga$cna_erneg_subtype <- (tcga$ERBB2_CNA_5copies == 1 & tcga$ER == 0)*1

# read associations
tplot_data <- rbind(
	test_association(tcga, 'pam50_subtype'),
	test_association(tcga, 'cna_subtype'),
	test_association(tcga, 'ihc_subtype'),
	test_association(tcga, 'ihc_erneg_subtype'),
        test_association(tcga, 'ihc_erpos_subtype'),
        test_association(tcga, 'cna_erneg_subtype')
	)
tplot_data <- tplot_data[order(-tplot_data$coef),]
tplot_data$index <- 1:nrow(tplot_data)

# create forest plot
create.scatterplot(
        index ~ coef,
        data = tplot_data,
        horizontal = TRUE,
        xlimits = c(-2.5,2.5),
        filename = paste0(date, '_TCGA_her2_definitions_scatterplot.png'),
        xat = log(c(0.2,0.5, 1, 2, 5)),
        xaxis.lab = c('0.2','0.5', '1.0', '2.0', '5.0'),
        xlab.label = 'Odds Ratio',
        ylab.label = 'HER2 Definition',
        yaxis.lab = rev(c(
        	'PAM50 HER2+','CNA (>5) HER2+/ER-','CNA (>5) HER2+','Clinical HER2+/ER-','Clinical HER2+', 'Clinical HER2+/ER+'
        	)), 
        yat = 1:nrow(tplot_data),
        ylimits = c(0.5, nrow(tplot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        add.text = TRUE,
        text.labels = paste0('n=', tplot_data$num),
        text.y = 1:nrow(tplot_data),
        text.x = rep(1.65,nrow(tplot_data)),
        text.cex = 1.2,
        x.error.right = tplot_data$u95-tplot_data$coef,
        x.error.left = tplot_data$coef-tplot_data$l95,
        height = 3.5,
        width = 8.5,
        top.padding = 2,
        left.padding = 5,
        resolution = 300
        )