### CREATE SUPPLEMENTARY FIGURE 1Q ################################################################
# create supplementary figures 1q

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### CREATE CONTINGENCY MULTIPLOT ##################################################################
create_contingency_multiplot <- function(df, filename, ylimits = c(0,0.22), yat = seq(0,0.2,0.05),
        ylab.label = 'Ratio of HER2+/HER2-\n', text.y = c(0.19,0.18)) {
        df_tmp <- df
        df_tmp$bds <- sign(df_tmp$bds)
        cont_table <- table(df_tmp[,c('bds','subtype')])
        bplot_data <- data.frame(
                ratio = c(
                        cont_table[1,2]/cont_table[1,1],
                        cont_table[2,2]/cont_table[2,1],
                        cont_table[3,2]/cont_table[3,1]
                        ),
                bds = c(-1,0,1)
                )
        fit <- glm(subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic, data = df, family = 'binomial')
        or <- exp(coef(fit)[['bds']])
        pvalue_sci <- scientific.notation(summary(fit)$coefficients['bds',4], digits = 2,type = 'list');
        main1 <- as.expression(substitute(
                                        base, 
                                        list(base = paste0('OR = ', round(or, digits = 2)))
                                        ))
        main2 <- as.expression(substitute(
                                        base *' x '* 10^exponent, 
                                        list(base = paste0('P = ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
                                        ))

                
        # create barplot
        create.barplot(
                ratio ~ bds,
                ylimits = ylimits,
                yat = yat,
                filename = filename,
                xaxis.lab = c('<WT','WT','>WT'),
                ylab.label = ylab.label,
                xlab.label = 'Epitope Burden',
                data = bplot_data,
                #main = main,
                #main.cex = 2,
                add.text = TRUE,
                text.labels = c(main1, main2),
                text.x = c(2.9, 3.05),
                text.y = text.y,
                #text.y = c(0.145,0.135),
                text.cex = 1.2,
                resolution = 300 
                )
}

### MAIN #######################################################################################
# read in summary data
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables','tcga_megatable.txt'), 
	as.is = TRUE
	)

# create supplementary figure 1q
tcga$bds <- sign(tcga$ZNF703)+sign(tcga$FGFR1)+sign(tcga$LETM2)+sign(tcga$EIF4EBP1)
tcga$subtype <- (tcga$ZNF703_CNA == 1 & tcga$ER == 1)*1
create_contingency_multiplot(tcga, filename = paste0(date, '_supplementary_figure1q.png'),
        ylimits = c(0,0.3), ylab.label = 'Ratio of IC6/Other\n', yat = seq(0,0.3,0.1),
        text.y = c(0.29,0.27))


