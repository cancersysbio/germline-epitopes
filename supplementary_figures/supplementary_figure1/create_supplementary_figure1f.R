### SUPPLEMENTARY FIGURE 1F #######################################################################
# test continuous her2 association

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

date <- Sys.Date()
### TCGA ##########################################################################################
# read in summary data
tcga <- read.delim(
        'tcga_megatable.txt',
	as.is = TRUE
	)
tcga$bds <- tcga$ERBB2
tcga$subtype <- (tcga$pam50 == 'Her2')*1

fit <- glm(
        subtype ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + somatic,
        data = tcga,
        family = 'binomial'
        )
beta <- round(coef(fit)[['bds']], digits = 2)

tcga$subtype <- factor(tcga$subtype)
# create boxplot
create.boxplot(
        bds ~ subtype,
        data = tcga,
        add.stripplot = TRUE,
        filename = paste0(date, '_HER2_boxplot.png'),
        ylab.label = 'Average Binders',
        xlab.label = 'Subtype',
        xaxis.lab = c('HER2-','HER2+'),
        ylimits = c(-1.2, 1.2),
        yat = seq(-1, 1, 0.5),
        legend = list(
                inside = list(
                         fun = draw.key,
                         args = list(
                             key = list(
                                 text = list(
                                     lab = c(
                                        expression(beta*'= -1.62'),
                                        as.expression(substitute(
                                                base *' x '* 10^exponent, 
                                                list(base = paste('P=', '2.77'), exponent = '-3')
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
                ),
        resolution = 300
        )
