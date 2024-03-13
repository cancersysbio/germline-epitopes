### CREATE SUPPLEMENTARY FIGURE 4IJK ##############################################################
# create supplementary figure 4ijk


### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general) 

### MAIN ##########################################################################################
# read in summary data
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)
# group IC4
metabric[metabric$INTCLUST %in% c('4ER+','4ER-'),'INTCLUST'] <- 4

## HER2 ###########################################################################################
meta_her2 <- metabric[which(metabric$CLAUDIN_SUBTYPE == 'Her2'),]
meta_her2$bds <- (sign(meta_her2$ERBB2) >= 0)*1
meta_her2$ie <- (meta_her2$TME_subtype %in% c('IE','IE/F'))*1

fit <- glm(ie ~ bds + PC1 + PC2 + PC3 + INTCLUST + stage + age, data = meta_her2, family = 'binomial')
or <- round(exp(coef(fit)[['bds']]), digits = 2)
p <- round(summary(fit)$coefficients['bds',4], digits = 2)

# create plot 
plot_data <- as.data.frame(table(meta_her2[,c('bds','ie')]))
plot_data$prop <- NA
plot_data[plot_data$bds == 0,'prop'] <- plot_data[plot_data$bds == 0,'Freq']/sum(plot_data[plot_data$bds == 0,'Freq'])
plot_data[plot_data$bds == 1,'prop'] <- plot_data[plot_data$bds == 1,'Freq']/sum(plot_data[plot_data$bds == 1,'Freq'])

create.barplot(
	prop ~ bds,
	data = plot_data,
	groups = plot_data$ie,
	xaxis.lab = c('Low','High'),
	xlab.label = 'GEB',
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	main = paste0('METABRIC HER2+\n(n=', nrow(meta_her2), ')'),
	main.cex = 1.5,
	ylab.label = 'Proportion of Samples',
	filename = paste0(date, '_METABRIC_HER2_TME_barplot.pdf'),
	stack = TRUE,
	key = NULL,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('grey','black')
                             ),
                         text = list(
                             lab = c('Immune Depleted','Immune Enriched')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     ),
                 x = -0.2,
                 y = 1.175
                 )
             ),
	add.text = TRUE,
	text.labels = c(paste0('OR=',or), paste0('P=',p)),
	text.x = c(1,1),
	text.y = c(0.1,0.05),
	text.col = 'black',
	text.cex = 1.2,
	text.fontface = 'bold',
	top.padding = 10,
	width = 6,
	height = 7,
	col = c('grey','black'),
	resolution = 300
	)

## ER #############################################################################################
meta_er <- metabric[which(metabric$INTCLUST %in% c(1,2,6,9) & metabric$stage > 1),]
meta_er$bds <- sign(colSums(rbind(
	rowSums(sign(meta_er[,c('MYC','SQLE','FBXO32')]))*meta_er$MYC_CNA,
	rowSums(sign(meta_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_er$RPS6KB1_CNA,
	rowSums(sign(meta_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_er$RSF1_CNA,
	rowSums(sign(meta_er[,c('ZNF703','FGFR1','LETM2')]))*meta_er$ZNF703_CNA
	)))
meta_er$ie <- (meta_er$TME_subtype %in% c('IE','IE/F'))*1
meta_er$bds <- (meta_er$bds > 0)*1

fit <- glm(ie ~ bds + PC1 + PC2 + PC3 + INTCLUST + stage + age, data = meta_er, family = 'binomial')
or <- round(exp(coef(fit)[['bds']]), digits = 2)
p <- round(summary(fit)$coefficients['bds',4], digits = 2)

# create plot 
plot_data <- as.data.frame(table(meta_er[,c('bds','ie')]))
plot_data$prop <- NA
plot_data[plot_data$bds == 0,'prop'] <- plot_data[plot_data$bds == 0,'Freq']/sum(plot_data[plot_data$bds == 0,'Freq'])
plot_data[plot_data$bds == 1,'prop'] <- plot_data[plot_data$bds == 1,'Freq']/sum(plot_data[plot_data$bds == 1,'Freq'])

create.barplot(
	prop ~ bds,
	data = plot_data,
	groups = plot_data$ie,
	xaxis.lab = c('Low','High'),
	xlab.label = 'GEB',
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	main = paste0('METABRIC ER+\n(n=', nrow(meta_er), ')'),
	main.cex = 1.5,
	ylab.label = 'Proportion of Samples',
	filename = paste0(date, '_METABRIC_ER_TME_barplot.pdf'),
	stack = TRUE,
	key = NULL,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('grey','black')
                             ),
                         text = list(
                             lab = c('Immune Depleted','Immune Enriched')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     ),
                 x = -0.2,
                 y = 1.175
                 )
             ),
	add.text = TRUE,
	text.labels = c(paste0('OR=',or), paste0('P=',p)),
	text.x = c(1,1),
	text.y = c(0.1,0.05),
	text.col = 'black',
	text.cex = 1.2,
	text.fontface = 'bold',
	top.padding = 10,
	width = 6,
	height = 7,
	col = c('grey','black'),
	resolution = 300
	)

### SUPPLEMENTARY FIGURE 4K #######################################################################
res_her2 <- list()
for (i in c('T_cell_traffic','Immune')) {
	# test her2
	fit_her2 <- lm(as.formula(paste(i, '~ bds + PC1 + PC2 + PC3 + INTCLUST + stage + age')), data = meta_her2)
	ci_her2 <- confint(fit_her2)
	res_her2[[i]] <- data.frame(
		feature = i,
		subtype = 'HER2',
		coef  = coef(fit_her2)[['bds']], 
		p = summary(fit_her2)$coefficients['bds',4],
		l95 = ci_her2['bds',1], 
		u95 = ci_her2['bds',2]
		)
} 
res_her2 <- do.call(rbind, res_her2)

res_er <- list()
for (i in c('T_cell_traffic','Immune')) {
	# test her2
	fit_er <- lm(as.formula(paste(i, '~ bds + PC1 + PC2 + PC3 + INTCLUST + stage + age')), data = meta_er)
	ci_er <- confint(fit_er)
	res_er[[i]] <- data.frame(
		feature = i,
		subtype = 'ER',
		coef  = coef(fit_er)[['bds']], 
		p = summary(fit_er)$coefficients['bds',4],
		l95 = ci_er['bds',1], 
		u95 = ci_er['bds',2]
		)
} 
res_er <- do.call(rbind, res_er)
# create plot_data
plot_data <- rbind(res_er,res_her2)
plot_data$index <- 1:nrow(plot_data)

cov <- list(
        rect = list(
                col = 'transparent',
                fill = rep(c('dodgerblue3', 'slateblue4'), 2)
                )
        );

cov.grob <- covariates.grob(
        covariates = cov,
        ord = c(1:nrow(plot_data_subtype)),
        side = 'right',
        size = 1
        );

cov.legend <- list(
        legend = list(
                colours =  c('dodgerblue3', 'slateblue4'),
                labels = c('Immune cells','T cells traffic'),
                title = 'Signature',
                border = 'transparent'
                )
        );

cov.legend.grob <- legend.grob(
        legends = cov.legend
        );


create.scatterplot(
        index ~ coef,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-0.8,0.8),
        filename = paste0(date, '_metabric_immune_sig_scatterplot.pdf'),
        xlab.label = 'Coefficient',
        ylab.label = '',
        yaxis.lab = c(paste0('ER+ \n(n=', nrow(meta_er), ')'),paste0('HER2+ \n(n=', nrow(meta_her2), ')')),
        yat = c(1.5,3.5),
        ylimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        x.error.right = plot_data$u95-plot_data$coef,
        x.error.left = plot_data$coef-plot_data$l95,
        add.rectangle = TRUE,
        xleft.rectangle = -2.5,
        ybottom.rectangle = 0,
        xright.rectangle = 2.5,
        ytop.rectangle = 2.5,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
        legend = list(
                right = list(fun = cov.grob),
                inside = list(fun = cov.legend.grob, corner = c(1,0), x = 0.99, y = 0.05)
                ),
        width = 8,
        height = 3,
        top.padding = 2,
        resolution = 300
        )