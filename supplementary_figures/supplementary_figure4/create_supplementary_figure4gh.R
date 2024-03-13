### CREATE SUPPLEMENTARY FIGURE 4GH ###############################################################
# create supplementary figure 4gh

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general) 

### MAIN ##########################################################################################
# read in summary data
tcga <- read.delim(
	'tcga_megatable.txt', 
	as.is = TRUE
	)

## HER2 ## 
meta_her2 <- tcga[which(tcga$pam50 == 'Her2'),]
meta_her2$bds <- (sign(meta_her2$ERBB2) >= 0)*1
meta_her2$ie <- (meta_her2$TME_subtype %in% c('IE','IE/F'))*1

fit <- glm(ie ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + ic10 + age + stage,
	data = meta_her2, family = 'binomial')
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
	main = paste0('TCGA HER2+\n(n=', nrow(meta_her2), ')'),
	main.cex = 1.5,
	ylab.label = 'Proportion of Samples',
	filename = paste0(date, '_HER2_TME_barplot.pdf'),
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

## ER ## 
meta_er <- tcga[which(tcga$ic10 %in% c(1,2,6,9)),]
meta_er$bds <- sign(colSums(rbind(
	rowSums(sign(meta_er[,c('MYC','SQLE','FBXO32')]))*meta_er$MYC_CNA,
	rowSums(sign(meta_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_er$RPS6KB1_CNA,
	rowSums(sign(meta_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_er$RSF1_CNA,
	rowSums(sign(meta_er[,c('ZNF703','FGFR1','LETM2')]))*meta_er$ZNF703_CNA
	)))
meta_er$ie <- (meta_er$TME_subtype %in% c('IE','IE/F'))*1
meta_er$bds <- (meta_er$bds > 0)*1

fit <- glm(ie ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + ic10 + age + stage, 
	data = meta_er, family = 'binomial')
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
	main = paste0('TCGA ER+\n(n=', nrow(meta_er), ')'),
	main.cex = 1.5,
	ylab.label = 'Proportion of Samples',
	filename = paste0(date, '_ER_TME_barplot.pdf'),
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


