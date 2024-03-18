### CREATE SUPPLEMENTARY FIGURE 4A-F,M-N ##########################################################
# create boxplots of immune scores stratified by GEB in TCGA

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### IMMUNE ########################################################################################
format_immune_cell_plot_data <- function(dtf, subtype, gene = NULL, feature, filename) {
	if (subtype == 'HER2') {
		dtf$bds <- sign(dtf$ERBB2)
		dtf_st <- dtf[dtf$pam50 == 'Her2',]
	} else {
		dtf_st <- dtf[dtf$ic10 %in% c(1,2,6,9),]
		dtf_st$bds <- sign(colSums(rbind(
			rowSums(sign(dtf_st[,c('MYC','SQLE','FBXO32')]))*dtf_st$MYC_CNA,
			rowSums(sign(dtf_st[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*dtf_st$RPS6KB1_CNA,
			rowSums(sign(dtf_st[,c('CCND1','RSF1','PAK1',"NARS2")]))*dtf_st$RSF1_CNA,
			rowSums(sign(dtf_st[,c('ZNF703','FGFR1','LETM2','EIF4EBP1')]))*dtf_st$ZNF703_CNA
			)))
	}

	plot_data <- dtf_st[,c(feature,'bds')]
	counts <- table(plot_data$bds)
	if (counts['-1'] > counts['1']) {
		plot_data$burden <- (plot_data$bds >= 0)*1
	} else {
		plot_data$burden <- (plot_data$bds > 0)*1
	}
	plot_data$subtype <- subtype
	return(plot_data)
	}

create_immune_boxplot <- function(dtf, feature, ylimits, yat, ylab.label, text.y) {
	plot_data_cell_boxplot <- rbind(
		format_immune_cell_plot_data(tcga, subtype = 'HER2', feature = feature),
		format_immune_cell_plot_data(tcga, subtype = 'ER',feature = feature)
		)

	plot_data_cell_boxplot$bds_code <- NA
	plot_data_cell_boxplot[which(plot_data_cell_boxplot$burden == '1'),'bds_code'] <- 'B'
	plot_data_cell_boxplot[which(plot_data_cell_boxplot$burden == '0'),'bds_code'] <- 'A'
	plot_data_cell_boxplot <- plot_data_cell_boxplot[!is.na(plot_data_cell_boxplot$bds_code),]
	plot_data_cell_boxplot$group <- paste(
		plot_data_cell_boxplot$subtype,
		plot_data_cell_boxplot$bds_code,
		sep = '_'
		)

	p <- list()
	fc <- list()
	for (i in c('HER2','ER')) {
		tmp <- plot_data_cell_boxplot[plot_data_cell_boxplot$subtype == i,]
		stats <- wilcox.test(
			tmp[tmp$bds_code == 'B', feature],
			tmp[tmp$bds_code == 'A', feature]
			)
		fc[[i]] <- median(tmp[tmp$bds_code == 'B',feature], na.rm = TRUE)-median(tmp[tmp$bds_code == 'A',feature], na.rm = TRUE)
		p[[i]] <- stats$p.value
	}

	create.boxplot(
		as.formula(paste(feature, '~ group')),
		data = plot_data_cell_boxplot,
		add.stripplot = TRUE,
		ylab.label = ylab.label,
		filename = filename,
		xaxis.lab = rep(c('Low','High'), 3),
		ylimits = ylimits,
		yat = yat,
		xlab.label = 'Epitope Burden',
		add.rectangle = TRUE,
		xleft.rectangle = 0,
		ybottom.rectangle = ylimits[1]-10,
		xright.rectangle = 2.5,
		ytop.rectangle = ylimits[2]+10,
		col.rectangle = 'grey50',
		alpha.rectangle = 0.5,
		add.text = TRUE,
		text.labels = c('ER+','HER2+',
			paste('ES=', round(fc[['ER']], digits = 2), '\nP=', round(p[['ER']], digits = 2)),
			paste('ES=', round(fc[['HER2']], digits = 2), '\nP=', round(p[['HER2']], digits = 2))
			),
		text.x = c(1.5,3.5,1.5,3.5),
		text.y = text.y,
		text.anchor = 'center',
		text.col = 'black',
		text.cex = 1.2,
		text.fontface = c(rep('bold',2), rep('plain',2)),
		resolution = 300
		)
}

### MAIN ##########################################################################################
# read in megatable 
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
	as.is = TRUE
	)

# create lymphocytes plot
create_immune_boxplot(dtf = tcga, 
	feature = 'Lymphocytes', ylimits = c(0,1), yat = seq(0,1,0.2), 
	ylab.label = 'Lymphocytes', text.y = c(rep(0.98,2), rep(0.9, 2)),
	filename = paste0(date, '_supplementary_figure4a.png'))

# create CD8+ T cells plot
create_immune_boxplot(dtf = tcga, 
	feature = 'T.cells.CD8', ylimits = c(-0.03,0.4), yat = seq(0,0.4,0.1), 
	ylab.label = 'CD8+ T Cells', text.y = c(rep(0.39,2), rep(0.35, 2)),
	filename = paste0(date, '_supplementary_figure4b.png'))

# create cytotoxic score plot
create_immune_boxplot(dtf = tcga, 
	feature = 'cytotoxic_score', ylimits = c(0,1000), yat = seq(0,1000,500), 
	ylab.label = 'Cytotoxic Score', text.y = c(rep(970,2), rep(900, 2)),
	filename = paste0(date, '_supplementary_figure4c.png'))

# create macrophages plot
create_immune_boxplot(dtf = tcga, 
	feature = 'Macrophages', ylimits = c(0,1), yat = seq(0,1,0.2), 
	ylab.label = 'Macrophages', text.y = c(rep(0.98,2), rep(0.9, 2)),
	filename = paste0(date, '_supplementary_figure4d.png'))

# create macrophages M2 plot
create_immune_boxplot(dtf = tcga, 
	feature = 'Macrophages.M2', ylimits = c(0,0.8), yat = seq(0,1,0.2), 
	ylab.label = 'Macrophages M2-like', text.y = c(rep(0.78,2), rep(0.71, 2)),
	filename = paste0(date, '_supplementary_figure4e.png'))

# create macrophages M1 plot
create_immune_boxplot(dtf = tcga, 
	feature = 'Macrophages.M1', ylimits = c(-0.03,0.24), yat = seq(0,0.3,0.1), 
	ylab.label = 'Macrophages M1-like', text.y = c(rep(0.23,2), rep(0.21, 2)),
	filename = paste0(date, '_supplementary_figure4f.png'))

# create mhc class I presentation
create_immune_boxplot(dtf = tcga, 
	feature = 'MHC1_21978456', ylimits = c(-2.5,1.7), yat = seq(-2,2,1), 
	ylab.label = 'MHC Class I Presentation', text.y = c(rep(1.6,2), rep(1.3, 2)),
	filename = paste0(date, '_supplementary_figure4m.png'))


# create mhc class I presentation
create_immune_boxplot(dtf = tcga, 
	feature = 'APM1', ylimits = c(0.3,0.57), yat = seq(0.3,0.5,0.1), 
	ylab.label = 'MHC Class I Presentation', text.y = c(rep(0.56,2), rep(0.54, 2)),
	filename = paste0(date, '_supplementary_figure4n.png'))


