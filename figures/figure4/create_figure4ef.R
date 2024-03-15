### create_figure4ef.R ############################################################################
# create dcis responders figure 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in summary data
dcis <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','dcis_megatable.txt'), 
	as.is = TRUE
	)

tbcrc_prog <- dcis[dcis$Cohort == 'TBCRC',]
tbcrc_prog$bds <- rowSums(sign(tbcrc_prog[,2:17]))
tbcrc_prog$burden <- (tbcrc_prog$bds > median(tbcrc_prog$bds))*1
# define progressors as those that have IBC recurrence
tbcrc_prog$progressors <- (tbcrc_prog$Diagnostic_Group == 'DCIS_with_IBC_recurrence')*1

# remove samples that recur as DCIS
plot_data <- as.data.frame(table(tbcrc_prog[tbcrc_prog$Diagnostic_Group != 'DCIS_with_DCIS_recurrence',c('burden','progressors')]))
plot_data$Prop <- plot_data$Freq/c(rep(sum(plot_data$Freq[1:2]),2), rep(sum(plot_data$Freq[3:4]),2))

# fit logistic regression
tmp_data <- tbcrc_prog[tbcrc_prog$Diagnostic_Group != 'DCIS_with_DCIS_recurrence',]
fit <- glm(progressors ~ burden + EA_PC1 + EA_PC2 + AA_PC1 + NA_PC1 + NA_PC2 + NA_PC3 + CO_PC1 + CO_PC2 + Her2_RNA + ER_RNA + grade,
	data = tmp_data,
	family = 'binomial')

create.barplot(
	Prop ~ progressors,
	groups = burden,
	data = plot_data,
	yat = seq(0,1,0.2),
	ylimits = c(0,1),
	ylab.label = 'Proportion',
	stack = TRUE,
	add.text = TRUE,
	text.labels = c('OR=0.37','P=0.035'),
	text.x = c(1,1),
	text.y = c(0.1,0.05),
	text.cex = 1.2,
	xaxis.lab = c('No IBC\nRecurrence', 'IBC\nRecurrence'),
	xlab.label = '',
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('lavenderblush','lightpink3')
                             ),
                         text = list(
                             lab = c('Low GEB','High GEB')
                             ),
                         padding.text = 5,
                         cex = 1.2
                         )
                     ),
                     # Positioning legend on plot
                     x = 1.01,
                     y = 0.6
                 )
             ),
	width = 7.5,
	right.padding = 15,
	col = c('lavenderblush','lightpink3'),
	filename = paste0(date, '_TBCRC_epitope_burden_barplot.pdf'),
	resolution = 300
	)

### FIGURE 4F #####################################################################################
# read in mibi ecad markers 
mibi <- read.delim(
  file.path(main_repo_path,'data','auxiliary_data','MIBI_ECAD_data.txt'),
	as.is = TRUE
	)

# add mibi markers
dcis <- merge(dcis, mibi, by = 'sample')
# calculate burden across all genes
dcis$bds <- rowSums(sign(dcis[,2:17]))
dcis$bds <- dcis$bds > median(dcis$bds)

# create plot
es <- round(
	median(dcis[dcis$bds == TRUE,'Myoep_cluster_fracECAD'])/median(dcis[dcis$bds == FALSE,'Myoep_cluster_fracECAD']),
	digits = 2
	)
p <- round(
	wilcox.test(dcis[dcis$bds == TRUE,'Myoep_cluster_fracECAD'], dcis[dcis$bds == FALSE,'Myoep_cluster_fracECAD'])$p.value,
	digits = 2
	)

dcis$bds <- factor(dcis$bds)

create.boxplot(
		Myoep_cluster_fracECAD ~ bds,
		add.stripplot = TRUE,
		xaxis.lab = c('Low','High'),
		xlab.label = 'GEB',
		ylab.label = '%ECAD in myoepithelium',
		data = dcis,
		filename = paste0(date, '_DCIS_bds_Myoep_cluster_fracECAD_boxplot.pdf'),
		legend = list(
			inside = list(
		                 fun = draw.key,
		                 args = list(
		                     key = list(
		                         text = list(
		                             lab = c(paste0('FC=', es), paste0('P=',p))  
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



