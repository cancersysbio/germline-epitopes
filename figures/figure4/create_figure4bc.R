### FIGURE 4BC #####################################################################################
# create heatmap of immune features

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(ConsensusClusterPlus)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()

### REFORMAT BINDERS IMMUNE AND RNA DATAFRAME #####################################################
reformat_bds_immune_rna_df <- function(dtf, subtype) {
	# subset down to only her2+
	if (subtype == 'HER2') {
		dtf_st <- dtf[dtf$pam50 == 'Her2',]
		dtf_st$bds <- sign(dtf_st$ERBB2)
	} else if (subtype == 'ER') {
		dtf_st <- dtf[which(dtf$ic10 %in% c(1,2,6,9)),]
		dtf_st$bds <- sign(colSums(rbind(
				rowSums(sign(dtf_st[,c('MYC','SQLE','FBXO32')]))*dtf_st$MYC_CNA,
				rowSums(sign(dtf_st[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*dtf_st$RPS6KB1_CNA,
				rowSums(sign(dtf_st[,c('CCND1','RSF1','PAK1',"NARS2")]))*dtf_st$RSF1_CNA,
				rowSums(sign(dtf_st[,c('ZNF703','FGFR1','LETM2','EIF4EBP1')]))*dtf_st$ZNF703_CNA
				)))
		}
	return(dtf_st)
}

### MAIN ##########################################################################################
# read in megatable
dtf <- read.delim(
  file.path(main_repo_path,'data','cohort_megatables','tcga_megatable.txt'),
	as.is = TRUE
	)

# set features to test
cells_cibersort <- c('Lymphocytes','T.cells.CD8','Macrophages','Macrophages.M2')
cells_gene <- c('PTPRC','FOXP3','CD68')
cytokines_sig <- c('Module3_IFN_score','Chemokine12_score','IL12_score_21050467',
	'STAT1_score','TGFB_score_21050467','IL8_21978456')
cytokines_gene <- c('GZMB','PRF1','IFNG','IL7','IL15','IL10','TGFB1')
stroma <- c('FAP','COL1A1','FN1')

### CREATE HER2 PLOT ##############################################################################
# reformat binders, immune and rna features 
bds <- reformat_bds_immune_rna_df(
	dtf = dtf,
	subtype = 'HER2'
	)
# reformat plot data
plot_data <- bds[,c(cells_cibersort, cells_gene, cytokines_sig, cytokines_gene, stroma)]
rownames(plot_data) <- bds[,'sample']

# ensure that input data matrix is equal to what the heatmap clustering produces
for (i in 1:ncol(plot_data)) {
	plot_data[,i] <- scale(as.numeric(as.character(plot_data[,i])))
}

# confirm clustering with consensus clustering
results <- ConsensusClusterPlus(t(plot_data),maxK=3,reps=50,pItem=0.8,pFeature=1,
	title='immune_clustering',clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
c2 <- data.frame(
	sample = names(results[[2]]$consensusClass),
	cluster = results[[2]]$consensusClass
	)
c2 <- merge(c2, bds[,c('sample','bds')], by = 'sample')
c2$bds <- (c2$bds >= 0)*1
fisher.test(table(c2[,c('bds','cluster')]))
c2 <- c2[order(c2$cluster, -c2$bds),]
plot_data <- plot_data[c2$sample,]

# set covariates
cov <- c2$bds
cov[cov == 1] <- 'lightpink3'
cov[cov == 0] <- 'lavenderblush'

# Dendrogram provided
dendrogram <- BoutrosLab.plotting.general::create.dendrogram(
         x = plot_data,
         cluster.dimension = 'row'
         );
# set colours for clusters 
clust1_cov <- rep(NA, nrow(plot_data))
clust1_cov[rownames(plot_data) %in% labels(dendrogram)[1:26]] <- 'plum4'
clust1_cov[rownames(plot_data) %in% labels(dendrogram)[27:nrow(plot_data)]] <- 'plum3'

# set colours for clusters from consensus clustering
clust2_cov <- rep(NA, nrow(plot_data))
clust2_cov[rownames(plot_data) %in% c2[c2$cluster == 1,'sample']] <- 'mediumpurple4'
clust2_cov[rownames(plot_data) %in% c2[c2$cluster == 2,'sample']] <- 'mediumpurple2'

covariate <- list(
         rect = list(
             col = 'transparent',
             fill = clust1_cov,
             lwd = 1.5
             ),
         rect = list(
             col = 'transparent',
             fill = clust2_cov,
             lwd = 1.5
             ),
         rect = list(
             col = 'transparent',
             fill = cov,
             lwd = 1.5
             )
         )

cluster1 <- c('TGFB1','FN1','COL1A1','FAP','TGFB_score_21050467',
	'Macrophages.M2','Macrophages')
cluster2 <-c('FOXP3','IFNG','PRF1','GZMB','STAT1_score','IL12_score_21050467',
	'Chemokine12_score','IL10','PTPRC','T.cells.CD8','Lymphocytes','IL7','IL15','IL8_21978456','Module3_IFN_score','CD68')

immune_cov <- colnames(plot_data)
immune_cov[immune_cov %in% cluster1] <- '#BD6656'
immune_cov[immune_cov %in% cluster2] <- '#387D9A'

immune_covariate <- list(
         rect = list(
             col = 'black',
             fill = immune_cov,
             lwd = 1.5
             )
         )

cov.legend <- list(
         legend = list(
             colours = c('lightpink3', 'lavenderblush'),
             labels = c('High','Low'),
             title = 'Epitope Burden'
             ),
         legend = list(
             colours = c('plum4', 'plum3'),
             labels = c('Cluster1','Cluster2'),
             title = 'Cluster (Diana)'
             ),
         legend = list(
             colours = c('mediumpurple4', 'mediumpurple2'),
             labels = c('Cluster1','Cluster2'),
             title = 'Cluster (ConsensusClustering)'
             ),
         legend = list(
             colours = c('#BD6656','#387D9A'),
             labels = c('Myleoid Predominant','Lymphocyte Predominant'),
             title = 'Immune Landscape'
             )
         );

create.heatmap(
	plot_data,
	yaxis.lab = paste0('   ', gsub('.', ' ', gsub('_score|_21050467|_21978456','', colnames(plot_data)), fixed = TRUE)),
	yaxis.cex = 1.2,
	cluster.dimensions = 'row',
	filename = paste0(date, '_figure4b.png'),
	covariates.top = covariate,
	covariates = immune_covariate,
	ylab.label = 'Immune Features',
	at = seq(-2,2,0.01),
	ylab.cex = 1.5,
	xlab.cex = 1.5,
	xlab.label = 'Samples',
	colour.scheme = c('dodgerblue2','white','firebrick2'),
	width = 9,
	height = 7,
	covariate.legend = cov.legend,
	legend.side = 'right',
	legend.cex = 1.2,
	colourkey.cex = 1.5,
	resolution = 300
	)

### CREATE ER PLOT ##############################################################################
# reformat binders, immune and rna features 
bds <- reformat_bds_immune_rna_df(
	dtf = dtf,
	subtype = 'ER'
	)
# reformat plot data
plot_data <- bds[,c(cells_cibersort,cells_gene, cytokines_sig, cytokines_gene, stroma)]
rownames(plot_data) <- bds[,'sample']
# remove samples without immune scores
plot_data <- plot_data[-unique(which(is.na(plot_data), arr.ind = TRUE)[,'row']),]

# ensure that input data matrix is equal to what the heatmap clustering produces
for (i in 1:ncol(plot_data)) {
	plot_data[,i] <- scale(as.numeric(as.character(plot_data[,i])))
}

# confirm clustering with consensus clustering
results <- ConsensusClusterPlus(t(plot_data),maxK=3,reps=50,pItem=0.8,pFeature=1,
	title='er_cna_immune_clustering',clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
c2 <- data.frame(
	sample = names(results[[2]]$consensusClass),
	cluster = results[[2]]$consensusClass
	)
c2 <- merge(c2, bds[,c('sample','bds')], by = 'sample')
c2$bds <- (c2$bds > 0)*1
fisher.test(table(c2[,c('bds','cluster')]))
c2 <- c2[order(c2$cluster, -c2$bds),]
plot_data <- plot_data[c2$sample,]

# set covariates
cov <- c2$bds
cov[cov == 1] <- 'lightpink3'
cov[cov == 0] <- 'lavenderblush'

# Dendrogram provided
dendrogram <- BoutrosLab.plotting.general::create.dendrogram(
         x = plot_data,
         cluster.dimension = 'row'
         );
# set colours for clusters 
clust1_cov <- rep(NA, nrow(plot_data))
clust1_cov[rownames(plot_data) %in% labels(dendrogram)[1:76]] <- 'plum4'
clust1_cov[rownames(plot_data) %in% labels(dendrogram)[77:nrow(plot_data)]] <- 'plum3'

# set colours for clusters from consensus clustering
clust2_cov <- rep(NA, nrow(plot_data))
clust2_cov[rownames(plot_data) %in% c2[c2$cluster == 1,'sample']] <- 'mediumpurple4'
clust2_cov[rownames(plot_data) %in% c2[c2$cluster == 2,'sample']] <- 'mediumpurple2'

covariate <- list(
         rect = list(
             col = 'transparent',
             fill = clust1_cov,
             lwd = 1.5
             ),
          rect = list(
             col = 'transparent',
             fill = clust2_cov,
             lwd = 1.5
             ),
         rect = list(
             col = 'transparent',
             fill = cov,
             lwd = 1.5
             )
         )

cluster1 <- c('TGFB1','FN1','COL1A1','FAP','TGFB_score_21050467',
	'Macrophages.M2','Macrophages')
cluster2 <-c('FOXP3','IFNG','PRF1','GZMB','STAT1_score','IL12_score_21050467',
	'Chemokine12_score','IL10','PTPRC','T.cells.CD8','Lymphocytes','IL7','IL15','IL8_21978456','Module3_IFN_score','CD68')

immune_cov <- colnames(plot_data)
immune_cov[immune_cov %in% cluster1] <- '#BD6656'
immune_cov[immune_cov %in% cluster2] <- '#387D9A'

immune_covariate <- list(
         rect = list(
             col = 'black',
             fill = immune_cov,
             lwd = 1.5
             )
         )

cov.legend <- list(
         legend = list(
             colours = c('lightpink3', 'lavenderblush'),
             labels = c('High','Low'),
             title = 'Epitope Burden'
             ),
         legend = list(
             colours = c('plum4', 'plum3'),
             labels = c('Cluster1','Cluster2'),
             title = 'Cluster (Diana)'
             ),
         legend = list(
             colours = c('mediumpurple4', 'mediumpurple2'),
             labels = c('Cluster1','Cluster2'),
             title = 'Cluster (ConsensusClustering)'
             ),
         legend = list(
             colours = c('#BD6656','#387D9A'),
             labels = c('Myeloid Predominant','Lymphocyte Predominant'),
             title = 'Immune Landscape'
             )
         );

create.heatmap(
	plot_data,
	yaxis.lab = paste0('   ', gsub('.', ' ', gsub('_score|_21050467|_21978456','', colnames(plot_data)), fixed = TRUE)),
	yaxis.cex = 1.2,
	cluster.dimensions = 'row',
	filename = paste0(date, '_figure4c.png'),
	covariates.top = covariate,
	covariates = immune_covariate,
	ylab.label = 'Immune Features',
	at = seq(-2,2,0.01),
	ylab.cex = 1.5,
	xlab.cex = 1.5,
	xlab.label = 'Samples',
	colour.scheme = c('dodgerblue2','white','firebrick2'),
	width = 9,
	height = 7,
	covariate.legend = cov.legend,
	legend.side = 'right',
	legend.cex = 1.2,
	colourkey.cex = 1.5,
	resolution = 300
	)


