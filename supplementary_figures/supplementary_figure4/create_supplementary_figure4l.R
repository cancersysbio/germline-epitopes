### CREATE SUPPLEMENTARY FIGURE 4L ################################################################
# test metabric IMC

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in METABRIC IMC data
# calculated proportions of each cell type
# grouped the following as cancer cells: 'Basal CKlow','HER2+','HR- CK7-','HR- CK7+','HR- CKlow CK5+',
	#'HR- Ki67+','HR+ CK7-','HR+ CK7- Ki67+',
	#'HR+ CK7- Slug+','HRlow CKlow','Hypoxia'
imc_prop <- read.delim(
	file.path(main_repo_path, 'data', 'auxiliary_data', 'metabric_imc_proportions.txt'),
	as.is = TRUE
	)
# read in summary data
metabric <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'metabric_megatable.txt'),
	as.is = TRUE
	)

## HER2 ## 
meta_her2 <- metabric[which(metabric$CLAUDIN_SUBTYPE == 'Her2'),]
meta_her2$bds <- (sign(meta_her2$ERBB2) >= 0)*1

meta_her2 <- merge(
	meta_her2[,c('sample','bds','PC1','PC2','PC3','PC4','PC5','PC6','INTCLUST','stage','age','grade','ER_IHC')],
	imc_prop[imc_prop$celltype %in% c('cancer','T cells'),],
	by = 'sample'
	)

her2_res <- list()
for (i in c('cancer','T cells')) {
	tmp <- meta_her2[meta_her2$celltype == i,]
	tmp$prop_scale <- scale(tmp$proportion, center = FALSE)
	fit <- lm(prop_scale ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + INTCLUST + age + stage, data = tmp)
	ci <- confint(fit)
	her2_res[[i]] <- data.frame(
		celltype = i,
		coef = coef(fit)[['bds']],
		p = summary(fit)$coefficients['bds',4],
		l95 = ci['bds',1],
		u95 = ci['bds',2]
		)
}
her2_res <- do.call(rbind, her2_res)

## ER ##
meta_er <- metabric[which(metabric$INTCLUST %in% c(1,2,6,9) & metabric$stage > 1),]
meta_er$bds <- colSums(rbind(
	rowSums(sign(meta_er[,c('MYC','SQLE','FBXO32')]))*meta_er$MYC_CNA,
	rowSums(sign(meta_er[,c('RPS6KB1','TUBD1','DHX40','BCAS3')]))*meta_er$RPS6KB1_CNA,
	rowSums(sign(meta_er[,c('CCND1','RSF1','PAK1',"NARS2")]))*meta_er$RSF1_CNA,
	rowSums(sign(meta_er[,c('ZNF703','FGFR1','LETM2')]))*meta_er$ZNF703_CNA
	))
meta_er$bds <- (meta_er$bds > 0)*1

meta_er <- merge(
	meta_er[,c('sample','bds','PC1','PC2','PC3','PC4','PC5','PC6','INTCLUST','stage','age','grade')],
	imc_prop[imc_prop$celltype %in% c('cancer','T cells'),],
	by = 'sample'
	)
# Because one sample had ~50% of cells as T-cells which seems unlikely, filtering out outliers
erq3 <- quantile(meta_er[meta_er$celltype == 'T cells','proportion'], 0.75)
erq1 <- quantile(meta_er[meta_er$celltype == 'T cells','proportion'], 0.25)
thres <- erq3+3*(erq3-erq1)

i <- meta_er[meta_er$celltype == 'T cells' & meta_er$proportion > thres,'sample']
meta_er <- meta_er[!meta_er$sample %in% i,]

er_res <- list()
for (i in c('cancer','T cells')) {
	tmp <- meta_er[meta_er$celltype == i,]
	tmp$prop_scale <- scale(tmp$proportion)
	fit <- lm(prop_scale ~ bds + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + stage + age + INTCLUST, tmp)
	ci <- confint(fit)
	er_res[[i]] <- data.frame(
			celltype = i,
			coef = coef(fit)[['bds']],
			p = summary(fit)$coefficients['bds',4],
			l95 = ci['bds',1],
			u95 = ci['bds',2]
			)
	}
er_res <- do.call(rbind, er_res)

plot_data <- rbind(er_res[er_res$celltype == 'T cells',], her2_res[her2_res$celltype == 'T cells',])
plot_data$type <- c('ER','HER2')
plot_data$index <- 1:nrow(plot_data)

create.scatterplot(
        index ~ coef,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-3,3),
        filename = paste0(date, '_metabric_imc_scatterplot.png'),
        xlab.label = 'Coefficient',
        ylab.label = '',
        yaxis.lab = c('ER+','HER2+'),
        yat = 1:2,
        ylimits = c(0.5, nrow(plot_data)+0.75),
        main.cex = 2,
        abline.v = 0,
        x.error.right = plot_data$u95-plot_data$coef,
        x.error.left = plot_data$coef-plot_data$l95,
        width = 8,
        height = 2.5,
        top.padding = 2,
        resolution = 300
        )
