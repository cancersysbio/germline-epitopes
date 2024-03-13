### CREATE SUPPLEMENTARY FIGURE 1B #################################################################
# test presence of HLAs that can present IISAVVGIL with HER2 subtype
# generate supplementary figure 1b
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

date <- Sys.Date()

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-i', '--idir', type = 'character', help = 'icgc directory');
parser$add_argument('-m', '--mdir', type = 'character', help = 'metabric directory');

args <- parser$parse_args();

### GET ICGC HLAs #################################################################################
get_icgc_hlas <- function(samples, basedir) {
	# read in hlas 
	hlas <- do.call(rbind, sapply(
		samples,
		function(id) {
			hlafile <- file.path(basedir, 'ERBB2', id, 'HLA.txt')
			if (file.exists(hlafile)) {
				tmp <- read.delim(
					hlafile,
					header = FALSE,
					as.is = TRUE
					)
				tmp <- unique(tmp$V1)
				counts <- table(factor(substr(tmp, 1, 1), levels = c('A','B','C')))
				data.frame(
					sample = id,
					hla = tmp
					)
				}
			},
		simplify = FALSE
		))
	return(hlas)
}

### GET METABRIC HLAs #############################################################################
get_metabric_hlas <- function(samples, basedir) {
	# read in hlas 
	metabric_hlas <- do.call(rbind, sapply(
		samples$V1,
		function(id) {
			tmp <- read.delim(
				file.path(basedir, 'ERBB2', id, 'HLA.txt'),
				header = FALSE,
				as.is = TRUE
				)
			tmp <- unique(tmp$V1)
			counts <- table(factor(substr(tmp, 1, 1), levels = c('A','B','C')))
			data.frame(
				sample = id,
				hla = tmp
				)
			},
		simplify = FALSE
		))
	# only keep alleles with 80% imputation accuracy
	alleles <- read.delim(
	        'tcga_cookhla_hla_alleles_imputation_accuracy.txt',
	        as.is = TRUE
	        )
	alleles <- alleles[alleles$accuracy > 0.8,]
	metabric_hlas <- metabric_hlas[metabric_hlas$hla %in% alleles$hla,]
	return(metabric_hlas)
}


### COUNT NUMBER OF BINDING ALLELES ###############################################################
count_number_binding_alleles <- function(dtf, alleles) {
	# bin patients by alleles that bind 
	hlasbinddf <- do.call(rbind, sapply(
		unique(dtf$sample),
		function(x) {
			tmp <- dtf[dtf$sample == x,]
			data.frame(
				sample = x,
				hla = sum(tmp$hla %in% alleles)
				)
			},
		simplify = FALSE
		))
	return(hlasbinddf)
}

### APPLY LOGISTIC MODEL ##########################################################################
apply_logistic_model <- function(hlas) {
	# apply logistic regression model and return results
	fit <- glm(subtype ~ hla + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = hlas, family = 'binomial')
	ci <- confint(fit)
	res <- data.frame(
		coef = coef(fit)[['hla']],
		p = summary(fit)$coefficients['hla',4],
		se = summary(fit)$coefficients['hla',2],
		l95 = ci['hla',1],
		u95 = ci['hla',2]
		)
	return(res)
	}

### MAIN ##########################################################################################
# read in gp2 and e75 alleles
gp2_e75_hlas <- read.delim(
	'gp2_e75_hlas.txt',
	as.is = TRUE,
	header = FALSE
	)
gp2_e75_hlas <- gp2_e75_hlas$V1

### TEST IN ICGC ##################################################################################
# read in icgc megatable 
icgc <- read.delim(
	'icgc_megatable.txt',
	as.is = TRUE
	)
# read in all samples that have epitope predictions
icgc_hlas <- get_icgc_hlas(samples = icgc$sample_name, basedir = args$idir)

# count number of binding alleles 
icgc_hlas_both <- count_number_binding_alleles(icgc_hlas, gp2_e75_hlas)

# annotate 
colnames(icgc_hlas_both) <- c('Experiment.System.ID','hla')
icgc_hlas_both <- merge(icgc_hlas_both, icgc, by = 'Experiment.System.ID')
# assign subtype
icgc_hlas_both$subtype <- (icgc_hlas_both$final.HER2 == 'positive')*1
icgc_hlas_both$hla <- (icgc_hlas_both$hla >= median(icgc_hlas_both$hla))*1
icgc_res <- apply_logistic_model(icgc_hlas_both)
icgc_res$num <- nrow(icgc_hlas_both)
icgc_res$cohort <- 'ICGC'

### TEST METABRIC #################################################################################
# read in metabric 
metabric <- read.delim(
	'metabric_megatable.txt',
	as.is = TRUE
	)
# read in all samples that have epitope predictions
metabric_hlas <- get_metabric_hlas(samples = metabric$sample, basedir = args$mdir)

# count number of binding alleles 
metabric_hlas_both <- count_number_binding_alleles(metabric_hlas, gp2_e75_hlas)

# annotate 
metabric_hlas_both  <- merge(metabric_hlas_both, metabric, by = 'sample')
# assign subtype
metabric_hlas_both$subtype <- (metabric_hlas_both$CLAUDIN_SUBTYPE == 'Her2')*1
metabric_hlas_both$hla <- (metabric_hlas_both$hla > median(metabric_hlas_both$hla))*1
metabric_res <- apply_logistic_model(metabric_hlas_both)
metabric_res$num <- nrow(metabric_hlas_both)
metabric_res$cohort <- 'METABRIC'

### PLOT ##########################################################################################
# create plot data
plot_data <- rbind(
	icgc_res,
	metabric_res
	)
plot_data$index <- 1:nrow(plot_data)

meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                        data = plot_data)
res <- rma(yi, vi, data=meta_data)


create.scatterplot(
        index ~ coef,
        data = plot_data,
        horizontal = TRUE,
        xlimits = c(-2,2),
        filename = paste0(date, '_e75_gp2_scatterplot.png'),
        xat = log(c(0.5, 1, 2.5)),
        xaxis.lab = c(0.5, 1, 2.5),
        xlab.label = 'Odds Ratio',
        ylab.label = 'Cohort',
        yaxis.lab = plot_data$cohort,
        yat = 1:nrow(plot_data),
        ylimits = c(0.5,2.5),
        add.text = TRUE,
        text.labels = c('OR=0.73','P=0.04', paste('n=', icgc_res$num), paste('n=', metabric_res$num)),
        text.y = c(rev(c(0.75,1)),1.25, 2.25),
        text.x = c(1.5,1.5, plot_data$coef),
        text.cex = 1.2,
        main.cex = 2,
        abline.v = 0,
        width = 10,
        height = 4,
        x.error.right = plot_data$u95-plot_data$coef,
        x.error.left = plot_data$coef-plot_data$l95,
        top.padding = 3,
        resolution = 300
        )





