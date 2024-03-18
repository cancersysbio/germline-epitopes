### CREATE SUPPLEMENTARY FIGURE 1P ################################################################
# create supplementary figure 1p 
# plot showing the number of samples that fall in each subtype definition

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in tcga 
tcga <- read.delim(
	file.path(main_repo_path, 'data', 'cohort_megatables', 'tcga_megatable.txt'),
	as.is = TRUE
	)

# calculate number of samples in each subtype definition 
her2_pam50 <- sum(tcga$pam50 == 'Her2')
her2_cna <- sum(tcga$ERBB2_CNA == 1)
her2_ic5 <- sum(tcga$ic10 == 5)

ic1 <- sum(tcga$ic10 == 1)
ic1_cna <- sum(tcga$RPS6KB1_CNA == 1 & tcga$ER == 1)

ic2 <- sum(tcga$ic10 == 2)
ic2_cna <- sum(tcga$RSF1_CNA == 1 & tcga$ER == 1)

ic6 <- sum(tcga$ic10 == 6)
ic6_cna <- sum(tcga$ZNF703_CNA == 1 & tcga$ER == 1)

ic9 <- sum(tcga$ic10 == 9)
ic9_cna <- sum(tcga$MYC_CNA == 1 & tcga$ER == 1)

plot_data <- data.frame(
	subtype = c(rep('IC5',3), rep(c('IC1','IC2','IC6','IC9'), each = 2)),
	definition = c(c('PAM50','CNA','IC10'), rep(c('IC10','CNA'), 4)),
	count = c(her2_pam50, her2_cna, her2_ic5, ic1, ic1_cna, ic2, ic2_cna,
		ic6, ic6_cna, ic9, ic9_cna)
	)
plot_data <- plot_data[order(plot_data$subtype, -plot_data$count),]
plot_data$index <- 1:nrow(plot_data)

definition_col <- as.character(plot_data$definition)
definition_col[definition_col == 'CNA'] <- default.colours(5, palette.type = 'spiral.noon')[1]
definition_col[definition_col == 'IC10'] <- default.colours(5, palette.type = 'spiral.noon')[3]
definition_col[definition_col == 'PAM50'] <- default.colours(5, palette.type = 'spiral.noon')[5]

ic_cov <- as.character(plot_data$subtype)
ic_cov[ic_cov == 'IC6'] <- "gold"
ic_cov[ic_cov == 'IC9'] <- "#F8B4E3"
ic_cov[ic_cov == 'IC2'] <- "chartreuse4"
ic_cov[ic_cov == 'IC1'] <- "darkorange2"
ic_cov[ic_cov == 'IC5'] <- "#8B0000"

cov <- list(
        rect = list(
                col = 'transparent',
                fill = definition_col
                ),
        rect = list(
                col = 'transparent',
                fill = ic_cov
                )
        );

cov.grob <- covariates.grob(
        covariates = cov,
        ord = c(1:length(definition_col)),
        side = 'top',
        size = 1
        );

cov.legend <- list(
        legend = list(
                colours =  default.colours(5, palette.type = 'spiral.noon')[c(1,3,5)],
                labels = c('CNA','IC10','PAM50'),
                title = 'Subtype Definition',
                border = 'transparent'
                ),
        legend = list(
                colours =  unique(ic_cov),
                labels = c('IC1','IC2','IC5','IC6','IC9'),
                title = 'Subtype',
                border = 'transparent'
                )
        );

cov.legend.grob <- legend.grob(
        legends = cov.legend
        );

create.barplot(
	count ~ index,
	data = plot_data,
	filename = paste0(date, '_supplementary_figure1p.png'),
	xaxis.lab = rep('', nrow(plot_data)),
	ylab.label = 'Number of Samples',
	ylimits = c(0,160),
	xlab.label = '',
	yat = seq(0,150,50),
        add.rectangle = TRUE,
        xleft.rectangle = c(2.5, 7.5),
        ybottom.rectangle = 0,
        xright.rectangle = c(4.5, 9.5),
        ytop.rectangle = 160,
        col.rectangle = 'grey50',
        alpha.rectangle = 0.5,
	legend = list(
            bottom = list(fun = cov.grob),
            right = list(fun = cov.legend.grob)
            ),
	width = 9,
	height = 4,
	bottom.padding = 3,
	resolution = 300
	)