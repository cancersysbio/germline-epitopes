### CREATE SUPPLEMENTARY FIGURE 1T ################################################################
# create plot of benign vs pathogenic variants
# create supplementary figure 1T 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# set gene list
genes <- c('ERBB2','RPS6KB1','TUBD1','DHX40','BCAS3',
        'RSF1','CCND1','PAK1','NARS2','MYC','SQLE','FBXO32'
        )

# read in plot data 
plot_data <- read.delim(
        file.path(main_repo_path, 'data', 'controls', 'pathogenic_benign_variants.txt'),
        as.is = TRUE
        )

plot_data$index <- NA
for (i in 1:length(genes)) {
        plot_data[plot_data$gene == genes[i],'index'] <- letters[i]
}

create.barplot(
        Freq ~ index,
        groups = plot_data$pathogenic,
        data = plot_data,
        xaxis.lab = genes,
        xaxis.cex = 1,
        xaxis.rot = 45,
        width = 9,
        ylab.label = 'Number Unique Variants',
        xlab.label = 'Gene',
        ylimits = c(0,40),
        yat = seq(0,40,10),
        col = c('darkgrey','black'),
        add.rectangle = TRUE,
        xleft.rectangle = c(5.5,9.5),
        ybottom.rectangle = 0,
        xright.rectangle = c(1.5,13.5),
        ytop.rectangle = 40,
        col.rectangle = 'grey',
        alpha.rectangle = 0.5,
        add.text = TRUE,
        text.labels = c('HER2+','IC1','IC2','IC9'),
        text.x = c(1,3.5,7.5,11),
        text.y = 38,
        xaxis.fontface = 4,
        filename = paste0(date, '_supplementary_figure1t.png'),
        legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('darkgrey', 'black')
                             ),
                         text = list(
                             lab = c('Benign','Pathogenic')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     ),
                     # Positioning legend on plot
                     x = 0.8,
                     y = 0.85
                 )
             ),
        resolution = 300
        )




