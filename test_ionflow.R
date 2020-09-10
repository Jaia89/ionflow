#' wl-03-07-2020, Fri: Load packages here.
#' wl-03-07-2020, Fri: qgraph loads plenty of R packages
#' wl-06-07-2020, Mon: debug functions pre_processing and exploratory_analysis
#' wl-30-07-2020, Thu: test on subset of IonData
#' wl-02-09-2020, Wed: test preProcessing with new data set

## ==== General settings ====
rm(list = ls(all = T))
#' set.seed(123)

#' tool_dir <- "C:/R_lwc/r_data/icl/"
tool_dir <- "~/my_galaxy/ionflow/"
setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2", 
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====

#' data for annotations
lib_dir <- paste0(tool_dir, "libraries/")
data_GOslim <- read.table("./libraries/data_GOslim.tsv", sep = "\t", header = T)
data_ORF2KEGG <- read.table("./libraries/data_ORF2KEGG.tsv", sep = "\t", header = T)

#' Load data set
#' ion_data <- read.table("./test-data/ionome_oe_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/ionome_ko_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/ionome_ko.tsv", header = T, sep = "\t")
ion_data <- read.table("./test-data/iondata_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/iondata.tsv", header = T, sep = "\t")

#' std_data <- read.table("./test-data/user_std.tsv", header = T, sep = "\t")
std_data <- NULL

## ==== Pre-processing ====
## names(ion_data)
##  [1] "orf_id_oe"   "plate_id_oe" "trays_id_oe" "batch_id_oe" "As"
##  [6] "Ca"          "Cd"          "Cl"          "Co"          "Cu"
## [11] "Fe"          "K"           "Mg"          "Mn"          "Mo"
## [16] "Na"          "Ni"          "P"           "S"           "Se"
## [21] "Zn"
pre_proc <- pre_processing(data = ion_data, stdev = std_data,
                           var_id = 1, batch_id = 2, data_id = 3)
## wl-02-09-2020, Wed: can select batch_id as 3 or 4

pre_proc$stats_raw_data
pre_proc$stats_outliers
pre_proc$stats_batch_data
pre_proc$stats_stand_data
head(pre_proc$data_long_bat)
head(pre_proc$data_long)
head(pre_proc$data_wide)
head(pre_proc$data_wide_symb)
pre_proc$plot_dot
pre_proc$plot_hist

## ==== Exploratory analysis ====
exp_anal <- exploratory_analysis(data = pre_proc$data_wide)

exp_anal$plot_pearson_correlation
exp_anal$plot_pca_individual
exp_anal$plot_heatmap
exp_anal$plot_pairwise_correlation_map
exp_anal$plot_correlation_network
head(exp_anal$data_pca_loadings)

## ==== Gene Clustering ====
gene_clust <- gene_clustering(data = pre_proc$data_wide,
                              data_symb = pre_proc$data_wide_symb,
                              thres_clus = 10, thres_anno = 5)
gene_clust$plot_profiles
gene_clust$stats_clusters
gene_clust$stats_kegg_goslim_annotation
gene_clust$stats_goterms_enrichment

## ==== Gene Network ====
gene_net <- gene_network(data = pre_proc$data_wide,
                         data_symb = pre_proc$data_wide_symb,
                         thres_clus = 5, thres_cor = 0.75)

gene_net$plot_pnet
gene_net$plot_impact_betweenness
gene_net$stats_impact_betweenness
gene_net$stats_impact_betweenness_tab
gene_net$stats_impact_betweenness_clus

