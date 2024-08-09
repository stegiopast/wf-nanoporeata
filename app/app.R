#!/usr/bin/env Rscript

################################################################################
##                                                                            ##
##                                NANOPOREATA                                 ##
##                                                                            ##
################################################################################

# ______________________________________________________________________________
# LIBRARIES ####
#options(repos = list(CRAN="http://cran.rstudio.com/"))
#source("install.R", local = T)
for (lib in names(readRDS("NanopoReaTA_Rpackages.RDS"))) { library(lib, character.only = TRUE) }

# Save list of required R packages + version
# package_versions = data.table::rbindlist(lapply(sessionInfo()$otherPkgs, function(i) data.frame("Version" = i$Version)), idcol = "pkg")
# write.table(package_versions, "NanopoReaTA_Rpackage_versions.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# ______________________________________________________________________________
# SETTINGS ####
options(shiny.maxRequestSize = 30*1024^2)

# ______________________________________________________________________________
# FUNCTIONS ####
## DEA ####
source("/nanoporeata/app/server/R_scripts/dea_function.R", local = TRUE)

## DTE ####
source("/nanoporeata/app/server/R_scripts/dte_function.R", local = TRUE)

## DTU #### 
source("/nanoporeata/app/server/R_scripts/dtu_function.R", local = TRUE)

## READ LENGTH DISTRIBUTION ####
source("/nanoporeata/app/server/R_scripts/read_length_distribution_plots.R", local = TRUE)

## GENE WISE ANALYSIS ####
source("/nanoporeata/app/server/R_scripts/gene_wise_analysis_function.R", local  = TRUE)

## INFER EXPERIMENT ####
source("/nanoporeata/app/server/R_scripts/infer_experiment_plots.R", local = TRUE)

## Gene Body coverage plots ####
source("/nanoporeata/app/server/R_scripts/plot_single_gene.R", local = TRUE)
#____________________________________________________________________________
# FRONTEND
source("/nanoporeata/app/ui/ui.R", local = TRUE)
# ______________________________________________________________________________
# BACKEND
source("/nanoporeata/app/server/server.R", local = TRUE)


metadata_path = "/data/metadata.tsv"
bed_file = "/data/basic_bed_file.bed"
gtf_path = "/data/converted_gtf.csv"
while (!file.exists(metadata_path) | !file.exists(gtf_path)){
  print("Files not present....Please wait!")
  Sys.sleep(10)
}

# ______________________________________________________________________________
# LAUNCH APP
APP <- shinyApp(ui, server)#, options=list(port=as.integer(Sys.getenv("PORT"))))
runApp(APP,host = '0.0.0.0',port=as.integer(Sys.getenv("PORT")))