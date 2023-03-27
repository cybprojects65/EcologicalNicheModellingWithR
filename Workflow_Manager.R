rm(list=ls(all=TRUE))
library(properties)

### ENM Workflow Manager ###
Workflow_Parameters_File <- paste0("./Workflow_Parameters.txt")
WPFprops <- read.properties(Workflow_Parameters_File)

EcologicalNicheModeler <- as.logical(WPFprops$"EcologicalNicheModeler")
if (EcologicalNicheModeler == TRUE){
  source("./EcologicalNicheModeler.R")
}

#MaxEnt_converter <- as.logical(WPFprops$"MaxEnt_converter")
#if (MaxEnt_converter == TRUE){
 # source("./MaxEnt_converter.R")
#}

#CSV_converter <- as.logical(WPFprops$"CSV_converter")
#if (CSV_converter == TRUE){
 # source("./CSV_converter.R")
#}

EnsambleModeler <- as.logical(WPFprops$"EnsambleModeler")
if (EnsambleModeler == TRUE){
  source("./EnsambleModeler.R")
}

AgreementIndex <- as.logical(WPFprops$"AgreementIndex")
if (AgreementIndex == TRUE){
  source("./AgreementIndex.R")
}

cat("\n All Done.\n")

