rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(beepr)

#GENERAL PARAMETERS
#AM_data_folder<-"./2019/AquaMaps"  #this models not available at 0.1 resolution
ME_data_folder<-"./2019/MaxEnt"
ANN_data_folder<-"./2019/ArtificialNeuralNetworks"
SVM_data_folder<-"./2019/SupportVectorMachines"
SVM_data_folderTS <-"./2019/SupportVectorMachines"

presence_data_folder<-"./Presence Records/vip"
csv_folder <- "./output/output_not_norm_historical/vip"
output_folder<-"./output/output_not_norm_historical_CSV/vip"





allspp<-list.files(presence_data_folder)
allspp<-gsub("Presence_","",allspp)
allspp<-gsub(".csv","",allspp)
#allspp<-gsub("_"," ",allspp)

for (species in allspp){
  
  cat("\nProcessing species ",species,"\n")
  #threshold SVM
  BMPFileSV<- paste0(SVM_data_folderTS,"/",species,"_metadata.txt")
  fileBMPSV<-readChar(BMPFileSV, file.info(BMPFileSV)$size)
  SVMTS<-gsub(".*Optimal decision threshold = ","",fileBMPSV)
  SVMTS<- as.double(gsub("\r\nOptimal kernel.*","",SVMTS))
  
  #threshold ANN
  BMPFileANN<- paste0(ANN_data_folder,"/",species,"_metadata.txt")
  fileBMPANN<-readChar(BMPFileANN, file.info(BMPFileANN)$size)
  ANNTS<-gsub(".*Optimal decision threshold = ","",fileBMPANN)
  ANNTS<- as.double(gsub("\nOptimal number of neurons.*","",ANNTS))
  
  #threshold ME
  BMPFileME<- paste0(ME_data_folder,"/",species,".html")
  fileBMPME<-readChar(BMPFileME, file.info(BMPFileME)$size)
  METS<-gsub(".*>Training omission rate</th><tr align=center><td>","",fileBMPME)
  METS<- gsub("</td><td>Fixed cumulative value 1.*","",METS)    
  METS<- as.double(gsub(".*</td><td>","",METS))
  
  #threshold AM  this models not available at 0.1 resolution
  AMTS<- 0.6
  
  data <- read.csv(paste0(csv_folder,"/ensemble_",species,".csv"))
  
  #ME column creation
  Thresholdfun <- function(val) {
    if(!is.na(val)){
    if (val < METS) return (as.double(0)) 
    else return (as.double(1))
    }else return(NA)
  }
  
  data$METS <- as.numeric(lapply(data$ME, function(val) sapply(val, Thresholdfun)))
  #ANN column creation
  Thresholdfun <- function(val) {
    if(!is.na(val)){
      if (val < ANNTS) return(as.double(0)) 
      else return(as.double(1))
    }else return(NA)
  }
  
  data$ANNTS <- as.numeric(lapply(data$ANN, function(val) sapply(val, Thresholdfun)))
  
  #SVM column creation
  Thresholdfun <- function(val) {
    if(!is.na(val)){
      if (val < SVMTS) return(as.double(0))
      else return(as.double(1))
    }else return(NA)
  }
  
  data$SVMTS <- as.numeric(lapply(data$SVM, function(val) sapply(val, Thresholdfun)))
  
  #AM column creation  this models not available at 0.1 resolution
  # Thresholdfun <- function(val) {
  #   if(!is.na(val)){
  #     if (val < AMTS) return(as.double(0)) 
  #     else return (as.double(1))
  #   }else return(NA)
  # }
  # 
  # data$AMTS <- as.numeric(lapply(data$AM, function(val) sapply(val, Thresholdfun)))
  

  #Agreement 
  #data$sum <- rowSums(data[,c("METS", "ANNTS","SVMTS","AMTS")], na.rm=TRUE)
  #data$sum <- rowSums(data[,c("METS", "ANNTS","SVMTS")], na.rm=TRUE)
  data$sum <- rowSums(data[,c("METS", "ANNTS","SVMTS")])
  data2<-data[,c("x","y","sum")]
  colnames(data2)<-c("longitude","latitude","sum_ts")
  write.csv(data2,paste0(output_folder,"/",species,".csv"),row.names = FALSE)

}
beep(4)