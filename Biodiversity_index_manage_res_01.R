rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(beepr)

Ensemble_csv <- "./output/output_not_norm_historical_CSV/vip"

resolution <- 0.1
TS <- 2
output_folder <- "./output/output_biodiversity_index"

allspp<-list.files(Ensemble_csv)
#allspp1<-gsub("ensembleTS_","",allspp1)
allspp<-gsub(".csv","",allspp)


xseq<-seq(from=-180,to=180,by=resolution)
yseq<-seq(from=-90,to=90,by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq)
grid_of_points$index <-0

Thresholdfun <- function(val) {
  if(!is.na(val)){
    if (val < TS) return(as.double(0)) 
    else return (as.double(1))
  }else return(NA)
}



for (species in allspp){
  cat("\nProcessing species ",species,"\n")
  datasp <- read.csv(paste0(Ensemble_csv,"/",species,".csv"))
  datasp$WTS <- as.numeric(lapply(datasp$sum_ts, function(val) sapply(val, Thresholdfun)))
  grid_of_points$index <- grid_of_points$index + datasp$WTS
  
}

colnames(grid_of_points)<-c("longitude","latitude","index")
write.csv(grid_of_points,paste0(output_folder,"/Biodiversity index TS ",TS," res 01.csv"),row.names = FALSE)

beep(4)