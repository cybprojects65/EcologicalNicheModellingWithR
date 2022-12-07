rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(beepr)

Ensemble_csv1 <- "./output/output_not_norm_historical_CSV/Fish"
Ensemble_csv2 <- "./output/output_not_norm_historical_CSV/Non Fish"
resolution <- 0.5
TS <- 3
output_folder <- "./output/output_biodiversity_index"

allspp1<-list.files(Ensemble_csv1)
allspp1<-gsub("ensembleTS_","",allspp1)
allspp1<-gsub(".csv","",allspp1)

allspp2<-list.files(Ensemble_csv2)
allspp2<-gsub("ensembleTS_","",allspp2)
allspp2<-gsub(".csv","",allspp2)



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



for (species in allspp1){
  cat("\nProcessing species ",species,"\n")
  datasp <- read.csv(paste0(Ensemble_csv1,"/ensembleTS_",species,".csv"))
  datasp$WTS <- as.numeric(lapply(datasp$sum_ts, function(val) sapply(val, Thresholdfun)))
  grid_of_points$index <- grid_of_points$index + datasp$WTS
  
}

for (species in allspp2){
  cat("\nProcessing species ",species,"\n")
  datasp <- read.csv(paste0(Ensemble_csv2,"/ensembleTS_",species,".csv"))
  datasp$WTS <- as.numeric(lapply(datasp$sum_ts, function(val) sapply(val, Thresholdfun)))
  grid_of_points$index <- grid_of_points$index + datasp$WTS
  
}
colnames(grid_of_points)<-c("longitude","latitude","index")
write.csv(grid_of_points,paste0(output_folder,"/Biodiversity_index_TS",TS,"historical.csv"),row.names = FALSE)

beep(4)