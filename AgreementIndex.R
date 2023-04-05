rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(properties)

Workflow_Parameters_File<- paste0("./Workflow_Parameters.txt")
WPFprops <- read.properties(Workflow_Parameters_File)

#input parameters
index <- WPFprops$"index"
cat(paste0("\n",index," calculation\n"))

#general parameters
Ensemble_folder <- "./Ensemble"
output_folder <- "."

model_projection <- as.logical(WPFprops$"model_projection")
if (model_projection == T){
  methods_folder<- gsub(pattern = "\\\"",replacement = "", x=WPFprops$output_folder_projection)
} else methods_folder<- gsub(pattern = "\\\"",replacement = "", x=WPFprops$output_folder)

all_IDs<-list.files(Ensemble_folder)
all_IDs<-gsub("ensemble_","",all_IDs)
all_IDs<-gsub("\\.asc","",all_IDs)

#set resolution from metadata
allmet<-list.dirs(methods_folder, recursive = FALSE)
method<- allmet[2]
BMPFile<- paste0(method,"/",all_IDs[1],"_metadata.txt")
BMPprops <- read.properties(BMPFile)
resolution<- as.double(BMPprops$"Spatial resolution")


#grid preparation
asc_sample <- raster(paste0(Ensemble_folder,"/ensemble_",all_IDs[1],".asc"))
min_x_in_raster<-asc_sample@extent[1]
max_x_in_raster<-asc_sample@extent[2]
min_y_in_raster<-asc_sample@extent[3]
max_y_in_raster<-asc_sample@extent[4]

xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq)

#set threshold and function
nativedistribution <- as.logical(WPFprops$"nativedistribution")
if (nativedistribution == T){
  TS<-length(allmet)-1
} else TS<-1

cat("\nAssessing presence by threshold = ",TS,"\n")
for (ID in all_IDs){
  cat("\nProcessing ",ID,"\n")
  dataID <- raster(paste0(Ensemble_folder,"/ensemble_",ID,".asc"))
  dataID <- as.data.frame(dataID, xy = TRUE)
  colnames(dataID) <- c('longitude','latitude','sum_ts')
  newcolumn<-sapply(1:length(dataID$sum_ts), function(idx){
    val<-dataID$sum_ts[idx]
    if(!is.na(val)){
      if (val ==-9999) return (-9999)
      if (val < TS) return(as.double(0)) 
      else return (as.double(1))
    }else return(-9999)
  },simplify = T)
  
  grid_of_points<- cbind(grid_of_points, newcolumn)
  
}

cat("\n",index," calculation\n")
dataIDs<-grid_of_points[,-c(1:2)]
sums<-sapply(1:dim(dataIDs)[1], function(idx){
  data_row<-dataIDs[idx,]
  naindex<-which(data_row == -9999)
  if (length(naindex) == length(data_row))
    return(-9999)
  else{
    data_row[naindex]<-0
    return (sum(data_row))  
  }
  
},simplify = T)

grid_of_points$index <- sums
grid_of_points<-grid_of_points[,c("x","y","index")]
colnames(grid_of_points)<-c("longitude","latitude","index")

#write the ASC file
#define the grid points to populate the dataset
ypoints<-unique(grid_of_points$latitude)
xpoints<-unique(grid_of_points$longitude)
ncol_r<-length(xpoints)
nrow_r<-length(ypoints)
#create a new raster with the same extent and resolution of the first layer
ro <- raster(ncol=ncol_r, nrow=nrow_r)
length(values(ro))

#populate the matrix
values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
row_counter<-0
for (y_c in 1:(nrow_r)){
  yp<-ypoints[y_c]
  row_rast<-grid_of_points[which(grid_of_points$latitude == yp),]
  row_rast<-row_rast[order(row_rast$longitude),]
  #values[(nrow_r-row_counter+1),]<-row_rast$index[1:(ncol_r)]
  values[(row_counter+1),]<-row_rast$index[1:(ncol_r)]
  row_counter<-row_counter+1
}
values_vec<-as.vector(t(values))
values(ro)<-values_vec
NAvalue(ro)<- -9999

#save the raster
writeRaster(ro, filename=paste0(output_folder,"/",index,"_ensembleThd_",TS,".csv"), format="ascii",overwrite=TRUE)