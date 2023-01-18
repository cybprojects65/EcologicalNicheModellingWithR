rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(properties)

#input parameters
method <- "ArtificialNeuralNetwork"
resolution <- 0.5
set_threshold <- 0.6
changing_thr <- T

#general parameters
main_folder<-"."
results_folder <- "./Results"
ENV_parameters_folder <- "./ENV"
dir.create(file.path(main_folder, "ENM"), showWarnings = FALSE)
ENM_folder <- "./ENM"
dir.create(file.path(ENM_folder, method), showWarnings = FALSE)

allenvpar <- list.files(ENV_parameters_folder)

all_IDs<-list.files(results_folder)
all_IDs<-gsub(".csv","",all_IDs)
all_IDs<-gsub(".txt","",all_IDs)
all_IDs<-unique(all_IDs)
ID <-"Alopias_volpinus"

#grid preparation
asc_sample <- raster(paste0(ENV_parameters_folder,"/",allenvpar[1]))
min_x_in_raster<-asc_sample@extent[1]
max_x_in_raster<-asc_sample@extent[2]
min_y_in_raster<-asc_sample@extent[3]
max_y_in_raster<-asc_sample@extent[4]

xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
grid_of_points<-expand.grid(x = xseq, y = yseq)

cat("\nGenerating .asc and metadata.txt files\n")
for(ID in all_IDs){
  cat("\nProcessing ",ID,"\n")

  #changing threshold
  if(changing_thr == T){
     BMPFile<- paste0(results_folder,"/",ID,".txt")
     BMPprops <- read.properties(BMPFile)
     threshold<- as.double(BMPprops$"Optimal decision threshold")
    } else threshold <- set_threshold

  #writing metadata.txt
  text_lines <- c(paste0("\nModel = ",method),
                  paste0("\nSpatial resolution = ",resolution),
                  paste0("\nOptimal decision threshold = ",threshold))
  writeLines(text_lines, paste0(ENM_folder,"/",method,"/",ID,"_metadata.txt"))
  
  datasp <- read.csv(paste0(results_folder,"/",ID,".csv"))
  datasp <- suppressMessages(full_join(grid_of_points, datasp))
  
  #write the ASC file
  #define the grid points to populate the dataset
  ypoints<-unique(grid_of_points$y)
  xpoints<-unique(grid_of_points$x)
  ncol_r<-length(xpoints)
  nrow_r<-length(ypoints)
  #create a new raster with the same extent and resolution of the first layer
  ro <- raster(ncol=ncol_r, nrow=nrow_r)
  length(values(ro))
  
  res(ro) <- resolution
  length(values(ro))
  extent(ro)<-extent(asc_sample)
  
  #populate the matrix
  values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
  row_counter<-1
  for (y_c in 1:(nrow_r)){
    yp<-ypoints[y_c]
    row_rast<-datasp[which(grid_of_points$y == yp),]
    row_rast<-row_rast[order(row_rast$x),]
    values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  values_vec<-as.vector(t(values))
  values(ro)<-values_vec
  NAvalue(ro)<- -9999
  #plot(ro)
  
  #save the raster
  writeRaster(ro, filename=paste0(ENM_folder,"/",method,"/",ID), format="ascii",overwrite=TRUE)
}
