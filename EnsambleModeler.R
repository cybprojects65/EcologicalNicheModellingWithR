rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(properties)

cat("\nEnsemble creation\n")

Workflow_Parameters_File<- paste0("./Workflow_Parameters.txt")
WPFprops <- read.properties(Workflow_Parameters_File)

#general parameters
main_output_folder<-"."
new_folder<- "Ensemble"
dir.create(file.path(main_output_folder, new_folder), showWarnings = FALSE)
output_folder<-"./Ensemble"
methods_folder<- "./ENM"

#extract methods list
allmet<-list.dirs(methods_folder)
presence_data_folder<- allmet[2]
allmet<- allmet[-1]
allmet<-gsub("./ENM/", "", allmet)

#extract ID list
all_IDs<-list.files(presence_data_folder)
all_IDs<-gsub("_metadata.txt","",all_IDs)
all_IDs<-gsub("\\.asc","",all_IDs)
all_IDs<-gsub("_ann.Rdata","",all_IDs)
all_IDs<-gsub("_svm.Rdata","",all_IDs)
all_IDs<-unique(all_IDs)

# ensamble data preparation
for (ID in all_IDs){
  
  cat("\nProcessing ",ID,"\n")
  cat("\nOccurrence enrichment and grid preparation\n")
  layers<- c()
  for(method in allmet){
    
    layer <-paste0("ENM","/",method,"/",ID,".asc")  
    layers<-append(layers, layer)
  }
  
  raster_out_file_name = paste0(output_folder,"/",ID,".asc") 
  if(file.exists(raster_out_file_name)){
    cat("Skipping\n")
    next
  }
  
  #set resolution from metadata
  BMPFile<- paste0("ENM","/",method,"/",ID,"_metadata.txt")
  BMPprops <- read.properties(BMPFile)
  resolution<- as.double(BMPprops$"Spatial resolution")
  
  #occurrence enrichment and grid preparation
  all_asc_files<-layers
  first_raster_data<-NA
  grid_of_points<-c()
  coordinate_at_res<-function(origin, coordinate, resolution){
    times<-round((coordinate-origin)/resolution)
    coordinate<-(times*resolution)+origin
    return (coordinate)
  }
  input_column_names<-c()
  models<-allmet
  
  model<-0
  for (file in all_asc_files){
    model<-model+1
    cat("Adding parameter:",models[model],"\n")
    asc_file<-raster(file)
    min_v_in_raster<-0
    max_v_in_raster<-cellStats(asc_file, max)
    min_x_in_raster<-asc_file@extent[1]
    max_x_in_raster<-asc_file@extent[2]
    min_y_in_raster<-asc_file@extent[3]
    max_y_in_raster<-asc_file@extent[4]
    
    #use the first input layer to setup the grid
    if (length(grid_of_points)==0){
      #generate a grid of points at the input resolution and with the extent of the first layer
      first_raster_data<-asc_file
      xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
      yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
      grid_of_points<-expand.grid(x = xseq, y = yseq) #combine the x and y coordinate to generate pairs
      
      #initialise the projection grid of parameters with just the coordinates
      grid_of_points_enriched<-grid_of_points
      }
      
    
    #extract raster values for the observations and the grid
    grid_of_points_extracted<-extract(x=asc_file,y=grid_of_points,method='simple')
    #enrich the columns of the training set and the grid
    column_name<-models[model]
    input_column_names<-c(input_column_names,column_name)
 
    grid_of_points_enriched[column_name]<-grid_of_points_extracted
  }

 #ensamble with threshold
  cat("\nAssessing presence by threshold and sum calculation\n")
  
  TSs<-c()
  for(method in allmet){
    
    #thresholds
    BMPFile<- paste0("ENM","/",method,"/",ID,"_metadata.txt")
    BMPprops <- read.properties(BMPFile)
    TS<- as.double(BMPprops$"Optimal decision threshold")
    
    TSs<-append(TSs, TS)
  }
  
  data <- grid_of_points_enriched
  
  counter=1
  for(threshold in TSs){
    
    #columns creation
    newcolumn <- sapply(1:length(data[,2+counter]), function(idx) {
      val<-data[,2+counter][idx]
      if(!is.na(val)){
        if (val ==-9999) return (-9999)
        if (val < threshold) return (as.double(0)) 
        else return (as.double(1))
      }else return(-9999)
    },simplify = T)
    data <- cbind(data, newcolumn)
    counter = counter+1
  }
  
  data_models<-data[,-c(1:(length(TSs)+2))]
  sums<-sapply(1:dim(data_models)[1], function(idx){
    data_row<-data_models[idx,]
    naindex<-which(data_row == -9999)
    if (length(naindex) == length(data_row))
      return(-9999)
    else{
      data_row[naindex]<-0
      return (sum(data_row))  
    }
    
  },simplify = T)
  
  data$sum <- sums
  data2<-data[,c("x","y","sum")]
  colnames(data2)<-c("longitude","latitude","sum_ts")
  
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
  extent(ro)<-extent(first_raster_data)
  
  #populate the matrix
  values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
  row_counter<-1
  for (y_c in 1:(nrow_r)){
    yp<-ypoints[y_c]
    row_rast<-data2[which(grid_of_points$y == yp),]
    row_rast<-row_rast[order(row_rast$longitude),]
    values[(nrow_r-row_counter+1),]<-row_rast$sum_ts[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  values_vec<-as.vector(t(values))
  values(ro)<-values_vec
  NAvalue(ro)<- -9999
  #plot(ro)
  
  #save the raster
  writeRaster(ro, filename=paste0(output_folder,"/ensemble_",ID), format="ascii",overwrite=TRUE) 

}