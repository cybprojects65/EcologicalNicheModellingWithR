rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(beepr)

#GENERAL PARAMETERS
AM_data_folder<-"./2100/AquaMaps/Fish"
ME_data_folder<-"./2100/MaxEnt/Fish"
ANN_data_folder<-"./2100/ArtificialNeuralNetworks/Fish"
SVM_data_folder<-"./2100/SupportVectorMachines/Fish"
presence_data_folder<-"./Presence Records/Fish"
output_folder_full<-"./output/output_not_norm_85_2100/Fish"


allspp<-list.files(presence_data_folder)
allspp<-gsub("Presence_","",allspp)
allspp<-gsub(".csv","",allspp)
#allspp<-gsub("_"," ",allspp)
species <-"Abudefduf_saxatilis"
for (species in allspp){
  
  AM_layer<- paste0(AM_data_folder,"/",species,".asc")
  ME_layer<-paste0(ME_data_folder,"/",species,"_native.asc")
  SVM_layer <- paste0(SVM_data_folder,"/",species,"_rcp85_2100.asc")
  ANN_layer <-paste0(ANN_data_folder,"/",species,".asc")
  
  layers<-c(ME_layer,SVM_layer,ANN_layer,AM_layer)
  
  
  cat("\nProcessing species ",species,"\n")
  raster_out_file_name = paste0(output_folder,"/",species,".asc")
  if(file.exists(raster_out_file_name)){
    cat("Skipping\n")
    next
  }
  

  resolution<-0.5 #deg

  
  
  #OCCURRENCE ENRICHMENT AND GRID PREPARATION
  cat("Step 1: Occurrence enrichment and grid preparation\n")

  all_asc_files<-layers
  first_raster_data<-NA
  grid_of_points<-c()
  coordinate_at_res<-function(origin, coordinate, resolution){
    
    times<-round((coordinate-origin)/resolution)
    coordinate<-(times*resolution)+origin
    return (coordinate)
  }
  input_column_names<-c()
  modelli<-c("ME","SVM","ANN","AM")
  modello<-0
  file<- all_asc_files[2]
  for (file in all_asc_files){
    modello<-modello+1
    cat("Adding parameter:",modelli[modello],"\n")
    asc_file<-raster(file)
    #plot(asc_file)
    #min_v_in_raster<-cellStats(asc_file, min)
    min_v_in_raster<-0
    max_v_in_raster<-cellStats(asc_file, max)
    min_x_in_raster<-asc_file@extent[1]
    max_x_in_raster<-asc_file@extent[2]
    min_y_in_raster<-asc_file@extent[3]
    max_y_in_raster<-asc_file@extent[4]
    
    #normalise the raster
    #asc_file_norm<-(asc_file - min_v_in_raster) / (max_v_in_raster - min_v_in_raster)
    
    #use the first input layer to setup the grid
    if (length(grid_of_points)==0){
      #generate a grid of points at the input resolution and with the extent of the first layer
      first_raster_data<-asc_file
      xseq<-seq(from=min_x_in_raster,to=max_x_in_raster,by=resolution)
      yseq<-seq(from=min_y_in_raster,to=max_y_in_raster,by=resolution)
      grid_of_points<-expand.grid(x = xseq, y = yseq)#combine the x and y coordinate to generate pairs
      
      #initialise the projection grid of parameters with just the coordinates
      grid_of_points_enriched<-grid_of_points
      }
      
    
    #extract raster values for the observations and the grid
    grid_of_points_extracted<-extract(x=asc_file,y=grid_of_points,method='simple')
    #grid_of_points_extracted<-extract(x=asc_file_norm,y=grid_of_points,method='simple')
    #enrich the columns of the training set and the grid
    column_name<-modelli[modello]
    input_column_names<-c(input_column_names,column_name)
 
    grid_of_points_enriched[column_name]<-grid_of_points_extracted
  }
  
  #adding normalized sum column
    # grid_of_points_enriched$sum <- rowSums(grid_of_points_enriched[,c("ME", "ANN","SVM","AM")], na.rm=TRUE)
    # max_in_sum <- max(grid_of_points_enriched$sum)
    # min_in_sum <- min(grid_of_points_enriched$sum)
    # grid_of_points_enriched$sum_norm <- (grid_of_points_enriched$sum - min_in_sum) / (max_in_sum - min_in_sum)
    

    
    # grid_of_points_only_normsum <-grid_of_points_enriched[,c("x","y","sum_norm")]
    # colnames(grid_of_points_only_normsum)<-c("longitude","latitude","ensamble_norm")
    # write.csv(grid_of_points_only_normsum,paste0(output_folder_csv,"/ensemble_",species,".csv"),row.names = FALSE)
  
    write.csv(grid_of_points_enriched,paste0(output_folder_full,"/ensemble_",species,".csv"),row.names = FALSE)
    

  
}
beep(4)