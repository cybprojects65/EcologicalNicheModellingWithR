rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(e1071)



RCP<- "85"
FNF <- "Non Fish"
year <- "2100"
FNF2 <-"NF"


support_vector_machines<-paste0("./ENMs/",FNF,"/SVM")   
presence_data_folder<-paste0("./Presence and Absence records/Presence Records/",FNF) 
base_environmental_layers<-"./Environmental Parameters/Res 05/Historical 2019"


projection_environmental_layers<-paste0("./Environmental Parameters/Res 05/RCP ",RCP,"/",year)
#projection_environmental_layers<-"./Environmental Parameters/Res 05/Historical 2019" #####test line ####

#output_folder
output_folder<-paste0("./ENMs_RCP",RCP,"_",year,"_",FNF2)  

#GENERAL PARAMETERS
resolution<-0.5 #deg
nativedistribution<-T
grid_of_points<-c()
input_column_names<-c()
all_asc_files<-list.files(path=projection_environmental_layers)

all_asc_files_original<-list.files(path=base_environmental_layers)

first_raster_data<-NA

coordinate_at_res<-function(origin, coordinate, resolution){
  
  times<-round((coordinate-origin)/resolution)
  coordinate<-(times*resolution)+origin
  return (coordinate)
}
env_p_original_idx<-1
#read all data
for (file in all_asc_files){
  cat("Adding parameter:",file,"\n")
  asc_file<-raster(paste0(projection_environmental_layers,"/",file))
  original_env_data_file<-raster(paste0(base_environmental_layers,"/",all_asc_files_original[env_p_original_idx]))
  
  #plot(asc_file)
  #min_v_in_raster<-cellStats(asc_file, min)
  #max_v_in_raster<-cellStats(asc_file, max)
  min_v_in_raster<-cellStats(original_env_data_file, min)
  max_v_in_raster<-cellStats(original_env_data_file, max)
  
  min_x_in_raster<-asc_file@extent[1]
  max_x_in_raster<-asc_file@extent[2]
  min_y_in_raster<-asc_file@extent[3]
  max_y_in_raster<-asc_file@extent[4]
  #normalise the raster
  asc_file_norm<-(asc_file - min_v_in_raster) / (max_v_in_raster - min_v_in_raster)
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
  grid_of_points_extracted<-extract(x=asc_file_norm,y=grid_of_points,method='simple')
  #enrich the columns of the training set and the grid
  column_name<-gsub(".asc","",file)
  input_column_names<-c(input_column_names,tolower(column_name))
  grid_of_points_enriched[column_name]<-grid_of_points_extracted
  
  env_p_original_idx<-env_p_original_idx+1
}
#rename the input columns of the training set and the grid to f1,..,fn
input_column_names_codes<-seq(1:(length(input_column_names)))
input_column_names_codes<-paste0("f",input_column_names_codes)
names(grid_of_points_enriched)<-c("x","y",input_column_names_codes)
grid_of_points_enriched_proj<-subset(grid_of_points_enriched, select = -c(x, y))

allspp<-list.files(presence_data_folder)
allspp<-gsub("Presence_","",allspp)
allspp<-gsub(".csv","",allspp)
allspp<-gsub("_"," ",allspp)
asc_files<-list()



for (species in allspp){
  
  cat("\nProcessing species ",species,"\n")
  raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) ,"_rcp",RCP,"_",year,".asc")
  svm_out_file_name = paste0(support_vector_machines,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_svm.Rdata")
  metadata_file<-paste0(support_vector_machines,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_metadata.txt")
  metadata_output_file<-paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_metadata.txt")
  
  load(svm_out_file_name)
  
  if(file.exists(raster_out_file_name)){
    cat("Skipping\n")
    next
  }
  
  #OCCURRENCE ENRICHMENT AND GRID PREPARATION
  cat("Step 1: Occurrence enrichment and grid preparation\n")
  presence_file<-paste0(paste("Presence",gsub(" ","_",species),sep = "_"),".csv")
  presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
  presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = presence$longitude, resolution = resolution)
  presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = presence$latitude, resolution = resolution)
  
  #project the optimal SVM on the grid
  probability_on_grid_of_points_enriched_proj<-predict(svmfit, grid_of_points_enriched_proj, probability = TRUE,na.action = na.exclude)
  
  min_v_in_probability<-min(probability_on_grid_of_points_enriched_proj,na.rm = TRUE)
  max_v_in_probability<-max(probability_on_grid_of_points_enriched_proj,na.rm = TRUE)
  probability_on_grid_of_points_enriched_proj_N<-(probability_on_grid_of_points_enriched_proj - min_v_in_probability) /(max_v_in_probability-min_v_in_probability)
  
  
  
  grid_of_pointssp<-grid_of_points
  grid_of_pointssp$probability<-probability_on_grid_of_points_enriched_proj_N
  grid_of_pointssp[rowSums(is.na(grid_of_pointssp)) > 0,3]<- -9999
  
  #WRITE METADATA AND OUTPUT
  fileConn<-file.copy(from = metadata_file, to = metadata_output_file)
  
  if (nativedistribution){
    cat("\nSTEP 2b: Adjusting the projection for native distribution...\n")
    #calculate average distances between the samples
    distances<-sapply(1:length(presence$longitude_res), function(i){
      p_long<-presence$longitude_res[i]
      p_lat<-presence$latitude_res[i]
      dx<-(presence$longitude_res-p_long)
      dy<-(presence$latitude-p_lat)
      sqrd<-(dx*dx)+(dy*dy)
      d<-mean(sqrt(sqrd))
      return (d)
    },simplify = T)
    sigma<-sd(distances)
    #get minimum distances between the grid and the presence samples
    mindistance<-sapply(1:length(grid_of_pointssp$probability), function(i){
      p_long<-grid_of_pointssp$x[i]
      p_lat<-grid_of_pointssp$y[i]
      dx<-(presence$longitude_res-p_long)
      dy<-(presence$latitude-p_lat)
      sqrd<-(dx*dx)+(dy*dy)
      d<-min(sqrt(sqrd))
      return (d)
    },simplify = T)
    #setup probability weights
    weights<-dnorm(mindistance,0,sigma)/dnorm(0,0,sigma)
    weights[which(mindistance<=sigma)]<-1
    #weight probability by distance from presence points
    na_points<-which(grid_of_pointssp$probability==-9999)
    grid_of_pointssp$probability[-na_points]<-grid_of_pointssp$probability[-na_points]*weights[-na_points]
  }
  
  #WRITE THE ASC file
  cat("\nSTEP 3: Projecting the ANN...\n")
  #raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , ".asc")
  #define the grid points to populate the dataset
  ypoints<-unique(grid_of_pointssp$y)
  xpoints<-unique(grid_of_pointssp$x)
  ncol_r<-length(xpoints)-1
  nrow_r<-length(ypoints)-1
  #create a new raster with the same extent and resolution of the first layer
  ro <- raster(ncol=ncol_r, nrow=nrow_r)
  length(values(ro))
  
  res(ro) <- resolution
  length(values(ro))
  extent(ro)<-extent(first_raster_data)
  #populate the matrix
  values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
  row_counter<-1
  for (y_c in 2:(nrow_r+1)){
    yp<-ypoints[y_c]
    row_rast<-grid_of_pointssp[which(grid_of_points$y == yp),]
    row_rast<-row_rast[order(row_rast$x),]
    values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  values_vec<-as.vector(t(values))
  values(ro)<-values_vec
  NAvalue(ro)<- -9999
  #plot(ro)
  
  #save the raster
  cat("Writing the output..\n")
  writeRaster(ro, filename=raster_out_file_name, format="ascii",overwrite=TRUE)
  cat("Species",species,"done.\n")
  #break
}