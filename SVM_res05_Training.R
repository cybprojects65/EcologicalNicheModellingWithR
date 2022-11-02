rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(e1071)

#GENERAL PARAMETERS
env_data_folder<-"./Environmental Parameters/Res 05/Historical 2019"
presence_data_folder<-"./Presence and Absence records/Presence Records/Non Fish"
absence_data_folder<-"./Presence and Absence records/Absence Records/Non Fish"
output_folder<-"./ENMs/Non Fish"

allspp<-list.files(presence_data_folder)
allspp<-gsub("Presence_","",allspp)
allspp<-gsub(".csv","",allspp)
allspp<-gsub("_"," ",allspp)



for (species in allspp){
  
  cat("\nProcessing species ",species,"\n")
  raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , ".asc")
  svm_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_svm.Rdata")
  if(file.exists(raster_out_file_name)){
    cat("Skipping\n")
    next
  }
  
  ths_sample<-100
  abs_sample<-10000
  resolution<-0.5 #deg
  costlist <- c(0.001,0.01,0.1,1,10,100)
  kernellist <- c("linear", "polynomial","radial", "sigmoid")
  kernellistwp <- c("linear", "radial", "sigmoid")
  coef0list <- c(0,1) #used in "polynomial"and "sigmoid"
  gammalist <-c(0.0001, 0.001, 0.01, 0.1, 1 ,10) #used in "polynomial","radial", "sigmoid"
  degreelist <- c(3,4)  #used in "polynomial"
  cross95 <- 20 # nfold in ANN
  do_the_shrinking <- FALSE #I will try to do not use data shrinking
  nativedistribution <- TRUE
  
  
  #OCCURRENCE ENRICHMENT AND GRID PREPARATION
  cat("Step 1: Occurrence enrichment and grid preparation\n")
  presence_file<-paste0(paste("Presence",gsub(" ","_",species),sep = "_"),".csv")
  absence_file<-paste0(paste("Absence",gsub(" ","_",species),sep = "_"),".csv")
  all_asc_files<-list.files(path=env_data_folder)
  first_raster_data<-NA
  grid_of_points<-c()
  coordinate_at_res<-function(origin, coordinate, resolution){
    
    times<-round((coordinate-origin)/resolution)
    coordinate<-(times*resolution)+origin
    return (coordinate)
  }
  input_column_names<-c()
  
  for (file in all_asc_files){
    cat("Adding parameter:",file,"\n")
    asc_file<-raster(paste0(env_data_folder,"/",file))
    #plot(asc_file)
    min_v_in_raster<-cellStats(asc_file, min)
    max_v_in_raster<-cellStats(asc_file, max)
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
      #read absence and presence file
      absence<-read.table(paste0(absence_data_folder,"/",absence_file),header = TRUE, sep=",")
      
      indexa    <- sample(1:dim(grid_of_points_enriched)[1], abs_sample)
      absence1<-grid_of_points_enriched[indexa,]
      absence1$species<-species
      absence1<-data.frame(species=absence1$species,longitude=absence1$x, latitude = absence1$y)
      absence<-absence1#rbind(absence,absence1)
      
      presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
      #put the observations at the input resolution
      absence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = absence$longitude, resolution = resolution)
      absence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = absence$latitude, resolution = resolution)
      presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = presence$longitude, resolution = resolution)
      presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = presence$latitude, resolution = resolution)
      #prepare an x,y string for fast indexing
      absencexy<-paste0(absence$longitude_res,",",absence$latitude_res)
      presencexy<-paste0(presence$longitude_res,",",presence$latitude_res)
      #take the intersection and remove the shared lines
      shared_coordinates_abs<-which(absencexy %in% presencexy)
      shared_coordinates_pres<-which(presencexy %in% absencexy)
      presence_disjoint<-presence
      absence_disjoint<-absence
      if (length(shared_coordinates_pres)>0)
        presence_disjoint<-presence[-shared_coordinates_pres,]
      if (length(shared_coordinates_abs)>0)
        absence_disjoint<-absence[-shared_coordinates_abs,]
      
      if (dim(absence_disjoint)[1]>ths_sample){
        index    <- sample(1:dim(absence_disjoint)[1], ths_sample)
        absence_disjoint <- absence_disjoint[index, ]
      }
      if (dim(presence_disjoint)[1]>ths_sample){
        index    <- sample(1:dim(presence_disjoint)[1], ths_sample)
        presence_disjoint <- presence_disjoint[index, ]
      }
      #prepare the presence-absence data frame as the training set
      coordinates_to_enrich<-rbind(data.frame(x=presence_disjoint$longitude,y=presence_disjoint$latitude),
                                   data.frame(x=absence_disjoint$longitude,y=absence_disjoint$latitude))
      #initialise the training set with just the coordinates
      training_set<-coordinates_to_enrich
    }
    
    #plot(asc_file_norm)
    
    #extract raster values for the observations and the grid
    enriched_with_asc<-extract(x=asc_file_norm,y=coordinates_to_enrich,method='simple')
    grid_of_points_extracted<-extract(x=asc_file_norm,y=grid_of_points,method='simple')
    #enrich the columns of the training set and the grid
    column_name<-gsub(".asc","",file)
    input_column_names<-c(input_column_names,tolower(column_name))
    training_set[column_name]<-enriched_with_asc
    grid_of_points_enriched[column_name]<-grid_of_points_extracted
  }
  #rename the input columns of the training set and the grid to f1,..,fn
  input_column_names_codes<-seq(1:(length(input_column_names)))
  input_column_names_codes<-paste0("f",input_column_names_codes)
  names(grid_of_points_enriched)<-c("x","y",input_column_names_codes)
  training_set$t<-0
  training_set$t[1:dim(presence_disjoint)[1]]<-1
  training_set<-na.omit(training_set) #delete rows that contain at least one NA
  names(training_set)<-c("x","y",input_column_names_codes,"t")
  
  
  cat("Step 2: SVM training with multiple parametrization\n")
  
  training_set_features_only<-subset(training_set, select = -c(x, y))
  grid_of_points_enriched_features_only<-subset(grid_of_points_enriched, select = -c(x, y),)
  
  grid_of_points_enriched
  f <- as.formula(paste("t", "~", paste(input_column_names_codes, collapse = " + ") ))
  set.seed(20)
  
  #initialise features that will record the optimal scores
  opt_cost<- "NA"
  opt_gamma<- "NA"
  opt_coef0 <- "NA"
  opt_degree <- "NA"
  opt_kernel <- ""
  
  opt_accuracy<-0
  opt_accuracy_self<-0
  opt_threshold<-0
 
  
  tuned <- tune(svm, t ~., data = training_set_features_only,ranges=list(kernel=kernellist, cost = costlist, gamma = gammalist, coef0=coef0list,degree=degreelist), scale = TRUE)
  print(tuned)
  
  
  
  opt_kernel <- tuned$best.parameters$kernel
  
  if(opt_kernel == "linear"){
    opt_kernel <- "linear"
    print(" linear kernel on!")
    opt_cost<- tuned$best.parameters$cost
    svmfit <- svm(t ~., data = training_set_features_only, kernel = "linear", cost = opt_cost , scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
  }
  
  if(opt_kernel == "polynomial"){
    opt_kernel <- "polynomial"
    print(" polynomial kernel on!")
    opt_cost<- tuned$best.parameters$cost
    opt_gamma<- tuned$best.parameters$gamma
    opt_degree <- tuned$best.parameters$degree
    opt_coef0 <- tuned$best.parameters$coef0
    
    svmfit <- svm(t ~., data = training_set_features_only, kernel = "polynomial", cost = opt_cost ,gamma= opt_gamma,degree=opt_degree,coef0= opt_coef0, scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
  }
  
  if(opt_kernel == "radial"){
    opt_kernel <-"radial"
    print(" radial kernel on!")
    opt_cost<- tuned$best.parameters$cost
    opt_gamma <- tuned$best.parameters$gamma
    
    svmfit <- svm(t ~., data = training_set_features_only, kernel = "radial", cost = opt_cost ,gamma= opt_gamma, scale = FALSE, cross=cross95,probability=TRUE,shrinking=do_the_shrinking) #model call
  }
  
  if(opt_kernel == "sigmoid"){
    opt_kernel <- "sigmoid"
    print(" sigmoid kernel on!")
    opt_cost<- tuned$best.parameters$cost
    opt_gamma<- tuned$best.parameters$gamma
    opt_coef0 <- tuned$best.parameters$coef0
    
    svmfit <- svm(t ~., data = training_set_features_only, kernel = "sigmoid", cost = opt_cost ,gamma= opt_gamma,coef0= opt_coef0, scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
  }
  
  print(svmfit)
  #add probability from model to test set
  training_set_features_only_w_probability <- training_set_features_only
  prediction_on_training_set <- predict(svmfit, training_set_features_only_w_probability, probability = TRUE)
  training_set_features_only_w_probability$probability<-prediction_on_training_set
  
  #project the optimal ANN on the grid
  grid_of_points_enriched_features_only_w_probability<-grid_of_points_enriched_features_only
  prediction_on_grid_of_point_enriched <- predict(svmfit, grid_of_points_enriched_features_only_w_probability, probability = TRUE, na.action = na.exclude)
  #min max normalization on grid prediction
  min_v_in_raster<-min(prediction_on_grid_of_point_enriched,na.rm = TRUE)
  max_v_in_raster<-max(prediction_on_grid_of_point_enriched,na.rm = TRUE)
  prediction_on_grid_of_point_enriched_N<-(prediction_on_grid_of_point_enriched - min_v_in_raster) /(max_v_in_raster-min_v_in_raster)
  
  grid_of_points_enriched_features_only_w_probability$probability<-prediction_on_grid_of_point_enriched_N
  grid_of_points_w_probability <- grid_of_points
  grid_of_points_w_probability$probability<-prediction_on_grid_of_point_enriched_N
  grid_of_points_w_probability[rowSums(is.na(grid_of_points_w_probability)) > 0,3]<- -9999
  
  
  #tuning the best threshold #######################################
  training_set_features_only_w_probability_N <- training_set_features_only
  prediction_on_training_set_N <- (prediction_on_training_set - min_v_in_raster) /(max_v_in_raster-min_v_in_raster)
  training_set_features_only_w_probability_N$probability<-prediction_on_training_set_N
  
  decision_thresholds<-as.numeric(quantile( training_set_features_only_w_probability_N$probability,probs=c(0.01,0.05,0.1,0.2,0.5,0.8))) #take 1% of the training set out
  opt_accuracy<-0
  opt_threshold<-0
  #loop on decision thresholds testing
  for (decision_threshold in decision_thresholds){
    training_set_features_only_w_probability_N$detected<-F
    training_set_features_only_w_probability_N$detected[which( ( training_set_features_only_w_probability_N$t ==0 & training_set_features_only_w_probability_N$probability<decision_threshold) | 
                                                                 (training_set_features_only_w_probability_N$t ==1 & training_set_features_only_w_probability_N$probability>decision_threshold) ) ]<-T
    accuracy<- length(which(training_set_features_only_w_probability_N$detected))/length(training_set_features_only_w_probability_N$detected)
    #record the best threshold and gained accuracy
    if (opt_accuracy<accuracy){
      opt_accuracy=accuracy
      opt_threshold<-decision_threshold
    }
    #cat("Accuracy selftest =", accuracy*100,"%","(thr:",decision_threshold,")","\n")
  }#end loop on decision thresholds###################################
  cat("Accuracy selftest best result =", opt_accuracy*100,"%","(thr:",opt_threshold,")","\n")
  
  
  
  #WRITE METADATA
   fileConn<-file(paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_metadata.txt"))
  writeLines(c(
    paste0("Species name = ",species),
    paste0("Spatial resolution = ",resolution),
    paste0("Number of folders of cross validation = ",cross95),
    paste0("Shrinking-heuristics = ",do_the_shrinking),
    paste0("Optimal overall accuracy selftest = ", opt_accuracy*100),
    paste0("Optimal decision threshold = ",opt_threshold),
    paste0("Optimal kernel = ",opt_kernel),
    paste0("Optimal cost = ",opt_cost),
    paste0("Optimal gamma = ",opt_gamma),
    paste0("Optimal coef0 = ",opt_coef0),
    paste0("Optimal degree = ",opt_degree)
  ), fileConn)
  close(fileConn)
  
  # Native Distribution
  if (nativedistribution){
    cat("\nSTEP 2b: Adjusting the projection for native distribution...\n")
    #calculate average distances between the samples
    distances<-sapply(1:length(presence$longitude_res), function(i){
      p_long<-presence$longitude_res[i]
      p_lat<-presence$latitude_res[i]
      dx<-(presence$longitude_res-p_long)
      dy<-(presence$latitude_res-p_lat)
      sqrd<-(dx*dx)+(dy*dy)
      d<-mean(sqrt(sqrd))
      return (d)
    },simplify = T)
    sigma<-sd(distances)
    #get minimum distances between the grid and the presence samples
    mindistance<-sapply(1:length(grid_of_points_w_probability$probability), function(i){
      p_long<-grid_of_points_w_probability$x[i]
      p_lat<-grid_of_points_w_probability$y[i]
      dx<-(presence$longitude_res-p_long)
      dy<-(presence$latitude_res-p_lat)
      sqrd<-(dx*dx)+(dy*dy)
      d<-min(sqrt(sqrd))
      return (d)
    },simplify = T)
    #setup probability weights
    weights<-dnorm(mindistance,0,sigma)/dnorm(0,0,sigma)
    weights[which(mindistance<=sigma)]<-1
    #weight probability by distance from presence points
    na_points<-which(grid_of_points_w_probability$probability==-9999)
    grid_of_points_w_probability$probability[-na_points]<-grid_of_points_w_probability$probability[-na_points]*weights[-na_points]
  }
  
  #WRITE THE ASC file
  cat("\nSTEP 3: Projecting the SVN...\n")
  
  ypoints<-unique(grid_of_points$y)
  xpoints<-unique(grid_of_points$x)
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
    row_rast<-grid_of_points_w_probability[which(grid_of_points$y == yp),]
    row_rast<-row_rast[order(row_rast$x),]
    values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
    row_counter<-row_counter+1
  }
  values_vec<-as.vector(t(values))
  
  values(ro)<-values_vec
  NAvalue(ro)<- -9999
  
  #save the raster
  cat("Writing the output..\n")
  writeRaster(ro, filename=raster_out_file_name, format="ascii",overwrite=TRUE)
  
  save(svmfit, file=svm_out_file_name)
  cat("Species",species,"done.\n")
  #break
}