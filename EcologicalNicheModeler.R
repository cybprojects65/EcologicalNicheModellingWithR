rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)
library(e1071)
library(properties)

Workflow_Parameters_File<- paste0("./Workflow_Parameters.txt")
WPFprops <- read.properties(Workflow_Parameters_File)

#GENERAL PARAMETERS
env_data_folder<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$env_data_folder)
presence_data_folder<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$presence_data_folder)
training_folder<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$output_folder)
output_folder<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$output_folder)
projection_environmental_layers<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$projection_environmental_layers)

dir.create(output_folder, showWarnings = FALSE)

ANN_Active <- as.logical(WPFprops$"ANN_Active")
SVM_Active <- as.logical(WPFprops$"SVM_Active")
AQUAMAPS_Active <- as.logical(WPFprops$"AQUAMAPS_Active")
MAXENT_Active <- as.logical(WPFprops$"MAXENT_Active")
Create_ANN_Rdata <- as.logical(WPFprops$"Create_ANN_Rdata")  
Create_SVM_Rdata <- as.logical(WPFprops$"Create_SVM_Rdata")
Create_AquaMaps_Rdata <- as.logical(WPFprops$"Create_AquaMaps_Rdata")
Create_MAXENT_Rdata <- as.logical(WPFprops$"Create_MAXENT_Rdata")
model_projection <- as.logical(WPFprops$"model_projection")
maxent_prevalence<-as.numeric(WPFprops$"maxent_prevalence")

if (model_projection){
  output_folder<-gsub(pattern = "\\\"",replacement = "", x=WPFprops$output_folder_projection)
  dir.create(output_folder, showWarnings = FALSE)
  Create_ANN_Rdata <- F
  Create_SVM_Rdata <- F
  Create_AquaMaps_Rdata <- F
  Create_MAXENT_Rdata <- F
}

support_vector_machines_folder <-paste0(training_folder,"/SVM/")
trained_neural_networks <-paste0(training_folder,"/ANN/")
trained_aquamaps <-paste0(training_folder,"/AquaMaps/")
trained_maxent <-paste0(training_folder,"/MaxEnt/")

#general models parameters
ths_sample<- as.double(WPFprops$"ths_sample")
abs_sample<- as.double(WPFprops$"abs_sample")
nativedistribution <- as.logical(WPFprops$"nativedistribution") 

#parameters SVM
costlist <- as.numeric(unlist(strsplit(WPFprops$"costlist", split=",")))
kernellist <- unlist(strsplit(WPFprops$"kernellistx", split=",")) 
kernellistwp <- unlist(strsplit(WPFprops$"kernellistwp", split=","))
coef0list <- as.numeric(unlist(strsplit(WPFprops$"coef0list", split=","))) #used in "polynomial"and "sigmoid"
gammalist <-as.numeric(unlist(strsplit(WPFprops$"gammalist", split=","))) #used in "polynomial","radial", "sigmoid"
degreelist <- as.numeric(unlist(strsplit(WPFprops$"degreelist", split=",")))  #used in "polynomial"
cross95 <- as.numeric(WPFprops$"cross95") # nfold in ANN
do_the_shrinking <- as.logical(WPFprops$"do_the_shrinking")  #I will try to do not use data shrinking

#parameters ANN
rp <- as.numeric(WPFprops$"rp")     # number of repetitions for the training - #20 gives more stability to the network and more independence of the initial training step
thld <- as.numeric(WPFprops$"thld")   # threshold for minimum decrease in overall error, default 0.01 = 1%
stp  <- as.numeric(WPFprops$"stp")  # the maximum steps for the training of the neural network, default 1e+05
alg  <- WPFprops$"alg"  # possible: backprop, rprop+, rprop-, sag, or slr
act.fct <- WPFprops$"act.fct" #"logistic" # possible: "logistic" or "tanh"; linear.output must be 
nfold <- as.numeric(WPFprops$"nfold") #20
hiddens <- as.numeric(unlist(strsplit(WPFprops$"hiddens", split=",")))

#creation of output folder for any model
if(SVM_Active == TRUE){
  path<-paste0(output_folder,"/SVM")
  if(!dir.exists(path)) dir.create(path)
}

if(ANN_Active==TRUE){
  path<-paste0(output_folder,"/ANN")
  if(!dir.exists(path)) dir.create(path)
}

if(AQUAMAPS_Active==TRUE){
  path<-paste0(output_folder,"/AquaMaps")
  if(!dir.exists(path)) dir.create(path)
}

if(MAXENT_Active==TRUE){
  path<-paste0(output_folder,"/MaxEnt")
  if(!dir.exists(path)) dir.create(path)
}

#list of ID
all_IDs<-list.files(presence_data_folder)
all_IDs<-gsub("Presence_","",all_IDs)
all_IDs<-gsub(".csv","",all_IDs)
all_IDs<-gsub("_"," ",all_IDs)

#####
#OCCURRENCE ENRICHMENT AND GRID PREPARATION
cat("Step 1: Occurrence enrichment and grid preparation\n")
maxentcache<-paste0(env_data_folder,"/maxent.cache")
if (dir.exists(maxentcache))
  unlink(maxentcache,recursive = T)
maxentcache<-paste0(projection_environmental_layers,"/maxent.cache")
if (dir.exists(maxentcache))
  unlink(maxentcache,recursive = T)

all_asc_files<-list.files(path=env_data_folder)

if(model_projection == TRUE){
  cat("Using projected environmental parameters\n")
  all_asc_files <-list.files(path=projection_environmental_layers)
  all_asc_files_original<-list.files(path=env_data_folder)
  env_p_original_idx<-1
}

first_raster_data<-NA
grid_of_points<-c()
input_column_names<-c()

coordinate_at_res<-function(origin, coordinate, resolution){
  times<-round((coordinate-origin)/resolution)
  coordinate<-(times*resolution)+origin
  return (coordinate)
}


for (file in all_asc_files){
  cat("Adding parameter:",file,"\n")
  asc_file<-raster(paste0(env_data_folder,"/",file))
  if(model_projection == TRUE){
    asc_file<-raster(paste0(projection_environmental_layers,"/",file))
    original_env_data_file<-raster(paste0(env_data_folder,"/",all_asc_files_original[env_p_original_idx]))
  }
  #plot(asc_file)
  min_v_in_raster<-cellStats(asc_file, min)
  max_v_in_raster<-cellStats(asc_file, max)
  if(model_projection == TRUE){
    min_v_in_raster<-cellStats(original_env_data_file, min)
    max_v_in_raster<-cellStats(original_env_data_file, max) 
  }
  
  min_x_in_raster<-asc_file@extent[1]
  max_x_in_raster<-asc_file@extent[2]
  min_y_in_raster<-asc_file@extent[3]
  max_y_in_raster<-asc_file@extent[4]
  resolution <- res(asc_file)[1]
  #normalise the raster
  asc_file_norm<-(asc_file - min_v_in_raster) / (max_v_in_raster - min_v_in_raster)
  #use the first input layer to setup the grid
  if (length(grid_of_points)==0){
    #generate a grid of points at the input resolution and with the extent of the first layer
    first_raster_data<-asc_file
    xseq<-seq(from=min_x_in_raster+(resolution/2),to=max_x_in_raster-(resolution/2),by=resolution)
    yseq<-seq(from=min_y_in_raster+(resolution/2),to=max_y_in_raster-(resolution/2),by=resolution)
    grid_of_points<-expand.grid(x = xseq, y = yseq)#combine the x and y coordinate to generate pairs
    #initialise the projection grid of parameters with just the coordinates
    grid_of_points_enriched<-grid_of_points
  }
  
  #plot(asc_file_norm)
  
  #extract raster values for the observations and the grid
  
  grid_of_points_extracted<-extract(x=asc_file_norm,y=grid_of_points,method='simple')
  #enrich the columns of the training set and the grid
  column_name<-gsub(".asc","",file)
  input_column_names<-c(input_column_names,tolower(column_name))
  
  grid_of_points_enriched[column_name]<-grid_of_points_extracted
  if(model_projection == TRUE){
    env_p_original_idx<-env_p_original_idx+1
  }
}
#rename the input columns of the training set and the grid to f1,..,fn
input_column_names_codes<-seq(1:(length(input_column_names)))
input_column_names_codes<-paste0("f",input_column_names_codes)
names(grid_of_points_enriched)<-c("x","y",input_column_names_codes)

grid_of_points_enriched_proj<-subset(grid_of_points_enriched, select = -c(x, y))



#####

for (ID in all_IDs){
  cat("\nProcessing ID ",ID,"\n")
  presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
  raster_out_file_name_ANN = paste0(output_folder,"/ANN/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
  raster_out_file_name_SVM = paste0(output_folder,"/SVM/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
  raster_out_file_name_AquaMaps = paste0(output_folder,"/AquaMaps/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
  raster_out_file_name_MaxEnt = paste0(output_folder,"/MaxEnt/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
  
  # if(file.exists(raster_out_file_name_ANN)||file.exists(raster_out_file_name_SVM)){
  #   cat("Skipping\n")
  #   next
  # }
  ####Grid preparation
  
  if(model_projection == FALSE){
	if (dim(grid_of_points_enriched)[1]<abs_sample)
      abs_sample=ths_sample											   
    indexa    <- sample(1:dim(grid_of_points_enriched)[1], abs_sample)
    absence1<-grid_of_points_enriched[indexa,]
    absence1$ID<-ID
    absence1<-data.frame(ID=absence1$ID,longitude=absence1$x, latitude = absence1$y)
    absence<-absence1#rbind(absence,absence1)
    presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
    presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster+(resolution/2),coordinate = presence$longitude, resolution = resolution)
    presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster+(resolution/2),coordinate = presence$latitude, resolution = resolution)
    
    #put the observations at the input resolution
    absence$longitude_res<-coordinate_at_res(origin = min_x_in_raster+(resolution/2),coordinate = absence$longitude, resolution = resolution)
    absence$latitude_res<-coordinate_at_res(origin = min_y_in_raster+(resolution/2),coordinate = absence$latitude, resolution = resolution)
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
    coordinates_to_enrich<-rbind(data.frame(x=presence_disjoint$longitude_res,y=presence_disjoint$latitude_res),
                                 data.frame(x=absence_disjoint$longitude_res,y=absence_disjoint$latitude_res))
    #initialise the training set with just the coordinates
    training_set<-left_join(coordinates_to_enrich, grid_of_points_enriched, by=c("x"= "x","y"="y"))
    training_set$t<-0
    training_set$t[1:dim(presence_disjoint)[1]]<-1
    training_set<-na.omit(training_set) #delete rows that contain at least one NA
    names(training_set)<-c("x","y",input_column_names_codes,"t")
	if (length(which(training_set$t==1))==0)
    {
      cat("ERROR: no positive examples in the training set!\n")    
      stop("ERROR: no positive examples in the training set!\n")
    }
  }
  
  #####
  if(model_projection == TRUE){
    metadata_file<-paste0(support_vector_machines_folder,gsub(pattern = " ",replacement = "_",x = ID), "_metadata.txt")
    metadata_output_file<-paste0(output_folder,"/SVM/",gsub(pattern = " ",replacement = "_",x = ID) , "_metadata.txt")
    file.copy(from = metadata_file, to = metadata_output_file)
  }
  
  if(SVM_Active == TRUE){
    if(model_projection == TRUE){
      cat("Step 2: SVM projection session\n")
      raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) ,".asc")
      svm_out_file_name = paste0(support_vector_machines_folder,paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_svm.Rdata")
      load(svm_out_file_name)
      
      if(file.exists(raster_out_file_name)){
        cat("Skipping\n")
        next
      }
      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
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
      
    }else{
      cat("Step 2: SVM training with multiple parametrization\n")
      
      training_set_features_only_SVM<-subset(training_set, select = -c(x, y))
      grid_of_points_enriched_features_only<-subset(grid_of_points_enriched, select = -c(x, y),)
      
      
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
      
      
      tuned <- tune(svm, t ~., data = training_set_features_only_SVM,ranges=list(kernel=kernellist, cost = costlist, gamma = gammalist, coef0=coef0list,degree=degreelist), scale = FALSE)
      print(tuned)
      
      
      
      opt_kernel <- tuned$best.parameters$kernel
      
      if(opt_kernel == "linear"){
        opt_kernel <- "linear"
        print(" linear kernel on!")
        opt_cost<- tuned$best.parameters$cost
        svmfit <- svm(t ~., data = training_set_features_only_SVM, kernel = "linear", cost = opt_cost , scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
      }
      
      if(opt_kernel == "polynomial"){
        opt_kernel <- "polynomial"
        print(" polynomial kernel on!")
        opt_cost<- tuned$best.parameters$cost
        opt_gamma<- tuned$best.parameters$gamma
        opt_degree <- tuned$best.parameters$degree
        opt_coef0 <- tuned$best.parameters$coef0
        
        svmfit <- svm(t ~., data = training_set_features_only_SVM, kernel = "polynomial", cost = opt_cost ,gamma= opt_gamma,degree=opt_degree,coef0= opt_coef0, scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
      }
      
      if(opt_kernel == "radial"){
        opt_kernel <-"radial"
        print(" radial kernel on!")
        opt_cost<- tuned$best.parameters$cost
        opt_gamma <- tuned$best.parameters$gamma
        
        svmfit <- svm(t ~., data = training_set_features_only_SVM, kernel = "radial", cost = opt_cost ,gamma= opt_gamma, scale = FALSE, cross=cross95,probability=TRUE,shrinking=do_the_shrinking) #model call
      }
      
      if(opt_kernel == "sigmoid"){
        opt_kernel <- "sigmoid"
        print(" sigmoid kernel on!")
        opt_cost<- tuned$best.parameters$cost
        opt_gamma<- tuned$best.parameters$gamma
        opt_coef0 <- tuned$best.parameters$coef0
        
        svmfit <- svm(t ~., data = training_set_features_only_SVM, kernel = "sigmoid", cost = opt_cost ,gamma= opt_gamma,coef0= opt_coef0, scale = FALSE,probability=TRUE, cross=cross95,shrinking=do_the_shrinking) #model call
      }
      
      print(svmfit)
      #add probability from model to test set
      training_set_features_only_w_probability <- training_set_features_only_SVM
      prediction_on_training_set <- predict(svmfit, training_set_features_only_w_probability, probability = TRUE)
      training_set_features_only_w_probability$probability<-prediction_on_training_set
      
      #project the optimal SVM on the grid
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
      
      
      
      #tuning the best threshold 
      training_set_features_only_w_probability_N <- training_set_features_only_SVM
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
        
      }#end loop on decision thresholds###################################
      cat("Accuracy selftest best result =", opt_accuracy*100,"%","(thr:",opt_threshold,")","\n")
      
      
    }  #end else training/projection 
    #WRITE METADATA SVM
    if(model_projection == FALSE){
      fileConn<-file(paste0(output_folder,"/SVM/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt"))
      writeLines(c(
        paste0("ID name = ",ID),
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
    }
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
      if(model_projection == TRUE){
        mindistance<-sapply(1:length(grid_of_pointssp$probability), function(i){
          p_long<-grid_of_pointssp$x[i]
          p_lat<-grid_of_pointssp$y[i]
          dx<-(presence$longitude_res-p_long)
          dy<-(presence$latitude-p_lat)
          sqrd<-(dx*dx)+(dy*dy)
          d<-min(sqrt(sqrd))
          return (d)
        },simplify = T) 
      }else{
        mindistance<-sapply(1:length(grid_of_points_w_probability$probability), function(i){
          p_long<-grid_of_points_w_probability$x[i]
          p_lat<-grid_of_points_w_probability$y[i]
          dx<-(presence$longitude_res-p_long)
          dy<-(presence$latitude_res-p_lat)
          sqrd<-(dx*dx)+(dy*dy)
          d<-min(sqrt(sqrd))
          return (d)
        },simplify = T)
      }
      #setup probability weights
      weights<-dnorm(mindistance,0,sigma)/dnorm(0,0,sigma)
      weights[which(mindistance<=sigma)]<-1
      #weight probability by distance from presence points
      if(model_projection == TRUE){
        na_points<-which(grid_of_pointssp$probability==-9999)
        grid_of_pointssp$probability[-na_points]<-grid_of_pointssp$probability[-na_points]*weights[-na_points]
      }else{
        na_points<-which(grid_of_points_w_probability$probability==-9999)
        grid_of_points_w_probability$probability[-na_points]<-grid_of_points_w_probability$probability[-na_points]*weights[-na_points]
      }
    }
    
    #WRITE THE ASC file
    cat("\nSTEP 3: Projecting the SVM...\n")
    if(model_projection == TRUE){
      ypoints<-unique(grid_of_pointssp$y)
      xpoints<-unique(grid_of_pointssp$x) 
    }else{ 
      ypoints<-unique(grid_of_points$y)
      xpoints<-unique(grid_of_points$x)
    }
    ncol_r<-length(xpoints)
    nrow_r<-length(ypoints)
    #create a new raster with the same extent and resolution of the first layer
    ro <- raster(ncol=ncol_r, nrow=nrow_r)
    #length(values(ro))
    
    #res(ro) <- resolution
    #length(values(ro))
    extent(ro)<-extent(first_raster_data)
    #populate the matrix
    values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
    row_counter<-1
    for (y_c in 1:(nrow_r)){
      yp<-ypoints[y_c]
      if(model_projection == TRUE){
        row_rast<-grid_of_pointssp[which(grid_of_points$y == yp),]
        
      }else{
        row_rast<-grid_of_points_w_probability[which(grid_of_points$y == yp),]
      }
      row_rast<-row_rast[order(row_rast$x),]
      values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
      row_counter<-row_counter+1
    }
    values_vec<-as.vector(t(values))
    
    values(ro)<-values_vec
    NAvalue(ro)<- -9999
    
    #save the raster
    cat("Writing the output..\n")
    writeRaster(ro, filename=raster_out_file_name_SVM, format="ascii",overwrite=TRUE)
    svm_out_file_name = paste0(output_folder,"/SVM/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_svm.Rdata")
    if(Create_SVM_Rdata){
      save(svmfit, file=svm_out_file_name)  #this line save also the ANN to be used on another set of Environmental parameter
    }
    cat("ID",ID," SVM done.\n")
  }
  
  if(ANN_Active==TRUE){
    grid_of_points_enriched_proj_ANN<-subset(grid_of_points_enriched, select = -c(x, y))
    
    if(model_projection==FALSE){
      cat("Step 4: ANN training with multiple hidden layers\n")
      training_set_features_only_ANN<-subset(training_set, select = -c(x, y))
      f <- as.formula(paste("t", "~", paste(input_column_names_codes, collapse = " + ") ))
      set.seed(20)
      
      #initialise features that will record the optimal scores
      super_hidden<-0
      super_accuracy<-0
      super_accuracy_self<-0
      super_threshold<-0
      super_nn<-NA    
      
      #for each hidden layer to test: train the ANN and self-test
      for (hidden in hiddens){
        hidden<-c(hidden)
        cat("Testing hidden neurons =",hidden,"\n")
        #train the ANN with the hidden neurons
        nn <- neuralnet(f,data = training_set_features_only_ANN,hidden = hidden, threshold = thld, stepmax = stp, rep = rp, act.fct = act.fct, linear.output = FALSE, lifesign = "minimal", algorithm = alg) 
        # Compute predictions on the training data
		if (length(input_column_names_codes)>1)
          pr.nn1<- compute(nn, training_set_features_only_ANN[,1:length(input_column_names_codes)])
        else
          pr.nn1<- compute(nn, training_set_features_only_ANN)												  
        training_set_features_only_ANN$pred<- pr.nn1$net.result
        #test multiple decision thresholds based on the values over the training set
        decision_thresholds<-as.numeric(quantile(training_set_features_only_ANN$pred,probs=c(0.01,0.05,0.1,0.2,0.5,0.8))) #take 1% of the training set out
        opt_accuracy<-0
        opt_threshold<-0
        #loop on decision thresholds testing
        for (decision_threshold in decision_thresholds){
          training_set_features_only_ANN$detected<-F
          training_set_features_only_ANN$detected[which( (training_set_features_only_ANN$t ==0 & training_set_features_only_ANN$pred<decision_threshold) | 
                                                           (training_set_features_only_ANN$t ==1 & training_set_features_only_ANN$pred>decision_threshold) ) ]<-T
          accuracy<- length(which(training_set_features_only_ANN$detected))/length(training_set_features_only_ANN$detected)
          #record the best threshold and gained accuracy
          if (opt_accuracy<accuracy){
            opt_accuracy=accuracy
            opt_threshold<-decision_threshold
          }
          #cat("Accuracy selftest =", accuracy*100,"%","(thr:",decision_threshold,")","\n")
        }#end loop on decision thresholds
        opt_accuracy_self<-opt_accuracy
        cat("Optimal accuracy selftest =", opt_accuracy*100,"%","(thr:",opt_threshold,")","hidden neurons:",hidden,"\n")
        
        
        #if nfold>0 do cross validation with 95%-5% approach
        if (nfold>0){
          cat("Cross-validating..\n")
          proportion <- 0.95 # Set to 0.995 for LOOCV
          accuracies <- NULL
          #for each fold, select a random (95%) subset for training and another (5%) for testing
          for(i in 1:nfold) {
            cat(i," ")
            #random selection of training and testing rows
            index    <- sample(1:nrow(training_set_features_only_ANN), round(proportion*nrow(training_set_features_only_ANN)))
            train_cv <- training_set_features_only_ANN[index, ]
            test_cv  <- training_set_features_only_ANN[-index, ]
            #ANN training
            nn_cv    <- neuralnet(f,data = train_cv,hidden = hidden,threshold = thld,stepmax = stp,rep = rp,act.fct = act.fct,linear.output = FALSE,algorithm = alg)
			if (length(input_column_names_codes)>1)
              pr.nn1     <- compute(nn, test_cv[,1:length(input_column_names_codes)]) 
            else
              pr.nn1     <- compute(nn, test_cv)
            test_cv$pred    <- pr.nn1$net.result
            decision_thresholds<-as.numeric(quantile(test_cv$pred,probs=c(0.01,0.05,0.1,0.2,0.5,0.8))) #take 1% of the training set out
            #take the best threshold and accuracy
            opt_accuracy_cv<-0
            #loop on decision thresholds testing
            for (decision_threshold in decision_thresholds){
              test_cv$detected<-F
              test_cv$detected[which( (test_cv$t ==0 & test_cv$pred<decision_threshold) | 
                                        (test_cv$t ==1 & test_cv$pred>decision_threshold) ) ]<-T
              accuracy_cv   <- length(which(test_cv$detected))/length(test_cv$detected)
              if (opt_accuracy_cv<accuracy_cv){
                opt_accuracy_cv=accuracy_cv
              }
            }#end loop on decision thresholds
            accuracies[i]   <- opt_accuracy_cv
          }#end loop on folds  
          mean_accuracy_cv <- mean(accuracies)
          cat("\nCrossvalidation accuracy =",mean_accuracy_cv,"(min",min(accuracies),", max",max(accuracies),")\n\n")
          opt_accuracy<-mean_accuracy_cv
        }#end nfold>0
        
        #if the current nfold is the optimal among all models record it
        if (super_accuracy<opt_accuracy){
          super_accuracy<-opt_accuracy
          super_threshold<-opt_threshold #the optimal threshold will be the one on the training set
          super_accuracy_self<-opt_accuracy_self
          super_hidden<-hidden #record the optimal number of hidden layers
          super_nn<-nn
        }
      }#end loop on hidden layers
      
      cat("Results for",ID,":",
          "Optimal overall accuracy nfold =", super_accuracy*100,"%","(nfold =",nfold,")",
          "Optimal overall accuracy selftest =", super_accuracy_self*100,"%","(thr:",super_threshold,")","hidden neurons:",super_hidden,"\n")
      
      #project the optimal ANN on the grid
      pr.nn_grid<- compute(super_nn, grid_of_points_enriched_proj_ANN)
      grid_of_points$probability<-pr.nn_grid$net.result
      grid_of_points[rowSums(is.na(grid_of_points)) > 0,3]<- -9999
      
      #WRITE METADATA AND OUTPUT
      fileConn<-file(paste0(output_folder,"/ANN/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt"))
      writeLines(c(
        paste0("ID name = ",ID),
        paste0("Optimal overall accuracy nfold = ",super_accuracy),
        paste0("Native distribution = ",nativedistribution),
        paste0("Number of folders of cross validation = ",nfold),
        paste0("Optimal overall accuracy selftest = ",super_accuracy_self),
        paste0("Optimal decision threshold = ",super_threshold),
        paste0("Optimal number of neurons = ",super_hidden),
        paste0("Spatial resolution = ",resolution),
        paste0("rp = ",rp),
        paste0("thld = ",thld),
        paste0("stp = ",stp),
        paste0("alg = ",alg),
        paste0("act.fct = ",act.fct)
      ), fileConn)
      close(fileConn)
    }else{
      cat("Step 4: ANN projection session\n")
      raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
      neural_out_file_name = paste0(trained_neural_networks,"/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_ann.Rdata")
      metadata_file<-paste0(trained_neural_networks,"/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt")
      metadata_output_file<-paste0(output_folder,"/ANN/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt")
      
      load(neural_out_file_name)
      
      if(file.exists(raster_out_file_name)){
        cat("Skipping\n")
        next
      }
      #OCCURRENCE ENRICHMENT AND GRID PREPARATION
      cat("Step 1: Occurrence enrichment and grid preparation\n")
      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
      presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
      presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = presence$longitude, resolution = resolution)
      presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = presence$latitude, resolution = resolution)
      
      #project the optimal ANN on the grid
      pr.nn_grid<- compute(super_nn, grid_of_points_enriched_proj_ANN)
      grid_of_pointssp<-grid_of_points
      grid_of_pointssp$probability<-pr.nn_grid$net.result
      grid_of_pointssp[rowSums(is.na(grid_of_pointssp)) > 0,3]<- -9999
      
      #WRITE METADATA AND OUTPUT
      fileConn<-file.copy(from = metadata_file, to = metadata_output_file)
      
      
    }
    if (nativedistribution){
      cat("\nSTEP 4b: Adjusting the projection for native distribution...\n")
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
      
      if(model_projection==FALSE){ 
        
        mindistance<-sapply(1:length(grid_of_points$probability), function(i){
          p_long<-grid_of_points$x[i]
          p_lat<-grid_of_points$y[i]
          dx<-(presence$longitude_res-p_long)
          dy<-(presence$latitude_res-p_lat)
          sqrd<-(dx*dx)+(dy*dy)
          d<-min(sqrt(sqrd))
          return (d)
        },simplify = T)
      }else{
        mindistance<-sapply(1:length(grid_of_pointssp$probability), function(i){
          p_long<-grid_of_pointssp$x[i]
          p_lat<-grid_of_pointssp$y[i]
          dx<-(presence$longitude_res-p_long)
          dy<-(presence$latitude-p_lat)
          sqrd<-(dx*dx)+(dy*dy)
          d<-min(sqrt(sqrd))
          return (d)
        },simplify = T)
        
      }
      #setup probability weights
      weights<-dnorm(mindistance,0,sigma)/dnorm(0,0,sigma)
      weights[which(mindistance<=sigma)]<-1
      #weight probability by distance from presence points
      if(model_projection==FALSE){ 
        na_points<-which(grid_of_points$probability==-9999)
        grid_of_points$probability[-na_points]<-grid_of_points$probability[-na_points]*weights[-na_points]
      }else{
        na_points<-which(grid_of_pointssp$probability==-9999)
        grid_of_pointssp$probability[-na_points]<-grid_of_pointssp$probability[-na_points]*weights[-na_points]     
      }
    }   
    
    #WRITE THE ASC file
    cat("\nSTEP 5: Projecting the ANN...\n")
    #raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
    #define the grid points to populate the dataset
    if(model_projection==FALSE){  
      ypoints<-unique(grid_of_points$y)
      xpoints<-unique(grid_of_points$x)
    }else{
      ypoints<-unique(grid_of_pointssp$y)
      xpoints<-unique(grid_of_pointssp$x)   
      
    }
    ncol_r<-length(xpoints)
    nrow_r<-length(ypoints)
    #create a new raster with the same extent and resolution of the first layer
    ro <- raster(ncol=ncol_r, nrow=nrow_r)
    #length(values(ro))
    
    #res(ro) <- resolution
    #length(values(ro))
    extent(ro)<-extent(first_raster_data)
    #populate the matrix
    values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
    row_counter<-1
    for (y_c in 1:(nrow_r)){
      yp<-ypoints[y_c]
      if(model_projection==FALSE){    
        row_rast<-grid_of_points[which(grid_of_points$y == yp),]
      }else{
        row_rast<-grid_of_pointssp[which(grid_of_points$y == yp),]
      }
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
    writeRaster(ro, filename=raster_out_file_name_ANN, format="ascii",overwrite=TRUE)
    neural_out_file_name = paste0(output_folder,"/ANN/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_ann.Rdata")
    if(Create_ANN_Rdata){
      save(super_nn, file=neural_out_file_name) #this line save also the ANN to be used on another set of Environmental parameter
    }
    cat("ID",ID,"ANN done.\n")
  }   
  
  if(AQUAMAPS_Active == TRUE){
    aqprob<-function(value,f_index,AQ_quantiles_per_feature){
      quantiles<-AQ_quantiles_per_feature[[f_index]]
      
      if (as.numeric(quantiles[1])==as.numeric(quantiles[5])){
        return (1)
      } 
      q1<-as.numeric(quantiles[1])
      q2<-as.numeric(quantiles[2])
      q3<-as.numeric(quantiles[3])
      q4<-as.numeric(quantiles[4])
      q5<-as.numeric(quantiles[5])
      q1<-q1-0.2*abs(q1)
      q2<-q2-0.2*abs(q2)
      q4<-q4+0.2*abs(q4)
      q5<-q5+0.2*abs(q5)
      
      if (value<q1){
        return (0)
      }
      if (value>q5){
        return (0)
      }
      if (value<q2){
        if (q2 == q1)
          y = 1
        else
          y<-(value-q1)/(q2-q1)
      }else if (value<=q4){
        return (1)
      }else if (value==q5){
        return (0)
      }else{
        if (q4==q5)
          y =1
        else
          y<-((q4-value)/(q5-q4)) +1
      }
    }
    
    if(model_projection == TRUE){
      cat("\nStep 6: AquaMaps projection session\n")
      raster_out_file_name_AquaMaps = paste0(output_folder,"/AquaMaps/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
      aquamaps_out_file_name = paste0(trained_aquamaps,paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_aquamaps.Rdata")
      load(aquamaps_out_file_name)
      
      if(file.exists(raster_out_file_name_AquaMaps)){
        stop("Cannot overwrite the original file\n")
      }
      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
      presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
      presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = presence$longitude, resolution = resolution)
      presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = presence$latitude, resolution = resolution)
      
    }else{
      cat("\nStep 6: AquaMaps training\n")
      
      training_set_features_only_AQ<-subset(training_set, select = -c(x, y))
      training_set_features_only_AQ<-training_set_features_only_AQ[which(training_set_features_only_AQ$t==1),]
      training_set_features_only_AQ$t<-NULL
      AQ_quantiles_per_feature<-list()
      for (q in 1:ncol(training_set_features_only_AQ)){
        qq<-quantile(training_set_features_only_AQ[,q])
        AQ_quantiles_per_feature[[q]]<-qq
      }
    }    #end else training/projection 
    probAQ<-list()
    listcount<-1
    grid_of_points_enriched_features_only_AQ<-subset(grid_of_points_enriched, select = -c(x, y),)
    for (i in 1:nrow(grid_of_points_enriched_features_only_AQ)){
      
      featuresrow<-grid_of_points_enriched_features_only_AQ[i,]
	  if (length(featuresrow)==1)
        featuresrow<-data.frame(f1=featuresrow)									   
      if (length(which(is.na(featuresrow)))>0){
        probAQ[[listcount]]<--9999
      }else{
        probs<-sapply(1:ncol(featuresrow), function(j){
          aqp<-aqprob(as.numeric(featuresrow[j]),j,AQ_quantiles_per_feature)
        },simplify = T)
        probAQ_f<-prod(probs)
        probAQ[[listcount]]<-probAQ_f
      }
      
      listcount<-listcount+1
    }
    
    prob_arr_aq<-array( unlist( probAQ ))
    grid_of_points_AQ<-grid_of_points
    grid_of_points_AQ$probability<-prob_arr_aq
    
    
    #WRITE METADATA AQ
    fileConn<-file(paste0(output_folder,"/AquaMaps/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt"))
    writeLines(c(
      paste0("ID name = ",ID),
      paste0("Spatial resolution = ",resolution),
      paste0("Optimal decision threshold = 0.25")
    ), fileConn)
    close(fileConn)  
    
    # Native Distribution
    if (nativedistribution){
      cat("\nSTEP 6b: Adjusting AquaMaps for native distribution...\n")
      presence_aq<-
        minx_aqm<-min(presence$longitude_res)-0.2*abs(min(presence$longitude_res))
      maxx_aqm<-max(presence$longitude_res)+0.2*abs(max(presence$longitude_res))
      miny_aqm<-min(presence$latitude_res)-0.2*abs(min(presence$latitude_res))
      maxy_aqm<-max(presence$latitude_res)+0.2*abs(max(presence$latitude_res))
      grid_of_points_AQ$probability[which((grid_of_points_AQ$x<minx_aqm) | (grid_of_points_AQ$x>maxx_aqm))]<--9999
      grid_of_points_AQ$probability[which((grid_of_points_AQ$y<miny_aqm) | (grid_of_points_AQ$y>maxy_aqm))]<--9999
    }
    
    #WRITE THE ASC file
    cat("\nSTEP 7: Projecting AquaMaps...\n")
    ypoints<-unique(grid_of_points$y)
    xpoints<-unique(grid_of_points$x)
    ncol_r<-length(xpoints)
    nrow_r<-length(ypoints)
    #create a new raster with the same extent and resolution of the first layer
    ro <- raster(ncol=ncol_r, nrow=nrow_r)
    #length(values(ro))
    #res(ro) <- resolution
    #length(values(ro))
    extent(ro)<-extent(first_raster_data)
    #populate the matrix
    values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
    row_counter<-1
    for (y_c in 1:(nrow_r)){
      yp<-ypoints[y_c]
      row_rast<-grid_of_points_AQ[which(grid_of_points_AQ$y == yp),]
      row_rast<-row_rast[order(row_rast$x),]
      values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
      row_counter<-row_counter+1
      
    }
    values_vec<-as.vector(t(values))
    values_vec[which(values_vec==0)]<--9999
    values(ro)<-values_vec
    NAvalue(ro)<- -9999
    #save the raster
    cat("Writing the output..\n")
    writeRaster(ro, filename=raster_out_file_name_AquaMaps, format="ascii",overwrite=TRUE)
    aquamaps_out_file_name = paste0(output_folder,"/AquaMaps/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_aquamaps.Rdata")
    if(Create_AquaMaps_Rdata){
      save(AQ_quantiles_per_feature, file=aquamaps_out_file_name)
    }
    cat("ID",ID," AquaMaps done.\n")
  } #end AquaMaps
  
  if(MAXENT_Active == TRUE){
    meprob<-function(maxentmodel,grid_of_points,ID){
      max_out<-paste0(paste0(maxentmodel,gsub(" ","_",ID)),".asc")
      asc_file<-raster(max_out)
      grid_of_points$probability<-NULL
      grid_of_points_extracted_me<-extract(x=asc_file,y=grid_of_points,method='simple')
      grid_of_points_me<-grid_of_points
      grid_of_points_me$probability<-grid_of_points_extracted_me
      return(grid_of_points_me)
    }
    
    get_thr_me<-function(maxentmodel,ID){
      max_out<-paste0(paste0(maxentmodel,gsub(" ","_",ID)),".html")
      fileBMPME<-readChar(max_out, file.info(max_out)$size)
      METS<-gsub(".*>Training omission rate</th><tr align=center><td>","",fileBMPME)
      METS<- gsub("</td><td>Fixed cumulative value 1.*","",METS)    
      METS<- as.double(gsub(".*</td><td>","",METS))
      return(METS)
    }
    
    metraining<-function(ID,env_data_folder,output_folder,prevalence){
      cat("...training MaxEnt with prevalence =",maxent_prevalence,"...\n")

      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
      presence_data<-paste0(presence_data_folder,"/",presence_file)
      maxent_out<-paste0(output_folder,"/MaxEnt/ME_",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "/")
      
      if(!dir.exists(maxent_out)) dir.create(maxent_out)
      
      command<-paste0("java -jar ./max_ent_cyb.jar \"",presence_data,"\" \"",env_data_folder, "/\" \"",maxent_out,"\" ",prevalence)
      
      maxent_execution<-system(command, intern = T,
                               ignore.stdout = FALSE, ignore.stderr = FALSE,
                               wait = TRUE, input = NULL, show.output.on.console = TRUE,
                               minimized = FALSE, invisible = TRUE)
      
      execution_success<-(length(which(grepl(pattern="OK MaxEnt",x=maxent_execution)))>0)
      cat("MaxEnt training OK=",execution_success,"\n")
      if (!execution_success){
        print(maxent_execution)
        stop("MaxEnt was not successful")
      }
      return(maxent_out)
    }
    
    meprojecting<-function(maxentmodel,ID,env_data_folder,output_folder,prevalence){
      
      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
      presence_data<-paste0(presence_data_folder,"/",presence_file)
      maxent_out_prj<-paste0(output_folder,"/MaxEnt/ME_",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_prj/")
      lambda_file<-paste0(paste0(maxentmodel,gsub(" ","_",ID)),".lambdas")
      if(!dir.exists(maxent_out_prj)) 
        dir.create(maxent_out_prj)
      
      maxent_out_prj_asc<-paste0(paste0(maxent_out_prj,gsub(" ","_",ID)),".asc")
      
      #java -cp max_ent_cyb.jar org.gcube.datanalysis.ecomod.MainProjecting ./Presence_Abudefduf_saxatilis.csv ./env_pars/ ./out_programatic/ ./out_programatic/Abudefduf_saxatilis.lambdas ./env_pars_2050/ ./test.asc
      
      command<-paste0("java -cp max_ent_cyb.jar org.gcube.datanalysis.ecomod.MainProjecting \"",
                      presence_data,"\" \"",env_data_folder, "/\" \"",maxentmodel,"\" \"",lambda_file,"\" \"",env_data_folder,"\" \"",maxent_out_prj_asc,"\" ",prevalence)
      
      maxent_execution<-system(command, intern = T,
                               ignore.stdout = FALSE, ignore.stderr = FALSE,
                               wait = TRUE, input = NULL, show.output.on.console = TRUE,
                               minimized = FALSE, invisible = TRUE)
      
      execution_success<-(length(which(grepl(pattern="OK MaxEnt",x=maxent_execution)))>0)
      cat("MaxEnt projection OK=",execution_success,"\n")
      if (!execution_success){
        print(maxent_execution)
        stop("MaxEnt was not successful")
      }
      
      max_out_previous<-paste0(paste0(maxentmodel,gsub(" ","_",ID)),".html")
      max_out_new<-paste0(paste0(maxent_out_prj,gsub(" ","_",ID)),".html")
      file.copy(max_out_previous, max_out_new)
      
      max_out_previous<-paste0(paste0(maxentmodel,gsub(" ","_",ID)),".lambdas")
      max_out_new<-paste0(paste0(maxent_out_prj,gsub(" ","_",ID)),".lambdas")
      file.copy(max_out_previous, max_out_new)
      
      return(maxent_out_prj)
      
    }
    
    if(model_projection == TRUE){
      cat("\nStep 8: MaxEnt projection session\n")
      raster_out_file_name_MaxEnt = paste0(output_folder,"/MaxEnt/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , ".asc")
      maxent_out_file_name = paste0(trained_maxent,paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_maxent.Rdata")
      load(maxent_out_file_name)
      
      presence_file<-paste0(paste("Presence",gsub(" ","_",ID),sep = "_"),".csv")
      presence<-read.table(paste0(presence_data_folder,"/",presence_file),header = TRUE, sep=",")
      presence$longitude_res<-coordinate_at_res(origin = min_x_in_raster,coordinate = presence$longitude, resolution = resolution)
      presence$latitude_res<-coordinate_at_res(origin = min_y_in_raster,coordinate = presence$latitude, resolution = resolution)
      
      if(file.exists(raster_out_file_name_MaxEnt)){
        stop("Cannot overwrite the original file\n")
      }
      maxentmodel<-meprojecting(maxentmodel,ID,projection_environmental_layers,output_folder,maxent_prevalence)
    }else{
      cat("\nStep 8: MaxEnt training\n")
      maxentmodel<-metraining(ID,env_data_folder,output_folder,maxent_prevalence)
    }#end else training/projection 
    
    probME<-list()
    listcount<-1
    grid_of_points_ME<-meprob(maxentmodel,grid_of_points,ID)
    decision_treshold<-get_thr_me(maxentmodel,ID)
    
    #WRITE METADATA ME
    
    fileConn<-file(paste0(output_folder,"/MaxEnt/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_metadata.txt"))
    writeLines(c(
      paste0("ID name = ",ID),
      paste0("Spatial resolution = ",resolution),
      #Fixed cumulative value 1.
      paste0("Optimal decision threshold = ",decision_treshold)
    ), fileConn)
    close(fileConn)  
    
    # Native Distribution
    if (nativedistribution){
      cat("\nSTEP 8b: Adjusting the projection for native distribution...\n")
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
      mindistance<-sapply(1:length(grid_of_points_ME$probability), function(i){
        p_long<-grid_of_points_ME$x[i]
        p_lat<-grid_of_points_ME$y[i]
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
      na_points<-which( (grid_of_points_ME$probability==-9999) | (is.na(grid_of_points_ME$probability)) )
      grid_of_points_ME$probability[-na_points]<-grid_of_points_ME$probability[-na_points]*weights[-na_points]     
    }
    
    #WRITE THE ASC file
    cat("\nSTEP 9: Projecting MaxEnt...\n")
    ypoints<-unique(grid_of_points$y)
    xpoints<-unique(grid_of_points$x)
    ncol_r<-length(xpoints)
    nrow_r<-length(ypoints)
    #create a new raster with the same extent and resolution of the first layer
    ro <- raster(ncol=ncol_r, nrow=nrow_r)
    #length(values(ro))
    #res(ro) <- resolution
    #length(values(ro))
    extent(ro)<-extent(first_raster_data)
    #populate the matrix
    values<-matrix(nrow = nrow_r,ncol = ncol_r,data = -9999)
    row_counter<-1
    for (y_c in 1:(nrow_r)){
      yp<-ypoints[y_c]
      row_rast<-grid_of_points_ME[which(grid_of_points_ME$y == yp),]
      row_rast<-row_rast[order(row_rast$x),]
      values[(nrow_r-row_counter+1),]<-row_rast$probability[1:(ncol_r)]
      row_counter<-row_counter+1
    }
    values_vec<-as.vector(t(values))
    values_vec[which(values_vec==0)]<--9999
    values(ro)<-values_vec
    NAvalue(ro)<- -9999
    #save the raster
    cat("Writing the output..\n")
    writeRaster(ro, filename=raster_out_file_name_MaxEnt, format="ascii",overwrite=TRUE)
    maxent_out_file_name = paste0(output_folder,"/MaxEnt/",paste0(gsub(pattern = " ",replacement = "_",x = ID)) , "_maxent.Rdata")
    if(Create_MAXENT_Rdata){
      save(maxentmodel, file=maxent_out_file_name)
    }
    cat("ID",ID," MaxEnt done.\n")
  } #end MaxEnt
  
  
}