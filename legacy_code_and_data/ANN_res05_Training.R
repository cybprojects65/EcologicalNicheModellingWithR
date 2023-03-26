rm(list=ls(all=TRUE))
library(dplyr)
library(raster)
library(neuralnet)

#GENERAL PARAMETERS
env_data_folder<-"/media/gp/228E389F8E386E03/All ENMs EcoScope/Environmental Parameters/Historical 2019"
presence_data_folder<-"/media/gp/228E389F8E386E03/All ENMs EcoScope/Presence Records/Non Fish"
absence_data_folder<-"/media/gp/228E389F8E386E03/All ENMs EcoScope/Absence Records/Non Fish"
output_folder<-"./ENM_2019_NF"

allspp<-list.files(presence_data_folder)
allspp<-gsub("Presence_","",allspp)
allspp<-gsub(".csv","",allspp)
allspp<-gsub("_"," ",allspp)

for (species in allspp){
  
  cat("\nProcessing species ",species,"\n")
  raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , ".asc")
  neural_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_ann.Rdata")
  if(file.exists(raster_out_file_name)){
    cat("Skipping\n")
    next
  }
  
  ths_sample<-100
  abs_sample<-10000
  resolution<-0.5 #deg
  rp            <- 10      # number of repetitions for the training - #20 gives more stability to the network and more independence of the initial training step
  thld          <- 0.001   # threshold for minimum decrease in overall error, default 0.01 = 1%
  stp           <- 1e+06  # the maximum steps for the training of the neural network, default 1e+05
  alg           <- "rprop-"  # possible: backprop, rprop+, rprop-, sag, or slr
  act.fct       <- "logistic"#"logistic" # possible: "logistic" or "tanh"; linear.output must be 
  nfold         <-0#20
  #hiddens<-c(10,25,50,75,100,150,200,250)
  nativedistribution<-T
  hiddens<-c(50,70,100,150,200)
  
  
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
  
  
  
  
  cat("Step 2: ANN training with multiple hidden layers\n")
  training_set_features_only<-subset(training_set, select = -c(x, y))
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
    nn <- neuralnet(f,data = training_set_features_only,hidden = hidden, threshold = thld, stepmax = stp, rep = rp, act.fct = act.fct, linear.output = FALSE, lifesign = "minimal", algorithm = alg) 
    # Compute predictions on the training data
    pr.nn1<- compute(nn, training_set_features_only[,1:length(input_column_names_codes)]) 
    training_set_features_only$pred<- pr.nn1$net.result
    #test multiple decision thresholds based on the values over the training set
    decision_thresholds<-as.numeric(quantile(training_set_features_only$pred,probs=c(0.01,0.05,0.1,0.2,0.5,0.8))) #take 1% of the training set out
    opt_accuracy<-0
    opt_threshold<-0
    #loop on decision thresholds testing
    for (decision_threshold in decision_thresholds){
      training_set_features_only$detected<-F
      training_set_features_only$detected[which( (training_set_features_only$t ==0 & training_set_features_only$pred<decision_threshold) | 
                                                   (training_set_features_only$t ==1 & training_set_features_only$pred>decision_threshold) ) ]<-T
      accuracy<- length(which(training_set_features_only$detected))/length(training_set_features_only$detected)
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
        index    <- sample(1:nrow(training_set_features_only), round(proportion*nrow(training_set_features_only)))
        train_cv <- training_set_features_only[index, ]
        test_cv  <- training_set_features_only[-index, ]
        #ANN training
        nn_cv    <- neuralnet(f,data = train_cv,hidden = hidden,threshold = thld,stepmax = stp,rep = rp,act.fct = act.fct,linear.output = FALSE,algorithm = alg)
        pr.nn1     <- compute(nn, test_cv[,1:length(input_column_names_codes)]) 
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
  
  cat("Results for",species,":",
      "Optimal overall accuracy nfold =", super_accuracy*100,"%","(nfold =",nfold,")",
      "Optimal overall accuracy selftest =", super_accuracy_self*100,"%","(thr:",super_threshold,")","hidden neurons:",super_hidden,"\n")
  
  #project the optimal ANN on the grid
  grid_of_points_enriched_proj<-subset(grid_of_points_enriched, select = -c(x, y))
  pr.nn_grid<- compute(super_nn, grid_of_points_enriched_proj)
  grid_of_points$probability<-pr.nn_grid$net.result
  grid_of_points[rowSums(is.na(grid_of_points)) > 0,3]<- -9999
  
  
  #WRITE METADATA AND OUTPUT
  fileConn<-file(paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , "_metadata.txt"))
  writeLines(c(
    paste0("Species name = ",species),
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
    mindistance<-sapply(1:length(grid_of_points$probability), function(i){
      p_long<-grid_of_points$x[i]
      p_lat<-grid_of_points$y[i]
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
    na_points<-which(grid_of_points$probability==-9999)
    grid_of_points$probability[-na_points]<-grid_of_points$probability[-na_points]*weights[-na_points]
  }
  
  #WRITE THE ASC file
  cat("\nSTEP 3: Projecting the ANN...\n")
  #raster_out_file_name = paste0(output_folder,"/",paste0(gsub(pattern = " ",replacement = "_",x = species)) , ".asc")
  #define the grid points to populate the dataset
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
    row_rast<-grid_of_points[which(grid_of_points$y == yp),]
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
  save(super_nn, file=neural_out_file_name)
  cat("Species",species,"done.\n")
  #break
}