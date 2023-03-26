rm(list=ls(all=TRUE))
library(dplyr) 
library(readxl)
library(raster)
library(properties)

Workflow_Parameters_File<- paste0("./Workflow_Parameters.txt")
WPFprops <- read.properties(Workflow_Parameters_File)

cat("\nME data format conversion\n")

main_folder<-"."
dir.create(file.path(main_folder, "ENM"), showWarnings = FALSE)
output_folder<-"./ENM"
ME_input_folder<-"./MaxEnt"
ME_ouput_folder<-paste0(output_folder,"/ME")
if(!dir.exists(ME_ouput_folder)) dir.create(ME_ouput_folder)

all_IDs<-list.files(ME_input_folder)


for (ID in all_IDs){
  
  if(grepl( ".asc", ID, fixed = TRUE)){
    print(ID)
    file.copy(from=paste0(ME_input_folder,"/",ID), to=paste0(ME_ouput_folder,"/",gsub("_native","",ID)), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
    #ro <- raster(paste0(ME_input_folder,"/",ID))
    #resolution <- NA
    #resolution <- unique(res(ro))
  }
  if(grepl( ".html", ID, fixed = TRUE)){
    BMPFileME<- paste0(ME_input_folder,"/",ID)
    fileBMPME<-readChar(BMPFileME, file.info(BMPFileME)$size)
    METS<-gsub(".*>Training omission rate</th><tr align=center><td>","",fileBMPME)
    METS<- gsub("</td><td>Fixed cumulative value 1.*","",METS)    
    METS<- as.double(gsub(".*</td><td>","",METS))
    print(METS)
    
    presence_record<-gsub(".*the run:<br>","",fileBMPME)
    presence_record<- as.double(gsub("presence records used for training..*","",presence_record))
    
    background_point<-gsub(".*training.<br>","",fileBMPME)
    background_point<- as.double(gsub("points used to determine the Maxent distribution .*","",background_point))
    
    training_AUC<-gsub(".*training AUC is ","",fileBMPME)
    training_AUC<- as.double(gsub(", unregularized training gain .*","",training_AUC))  
    
    
    ro <- raster(paste0(ME_input_folder,"/",unique(gsub(".html","_native.asc", ID))))
    resolution <- unique(res(ro))
    
    
    
    #WRITE METADATA
    fileConn<-file(paste0(ME_ouput_folder,"/",paste0(gsub(pattern = ".html",replacement = "",x = ID)) , "_metadata.txt"))
    writeLines(c(
      paste0("ID name = ",gsub(pattern = ".html",replacement = "",x = ID)),
      paste0("Optimal decision threshold = ",METS),
      paste0("Presence records used for training = ",presence_record),
      paste0("Points used to determine the Maxent distribution = ",background_point),
      paste0("Training AUC = ",training_AUC),
      paste0("Spatial resolution = ",resolution)
    ), fileConn)
    close(fileConn)
    
    }

  }
