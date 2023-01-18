library("dplyr") 
library("readxl")
ME_input_folder<-"./MaxEnt"
ME_ouput_folder<-"./MaxEntOut"


allspp<-list.files(ME_input_folder)

#allspp<-gsub(".csv","",allspp)
#allspp<-gsub("_"," ",allspp)

for (species in allspp){
  if(grepl( ".asc", species, fixed = TRUE)){
    print(species)
    file.copy(from=paste0(ME_input_folder,"/",species), to=paste0(ME_ouput_folder,"/",gsub("_native","",species)), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
  }
  if(grepl( ".html", species, fixed = TRUE)){
    BMPFileME<- paste0(ME_input_folder,"/",species)
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
    
    #WRITE METADATA
    fileConn<-file(paste0(ME_ouput_folder,"/",paste0(gsub(pattern = ".html",replacement = "",x = species)) , "_metadata.txt"))
    writeLines(c(
      paste0("Species name = ",gsub(pattern = ".html",replacement = "",x = species)),
      paste0("Optimal decision threshold = ",METS),
      paste0("Presence records used for training = ",presence_record),
      paste0("Points used to determine the Maxent distribution = ",background_point),
      paste0("Training AUC = ",training_AUC)
    ), fileConn)
    close(fileConn)
    
    }



   
    
  }




