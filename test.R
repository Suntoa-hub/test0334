library(logging)
logReset()
basicConfig(level='FINEST')
if(file.exists("./3mcor_error.log")){
  unlink("./3mcor_error.log", recursive = FALSE, force = FALSE)
}
addHandler(writeToFile, file="./3mcor_error.log", level='DEBUG')

result = tryCatch({
  # print("prt")
  # write.csv("task_path","task1.csv")
  source('../r_api/metabolitesPrepro.R')
  source('../r_api/microbesPrepro.R')
  
  library ("WGCNA") ### -	Clustering software. Previously reported work done using v1.34
  library ("psych") ### Clustering software
  library ("plyr") ### Data transformations
  library("phyloseq")
  # install.packages(c("FactoMineR", "factoextra")) 
  # library("FactoMineR")
  
  args <- commandArgs(T)
  meta_name = args[8]
  micr_name = args[9]
  task_path = args[10]
  metaImpu=args[11]
  micrImpu=args[12]
  ifphonetype=args[13]
  phonetypename=args[14]
  
  if (args[1] == "T") {
    args[1] = TRUE;
    # dir.create("11results")
    # system.out(args[1])
  }
  if (args[2] == "T") {
    args[2] = TRUE;
  }
  if (args[3] == "T") {
    args[3] = TRUE;
  }
  if (args[4] == "T") {
    args[4] = TRUE;
  }
  if (args[5] == "T") {
    args[5] = TRUE;
  }
  if (args[6] == "T") {
    args[6] = TRUE;
  }
  if (args[7] == "T") {
    args[7] = TRUE;
  }
  
  
  if (args[1] == "NULL") {
    args[1] = FALSE;
    # dir.create("11results")
    # system.out(args[1])
  }
  if (args[2] == "NULL") {
    args[2] = FALSE;
  }
  if (args[3] == "NULL") {
    args[3] = FALSE;
  }
  if (args[4] == "NULL") {
    args[4] = FALSE;
  }
  if (args[5] == "NULL") {
    args[5] = FALSE;
  }
  if (args[6] == "NULL") {
    args[6] = FALSE;
  }
  if (args[7] == "NULL") {
    args[7] = FALSE;
  }
  # write.csv(task_path,"task.csv")
  setwd(task_path) 
  
  if(phonetypename=="NULL"){
    phenotypes<-NA
  }else{
    if(ifphonetype!="no"){
      phone_input_path = paste0("raw/",phonetypename)
      phenotypes <- read.csv (phone_input_path, row.names = 1,header = T, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
      write.csv(phenotypes,"results/raw/phenotypes.csv")
      # print(args[16])
      #args[16]=colnames(phenotypes_m)[1]
      # print(args[16])
      
    }else{
      phenotypes<-NA
    }
    
  }

  if(meta_name != "F"){
    meta_input_path = paste0("raw/",meta_name)
    meta_input<-read.csv(meta_input_path,row.names = 1, header = T, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
    
    if(micr_name != "F"){
      micro_input_path = paste0("raw/",micr_name)
      micro_input<-read.csv(micro_input_path,row.names = 1, header = T, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
    }
    # setwd(task_path) 
    dir.create("results")
    dir.create("results/raw",recursive = T)
    #metabolitesPrepro(meta_input,args[1],args[2],args[3])
    metabolitesPrepro(meta_input, missPro = args[1], missMethod = metaImpu, scalingPro = args[2], transPro = args[3], phenoData = phenotypes,phenoDataType = ifphonetype)
    
    
    
    write.csv(meta_input,"results/raw/metabolome_raw.csv")
    if(micr_name != "F"){
      #microbesPrepro(micro_input,args[4],args[5],args[6],args[7])
      microbesPrepro (micro_input,missPro = args[5],missMethod = micrImpu,rarePro =args[4], scalingPro = args[6], transPro =args[7], kValue = 3,phenoData = phenotypes,phenoDataType =ifphonetype)
      
      write.csv(micro_input,"results/raw/microbiome_raw.csv")
    }
    
  }else{
    
    micro_input_path = paste0("raw/",micr_name)
    micro_input<-read.csv(micro_input_path,row.names = 1, header = T, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
    # setwd(task_path) 
    dir.create("results")
    dir.create("results/raw",recursive = T)
    microbesPrepro (micro_input,missPro = args[5],missMethod = micrImpu,rarePro =args[4], scalingPro = args[6], transPro =args[7], kValue = 3,phenoData = phenotypes,phenoDataType =ifphonetype)
    write.csv(micro_input,"results/raw/microbiome_raw.csv")
  }
  # logging::loginfo(logger='Run successfully!')
}, error = function(e) {
  # 出现error的处理逻辑
  # error-handler-code
  logging::logerror('Preprocessing failed!', logger='')
}
)











