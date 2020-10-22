setwd("/Users/vinitaperiwal/GrowthCurver/")
library(dplyr)
library(tibble)
library(reshape2)
library(growthcurver)

# source the script with functions
source("/Users/vinitaperiwal/GrowthCurver/Scripts/Functions.R")

#read all files with .tab extension
file.names <- list.files(path = '/Users/vinitaperiwal/GrowthCurver/', recursive = TRUE, pattern = "\\.tab$") #recursive reads through all subfolders and files
file.names

#loop for each .tab file, fits logistic curve, annotates
for(i in 1:length(file.names)){
  
  #read file/table
  input_file<-tools::file_path_sans_ext(file.names[i]) #read filename w/o extension
  input_file
  
  #create a new dataframe with additional annotations for each file
  
  annot <- tibble(Bug_ID = create_annot(input_file)[1], Replicate_no = create_annot(input_file)[2], Plate_no = create_annot(input_file)[3],
                  Drug_name = create_annot(input_file)[5])
  
  X<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  
  head(X)
  colnames(X)[1]<-"time"
  
  annot_file<-merge(annot,X)
  head(annot_file)
  
  #write annotated files
  annot_outfile<-paste0(input_file,".annot")
  write.table(annot_file, file = annot_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

#merge all replicates of a bug
# Merge files across all replicates of a strain (.annot files)
dir.names<-list.dirs(path = '/Users/vinitaperiwal/GrowthCurver', recursive = FALSE)
dir.names

for(d in 1:length(dir.names)){
  
  if(str_detect(dir.names[d], "/Users/vinitaperiwal/GrowthCurver/NT5")){
    
    all_rep<-do.call(rbind, lapply(a<-list.files(path=dir.names[d], pattern="\\.annot$", full.names = TRUE), function(i){
      read.table(i,header=TRUE, sep="\t", stringsAsFactors = FALSE)}))
    write.table(all_rep, file = paste0(dir.names[d],"_all_replicates.merged"), sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
  
}

##################################################################################################################################
##################################################################################################################################
# below code will use growthcurver to fit curves
# to compare replicates all curves should have equivalent time readings, time trim replicates to have same time points

bug_min_reads<-data.frame(Bug_ID = character(), Replicate_no = character(), count_tp = double(), stringsAsFactors = FALSE)

#create an empty data frame to store all raw points and fitted points
fitted_values<-data.frame(stringsAsFactors = FALSE)

#create an empty data frame to store all computed model parameters
gc_fit_params<-data.frame(stringsAsFactors = FALSE)

#read all files with .tab extension
file.merged<-list.files(path = '/Users/vinitaperiwal/GrowthCurver/', pattern = ".merged$") #recursive reads through all subfolders and files
file.merged

#loop for each .tab file, fits logistic curve, annotates
for(j in 1:length(file.merged)){
  
  #read file/table
  in_file<-tools::file_path_sans_ext(file.merged[j]) #read filename w/o extension
  in_file
  
  Y<-data.frame(read.table(file.merged[j], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  head(Y)
  
  #count time points in all replicates, all should have same time point reading else trim if needed
  
  time_points<-Y %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
    summarise(count_tp = max(unique(time)))
  head(time_points)
  time_points<-time_points %>% mutate(min_tp = min(count_tp))
  min_tp<-min(time_points$count_tp)
  min_tp
  
  bug_min_reads<-rbind(bug_min_reads, data.frame(time_points))
  bug_min_reads
  
  Z<-Y %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Drug_name) %>%
    filter(time <= min_tp)
  
  head(Z)
  nrow(Z)
  
  write.table(Z, file = paste0("/Users/vinitaperiwal/GrowthCurver/trim_points/",in_file,".trim"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  unique_rep<-unique(Z$Replicate_no)
  unique_rep
  
  unique_plate<-unique(Z$Plate_no)
  unique_plate
  
  #loop plate wise and then well wise (col_name) to fit growth curve on each well, later merge all data in a data frame
  for(j in 1:length(unique_plate)){
    
    dat<-Z %>% filter(Plate_no == unique_plate[j])
    dat
    
    if (dim(dat)[1] != 0){ #if missing plate
      
      for(k in 1:length(unique_rep)){
        
        plate_dat<-dat %>% filter(Replicate_no == unique_rep[k])
        plate_dat
        
        if (dim(plate_dat)[1] != 0){ #if missing replicate eg NT5022, plate16 has no rep1
          
          bug<-plate_dat[1,1]
          bug
          rep<-plate_dat[1,2]
          rep
          plate<-plate_dat[1,3]
          plate
          drug<-plate_dat[1,4]
          drug
          
          #plot_file<-paste0(ID,"_",rep,"_",plate,"_",drug)
          #plot_file
          
          data<-plate_dat[,-(1:4)]
          data
          
          #fit model to each well of a plate and save the fitted model
          for(col_name in names(data)){
            
            col_name
            
            if(col_name != "time"){
              
              current_well<-data[, c("time",col_name)]
              current_well
              min_OD<-min(current_well[, col_name])
              min_OD
              
              #do background correction using min OD of each well
              current_well[,col_name]<-current_well[,col_name] - min_OD
              current_well[,col_name]
              
              #each time create a new variable for each well
              gc_fit<-SummarizeGrowth(data_t = current_well[,"time"], data_n = current_well[,col_name])
              #saveRDS(gc_fit, file = paste0("/Users/vinitaperiwal/GrowthCurver/models/",ID,"_",rep,"_",plate,"_",drug,"_",col_name,".rds"))
              
              gc_fit
              
              #create a data frame of raw values and fitted values
              mod_t<-data.frame(matrix(unlist(gc_fit$data$t)))
              mod_t
              mod_N<-data.frame(matrix(unlist(gc_fit$data$N)))
              mod_N
              
              if(gc_fit$vals$k != 0 & gc_fit$vals$n0 != 0 & gc_fit$vals$r != 0){
                
                mod<-cbind(mod_t, mod_N, gc_fit$model$m$fitted()) #m is the model object with all fitted and residual values
                
              }else{
                
                mod<-cbind(mod_t, mod_N, gc_fit$vals$r)
                
              }
              
              colnames(mod)<-c("time","OD","fitted")
              mod
              
              #add annotation to each row
              annot_mod<-cbind(bug,rep,plate,drug,col_name,mod)
              annot_mod
              
              fitted_values<-rbind(fitted_values, annot_mod)
              
              fitted_param<-as.data.frame(cbind(bug,rep,plate,drug,col_name,gc_fit$vals$k,gc_fit$vals$k_se,gc_fit$vals$k_p,gc_fit$vals$n0,gc_fit$vals$n0_se,gc_fit$vals$n0_p,gc_fit$vals$r,gc_fit$vals$r_se,gc_fit$vals$r_p,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$t_gen,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note))
              colnames(fitted_param)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")
              fitted_param
              
              gc_fit_params<-rbind(gc_fit_params, fitted_param)
              tail(gc_fit_params)
              
            }
          }
          
        }
        
      }
    }
    
  }
  
}


head(fitted_values)
nrow(fitted_values) #1,292,544

View(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 53,184

write.table(fitted_values, file = "/Users/vinitaperiwal/GrowthCurver/GC_fitted_model_values", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gc_fit_params, file = "/Users/vinitaperiwal/GrowthCurver/GC_fit_params", sep = "\t", quote = FALSE, row.names = FALSE)
