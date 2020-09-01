setwd("/Users/vinitaperiwal/GrowthCurver/")
library(dplyr)
library(tibble)
library(reshape2)
library(janitor) # removes empty rows and columns
library(ggfortify)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(lattice)
library(gridExtra)
library(ggsci) #fancy color palettes
library(Cairo) #high resolution images
library(ggforce) #facet wrapping
library(ggpubr) #reg equation
library(metRology) # t distribution
library(cluster) #autoplot cluster
library(broom)
library(car) #for Boxplots with labels
library(compare)
library(corrplot)
library(Hmisc)
library(pracma)
library(dunn.test)
library(purrr)
library(rcompanion) #cldList
library(rstatix) #KW test and others
library(corrr)
library(psych)
library(pROC) #auc function in growthcurver
library(growthcurver)
library(MASS) #fitdistr
library(gtools) #mixedsort
library(qvalue) #estimates FDR from a list of input p-values
library(fitdistrplus)
library(data.table)
library(multtest) #needed for metap
library(metap) #fisher test sumlog

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


###################### Plotting

#read fitted values of models
fits<-data.frame(read.table("/Users/vinitaperiwal/GrowthCurver/GC_fitted_model_values", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(fits)

fits<-fits %>% filter(Plate_no != 'plate16')
nrow(fits) #1,278,336

#read parameters
params<-data.frame(read.table("/Users/vinitaperiwal/GrowthCurver/GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
View(params)

params<-params %>% filter(Plate_no != 'plate16')
nrow(params) #52,608

# #combine all_AUC with dose information of primary screen
# doses_primary_screen<-data.frame(read.table("/Users/vinitaperiwal/GrowthCurver/Scripts/Doses_primary_screen", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
# head(doses_primary_screen)
# nrow(doses_primary_screen) #599
# 
# params_doses<-merge(params, doses_primary_screen, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name"))
# head(params_doses)
# nrow(params_doses) #54,912
# 
# write.table(params_doses, file = "/Users/vinitaperiwal/GrowthCurver/params_doses", sep = "\t", quote = FALSE, row.names = FALSE)


################################################### 
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #23

params<-data.frame(read.table("/Users/vinitaperiwal/GrowthCurver/GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))

params<-params %>% filter(Plate_no != 'plate16' & !(Bug_ID %in% c('NT5003','NT5048')))
#View(params)
nrow(params) #47,232

#select specific columns
params<-params[,c(1:5,19)]
head(params)

#removed varied dose plates from SFall
nt5004_rep1_5fu<-params %>% filter(Bug_ID == 'NT5004' & Plate_no %in% c('plate10','plate11') & Replicate_no == 'Replicate1')
nrow(nt5004_rep1_5fu)

SFall<-anti_join(params, nt5004_rep1_5fu)
nrow(SFall) #47,136 (96 wells of plate10)

nt5009_rep1_lox_5fu<-SFall %>% filter(Bug_ID == 'NT5009' & Plate_no %in% c('plate4','plate5','plate10','plate11') & Replicate_no == 'Replicate1')
nrow(nt5009_rep1_lox_5fu)

SFall<-anti_join(SFall, nt5009_rep1_lox_5fu)
nrow(SFall) #46752 (384 wells of 4 plates)

SFall<-merge(SFall,bug_desc, by="Bug_ID")
head(SFall)
nrow(SFall) #46,752

write.table(SFall, file = "/Users/vinitaperiwal/GrowthCurver/SFall", sep = "\t", quote = FALSE, row.names = FALSE)


###########################
SFall<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/SFall", header = TRUE, sep = '\t', stringsAsFactors = FALSE)

#normalizing by plate median only for food plates
normAUCl_food<-SFall %>% filter(Drug_name == 'control') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
  mutate(normAUC = round((auc_l/median(auc_l)), digits = 3))

View(normAUCl_food)
nrow(normAUCl_food) #7776

#plot median normalized AUCs and SFs
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/control_wells_normAUC.svg", sep = ""), width = 5, height = 3, bg = "white")
normAUCl_food %>% dplyr::group_by(Bug_ID,Replicate_no) %>% filter(well == 'H12' & Plate_no %in% c('plate12','plate13')) %>%
  ggplot(aes(x=Sp_short,y=normAUC)) + geom_bar(aes(fill = Replicate_no), stat = "identity", position = "dodge2") +
  th + scale_fill_jama() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Defining reference wells: from food plates
food_ref<-normAUCl_food %>% filter(normAUC >= 0.7 & well != 'H12')
head(food_ref)
nrow(food_ref) #7562

### pooled replicates
ref_allreps<-food_ref %>% dplyr::group_by(Bug_ID,Plate_no,well) %>% mutate(count = length(well)) 
head(ref_allreps)
nrow(ref_allreps) #7,562
write.table(ref_allreps, file = "/Users/vinitaperiwal/GrowthCurver/normAUC_0.7", sep = "\t", quote = FALSE, row.names = FALSE)

# count wells
ref<-ref_allreps %>% filter(count >= 3)
head(ref)
nrow(ref) #7,372
ref_allreps_control<-unique(ref[,c(1,3,5)])
head(ref_allreps_control)
nrow(ref_allreps_control) #2,367
write.table(ref_allreps_control, file = "/Users/vinitaperiwal/GrowthCurver/normAUC_0.7_uniq", sep = "\t", quote = FALSE, row.names = FALSE)

#plot count of ref wells in all bugs
ref_allreps_control<-merge(ref_allreps_control, bug_desc, by="Bug_ID")
head(ref_allreps_control)

CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/ref_wells.svg", width = 4, height = 3, bg = "white") 
ref_allreps_control %>% dplyr::group_by(Bug_ID,Sp_short) %>% summarise(total = length(well)) %>%
  ggplot(aes(x=Sp_short,y=total)) + geom_bar(aes(fill=""), stat = "identity") + 
  scale_fill_jama() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + scale_y_continuous(name = "# Ref wells") +
  th + geom_text(aes(label=total), color="black", size = 2.5, nudge_y = 10)
dev.off()  

# compute mean of reference wells
reference_wells<-ref %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
  summarise(mean_ref_wells = mean(auc_l))
nrow(reference_wells) #41
head(reference_wells)

#normalize other plates by reference wells
SFall_others<-SFall %>% filter(Drug_name != 'control')
nrow(SFall_others) #38,976

ref<-merge(SFall_others, reference_wells, by = c('Bug_ID','Replicate_no')) %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
  mutate(normAUC = auc_l/mean_ref_wells)
ref<-ref[,c(1:9,11)]
View(ref)
nrow(ref) #38,976

#merge all normalized
all_norm<-rbind(ref,normAUCl_food)
nrow(all_norm) #46,752
head(all_norm)

write.table(all_norm, file = "/Users/vinitaperiwal/GrowthCurver/Figures/final_normAUCs", sep = "\t", quote = FALSE, row.names = FALSE)

##
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC.svg", sep = ""), width = 8, height = 4.6, bg = "white")
all_norm %>% dplyr::group_by(Bug_ID,Replicate_no) %>% filter(normAUC < 2) %>%
  ggplot(aes(x=Sp_short,y=normAUC)) + geom_boxplot(aes(fill = Plate_no), outlier.size = 0.01, lwd=0.15) +
  th + scale_fill_d3(palette = "category20") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

########################### correlation analysis (for wells having normAUC < 1.1)

mat_melt<-all_norm %>% filter(normAUC < 1.1) %>% reshape2::melt() %>% reshape(idvar = c("Bug_ID","Plate_no","Drug_name","well","Phyla","Species","Sp_short","variable"), timevar = c("Replicate_no"), direction = "wide")
head(mat_melt)  
nrow(mat_melt)
colnames(mat_melt)<-c("Bug_ID","Plate_no","Drug_name","well","Phyla","Species","Sp_short","variable","Rep1","Rep2","Rep3","Rep4")

write.table(mat_melt, file = "/Users/vinitaperiwal/GrowthCurver/Figures/mat", sep = "\t", quote = FALSE, row.names = FALSE)

#Select variable to compute correlation
mat<-mat_melt %>% filter(variable == 'normAUC')
head(mat)
nrow(mat) 

bugs<-unique(mat$Bug_ID)
bugs

correl_df<-data.frame(Bug_ID=character(),
                      Plate_no=character(),
                      comparisons=character(),
                      value=double(),
                      Parameter=character(),
                      stringsAsFactors=FALSE)

for(i in 1:length(bugs)){
  
  f<-mat %>% filter(Bug_ID == bugs[i])
  head(f)
  
  g<-data.frame(f[colSums(!is.na(f)) > 0])
  head(g)
  nrow(g)
  
  #g<-g %>% drop_na()
  
  plates<-unique(g$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    c<-g %>% filter(Plate_no == plates[j])
    head(c)
    c<-na.omit(c)
    
    if(nrow(c) > 4){
      
      #print(nrow(c))   
      Bug_ID<-c[1,1]
      Bug_ID
      Plate_no<-plates[j]
      Plate_no
      
      select_numeric_columns<-c %>% select_if(., is.numeric) 
      
      correlations<-correlate(select_numeric_columns, method = "pearson")
      correlations
      corr<-melt(correlations)
      corr<-na.omit(corr)
      corr
      corr$comp<-paste(corr$rowname,corr$variable,sep = "-")
      corr<-corr[,c(4,3)]
      
      combine_df<-cbind(data.frame(Bug_ID,Plate_no), corr)
      print(combine_df)
      correl_df<-rbind(correl_df, data.frame(combine_df))
      
    }
  }
  
}

head(correl_df)

corr_values<-correl_df
head(corr_values)
nrow(corr_values) #756
corr_values<-corr_values %>% filter(!comp %in% c("Rep3-Rep2","Rep4-Rep2","Rep4-Rep3","Rep2-Rep1","Rep3-Rep1","Rep4-Rep1"))
nrow(corr_values) #378

write.table(corr_values, file = paste0("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

median_corr_bug_wise<-corr_values %>% dplyr::group_by(Bug_ID) %>%
  dplyr::summarise(median_corr = median(value))

median_corr_bug_wise

write.table(median_corr_bug_wise, file = paste0("/Users/vinitaperiwal/GrowthCurver/median_correl_pearson"), sep = "\t", quote = FALSE, row.names = FALSE)

#bug names
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc)

all_bugs_desc<-merge(median_corr_bug_wise,bug_desc, by="Bug_ID")
head(all_bugs_desc)

all_bugs_corr_desc<-merge(corr_values,bug_desc, by="Bug_ID")
head(all_bugs_corr_desc)

#plots colored histogram of all correlation values
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC_correlation.svg", sep = ""), width = 6, height = 4, bg = "white")
ggplot(all_bugs_corr_desc, aes(value)) + geom_histogram(aes(fill=Sp_short), color="white") + theme_bw() + th + scale_fill_d3(palette = "category20")
dev.off()

#plots median correlation values of plates of each bug
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/median_correlation.svg", sep = ""), width = 4, height = 4, bg = "white")
all_bugs_desc %>% dplyr::group_by(Bug_ID,Sp_short) %>% 
  ggplot(aes(x=Sp_short,y=median_corr)) + geom_bar(stat = "identity") + theme_bw() + th + 
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) + scale_x_discrete(name = "Species")
dev.off()

###############################################################################

###### fit t-distribution per bug plate wise (all replicates pooled)
ref_all_wells<-ref_allreps %>% filter(count >= 3)
head(ref_all_wells)
nrow(ref_all_wells) #7372

ref_all_wells<-ref_all_wells[,c(1:10)]
ref_all_wells["label"]<-"reference"
head(ref_all_wells)
nrow(ref_all_wells) #7372
ncol(ref_all_wells) #11

###########################################
final_normAUCs<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/final_normAUCs", sep = "\t", header = TRUE)
head(final_normAUCs)
nrow(final_normAUCs) #46,752

data_foodplates<-final_normAUCs %>% filter(Drug_name == 'control')
nrow(data_foodplates) #7776
head(data_foodplates)

samples_data_foodplates<-anti_join(data_foodplates, ref_all_wells, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","auc_l","normAUC","Phyla","Species","Sp_short"))
nrow(samples_data_foodplates) #404
head(samples_data_foodplates)

samples_data_foodplates["label"]<-"sample"
head(samples_data_foodplates)
ncol(samples_data_foodplates)

total_list<-rbind(data.frame(ref_all_wells),data.frame(samples_data_foodplates))
head(total_list)
nrow(total_list) #7776
View(total_list)
ncol(total_list) #11

#### computing p-values from student distribution
# total_list<-total_list %>% filter(normAUC < 1.5)
# nrow(total_list) #7442
#Replicate wise
hits_rep<-total_list %>% dplyr::group_by(Bug_ID,Replicate_no) %>% do(compute_pval(.))
head(hits_rep)
nrow(hits_rep) #7776

#Replicate and Plate wise
hits_rep_plate<-total_list %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% do(compute_pval(.))
head(hits_rep_plate)
nrow(hits_rep_plate) #7776

# take max of both pvals (conservative estimate)
max_pvals<-merge(hits_rep,hits_rep_plate,by=c("Bug_ID","Replicate_no","Plate_no",
                                                  "Drug_name","well","auc_l","Phyla",
                                                  "Species","Sp_short","normAUC","label"))


max_pvals$pv<-pmax(max_pvals$pv.x,max_pvals$pv.y)
View(max_pvals)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_nullhypothesis.svg", sep = ""), width = 3, height = 2, bg = "white")
max_pvals %>% filter(label=='reference') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + th
dev.off()

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_unusedcontrols.svg", sep = ""), width = 3, height = 2, bg = "white")
max_pvals %>% filter(label=='sample' & well == 'H12') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + th
dev.off()

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_alternatehypothesis.svg", sep = ""), width = 3, height = 2, bg = "white")
max_pvals %>% filter(label=='sample' & well != 'H12') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + th
dev.off()

write.table(max_pvals, file = "/Users/vinitaperiwal/GrowthCurver/max_pvals", sep = "\t", quote = FALSE, row.names = FALSE)

### combining p values across replicates using fisher's method

head(max_pvals)
nrow(max_pvals) #7,776

combined_pval<-max_pvals %>% dplyr::group_by(Bug_ID,Plate_no,Drug_name,well,Phyla,Species,Sp_short) %>%
  summarise(combined_pv = sumlog(pv)[["p"]], meannormAUC = mean(normAUC))

head(combined_pval)
nrow(combined_pval) #2,496

# significant wells
hits<-combined_pval %>% filter(combined_pv < 0.05)
head(hits)
nrow(hits) #250

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_combined.svg", sep = ""), width = 4, height = 2, bg = "white")
hits %>% ggplot(aes(combined_pv)) + geom_histogram(aes(fill=Plate_no)) + 
  th + scale_fill_jama()
dev.off()

### multiple hypotheses testing: error correction

p_bh<-combined_pval %>% dplyr::group_by(Bug_ID) %>%
  mutate(p_bh = p.adjust(combined_pv,method = "BH"))
View(p_bh)

hits_bh<-p_bh %>% filter(p_bh < 0.05)
nrow(hits_bh)
View(hits_bh)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_bh.svg", sep = ""), width = 4, height = 2, bg = "white")
hits_bh %>% ggplot(aes(p_bh)) + geom_histogram(aes(fill=Plate_no)) + 
  th + scale_fill_jama()
dev.off()

## plot hit wells
# read food description plate wise (Misc excel sheet)
annot<-as_tibble(read.table("/Users/vinitaperiwal/GrowthCurver/Figures/food_desc", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(annot)
nrow(annot) #1,330

#merge with hit wells
hits_bh_annot<-merge(hits_bh, annot, by=c("Drug_name","Plate_no","well"))
View(hits_bh_annot)
nrow(hits_bh_annot) #76

write.table(hits_bh_annot, file = "/Users/vinitaperiwal/GrowthCurver/Figures/hits_bh_annot", sep = "\t", quote = FALSE, row.names = FALSE)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/p_bh_hits_heatmap.svg", sep = ""), width = 6, height = 7, bg = "white")
hits_bh_annot %>% ggplot(aes(x=Sp_short,y=Product.name)) + geom_tile(aes(fill=p_bh)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_viridis_b() + scale_y_discrete(name = "Food compound") 
dev.off()

# ggplot2::qplot(rank(pvals), pvals, xlab = "p-value rank", ylab = "p-values") + 
#   geom_abline(intercept = 0, slope = alpha/2496, aes(color="red")) +
#   ylim(c(0,0.2))


############## Bliss interactions

#normAUCs food compounds
SFa<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well != 'H12' & Plate_no %in% c("plate12","plate13")) %>%
  mutate(SFa = normAUC)
head(SFa)  
nrow(SFa) #7695
#merge(params, control_wells, by = c('Bug_ID','Replicate_no')) %>% dplyr::group_by(Bug_ID) %>% filter(well != 'H12' & Plate_no %in% c("plate12","plate13")) %>%
#mutate(SFa = auc_l/max_aucl_h12_wells)

#normAUCs drugs
SFq<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well == 'H12' & !(Plate_no %in% c("plate12","plate13","plate16"))) %>%
  mutate(SFq = normAUC)
head(SFq)
nrow(SFq) #406

###################### boxplot of effect of HTD drugs alone (compared with bug alone)
SF_drugs_bugs<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well == 'H12' & Plate_no != "plate16") %>%
  mutate(SFq = normAUC)
View(SF_drugs_bugs)
nrow(SF_drugs_bugs) #487

write.table(SF_drugs_bugs, file = "/Users/vinitaperiwal/GrowthCurver/SF_drugs_bugs_desc", sep = "\t", quote = FALSE, row.names = FALSE)

# box plot for all H12 wells by h12-control
CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/SF_drugs_bugs_desc.svg", width = 11, height = 6, bg = "white")
SF_drugs_bugs %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
  ggplot(aes(x=Species,y=SFq)) + geom_boxplot(aes(fill = Drug_name), outlier.size=1, lwd=0.3) +
  facet_wrap("Bug_ID", scales = "free", nrow = 3) + th + scale_fill_nejm() + scale_x_discrete(name = "Species") +
  scale_y_continuous(name = "normalized AUC")
dev.off()

############################################################
#normAUCs food+drug
SFaq<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well != 'H12' & Drug_name != 'control') %>%
  mutate(SFaq = normAUC)
head(SFaq)
nrow(SFaq) #38,570

##### start joining SFaq and SFa

#merge promazine
#change plate numbers 12 and 13 to plate 2 and 3 respectively

p<-SFa
head(SFa)
p[p=="plate12"]<-"plate2"
p[p=="plate13"]<-"plate3"
colnames(p)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(p) #7,695

merge_proma<-merge(SFaq, p, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_proma)
nrow(merge_proma) #7,695

#merge loxapine
#change plate numbers 12 and 13 to plate 4 and 5 respectively

l<-SFa
l[l=="plate12"]<-"plate4"
l[l=="plate13"]<-"plate5"
colnames(l)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(l) #7,695

merge_loxa<-merge(SFaq, l, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_loxa)
nrow(merge_loxa) #7,505

#merge omeprazole
#change plate numbers 12 and 13 to plate 6 and 7 respectively

o<-SFa
o[o=="plate12"]<-"plate6"
o[o=="plate13"]<-"plate7"
colnames(o)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(o) #7,695

merge_ome<-merge(SFaq, o, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_ome)
nrow(merge_ome) #7,695

#merge felodipine
#change plate numbers 12 and 13 to plate 8 and 9 respectively

f<-SFa
f[f=="plate12"]<-"plate8"
f[f=="plate13"]<-"plate9"
colnames(f)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(f) #7,695

merge_felo<-merge(SFaq, f, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_felo)
nrow(merge_felo) #7,695

#merge 5fu
#change plate numbers 12 and 13 to plate 10 and 11 respectively

fu<-SFa
fu[fu=="plate12"]<-"plate10"
fu[fu=="plate13"]<-"plate11"
colnames(fu)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(fu) #7,69

merge_fu<-merge(SFaq, fu, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_fu)
nrow(merge_fu) #7,410

#merge duloxetine
#change plate numbers 12 and 13 to plate 10 and 11 respectively

dulo<-SFa
dulo[dulo=="plate12"]<-"plate14"
dulo[dulo=="plate13"]<-"plate15"
colnames(dulo)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","SFa_auc_l","Phyla","Species","Sp_short","SFa_normAUC","SFa")
nrow(dulo) #7,695

merge_dulo<-merge(SFaq, dulo, by=c("Bug_ID","Replicate_no","Plate_no","well","Phyla","Species","Sp_short"))
head(merge_dulo)
nrow(merge_dulo) #570


merge_SFaq_SFa<-rbind(merge_proma,merge_loxa,merge_ome,merge_felo,merge_fu,merge_dulo)
colnames(merge_SFaq_SFa)[8]<-"Drug_name"
head(merge_SFaq_SFa)
nrow(merge_SFaq_SFa) #38,570

merge_SFaq_SFa_SFq<-merge(merge_SFaq_SFa, SFq, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name","Phyla","Species","Sp_short"))
head(merge_SFaq_SFa_SFq)
nrow(merge_SFaq_SFa_SFq) #38,570

colnames(merge_SFaq_SFa_SFq)<-c("Bug_ID","Replicate_no","Plate_no","Drug_name",
                                "Phyla","Species","Sp_short","SFaq_well","SFaq_aucl",
                                "SFaq_normAUC",
                                "SFaq","SFa_Drug_name","SFa_aucl",
                                "SFa_normAUC",
                                "SFa","SFq_well","SFq_aucl",
                                "SFq_normAUC",
                                "SFq")
head(merge_SFaq_SFa_SFq)

#bliss<-merge_SFaq_SFa_SFq %>% mutate(bliss_ex = SFa*SFq, bliss_score = SFaq-bliss_ex, label = ifelse(SFaq < bliss_ex, "synergy","antagonism"))
bliss<-merge_SFaq_SFa_SFq %>% mutate(bliss_ex = SFa*SFq, bliss_score = SFaq-bliss_ex)
head(bliss)

write.table(bliss, file = "/Users/vinitaperiwal/GrowthCurver/Figures/bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

############# reading annotations to compound wells #########################
# read food description plate wise (Misc excel sheet)
bliss<-read.table("/Users/vinitaperiwal/GrowthCurver/Figures/bliss_scores", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
annot<-as_tibble(read.table("/Users/vinitaperiwal/GrowthCurver/Figures/food_desc", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(annot)
nrow(annot) #1,330

colnames(bliss)[8]<-"well"
nrow(bliss) #38,570
head(bliss)
#merge with bliss desc
bliss_annot<-merge(bliss, annot, by=c("Drug_name","Plate_no","well"))
View(bliss_annot)
nrow(bliss_annot) #38,570

write.table(bliss_annot, file = "/Users/vinitaperiwal/GrowthCurver/Figures/bliss_annot", sep = "\t", quote = FALSE, row.names = FALSE)

#statistical determination of synergy/antagonism

short_blissannot<-bliss_annot[,c(1:8,11,15,19,22:25)]
View(short_blissannot)
nrow(short_blissannot) #115710

A<-short_blissannot %>% dplyr::group_by(Bug_ID,Plate_no,Drug_name,well) %>% 
  mutate(n1 = length(SFa), lSFa = log(SFa), y1=mean(lSFa), s1=var(lSFa)*(n1-1),
         n2 = length(SFq), lSFq = log(SFq), y2=mean(lSFq), s2=var(lSFq)*(n2-1),
         n3 = length(SFaq), lSFaq = log(SFaq), y3=mean(lSFaq), s3=var(lSFaq)*(n3-1),
         sy=s1+s2+s3, dft=n1+n2+n3-3, denf=1/n1+1/n2+1/n3,
         tss=(y1+y2-y3)/sqrt(sum(sy)/dft)/denf,
         bliss=
         pv=2*(1-pt(abs(tss),df=dft)),
         pvP=1-pt(tss,df=dft)) #n - no of observations (replicates), sy - total sum of squares, dft,denf - df, tss - t-statistic, pv - p-value for Bliss independence hypothesis, pvP - One-sided p-value

View(A)





# density to determine bliss score distribution
CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/bliss_density.svg", width = 9, height = 4, bg = "white")
bliss_annot %>% ggplot(aes(bliss_score)) + geom_density(aes(fill=Drug_name, color=Drug_name), alpha=0.8) +
  facet_wrap("Sp_short", nrow = 3, scales = "free") + theme_bw() +
  scale_fill_jama() + scale_color_jama() + th
dev.off()

#statistical determination of synergy/antagonism



log_bliss<-bliss_annot %>% dplyr::group_by(Bug_ID) %>%
  mutate(logSFa = log10(SFa), logSFq = log10(SFq), logSFaq = log10(SFaq))

View(log_bliss)


# expected vs observed viability
CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/bliss_obs_exp.svg", width = 9, height = 6, bg = "white")
bliss_annot %>% ggplot(aes(x=SFaq,y=bliss_ex)) + geom_point(aes(color=Drug_name), size=0.5) + geom_abline() + theme_bw() +
  scale_y_continuous(name = "Expected viability (V1*V2)") + 
  scale_x_continuous(name = "Observed viability (V12)") + scale_color_jama() + 
  facet_wrap(~Sp_short, scales = "free") + th
dev.off()

#heatmap
bliss_annot %>% ggplot(aes(x=Drug_name,y=Product.name)) + facet_grid(~Sp_short) + geom_tile(aes(fill=bliss_score)) + theme_bw() +
  th + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_gradient2(low="#67001f",high="#1a1a1a") + scale_y_discrete(name = "Food compounds") 

#####references
View(bliss_annot)

ref_bliss_annot<-bliss_annot %>% filter(bliss_score > -0.1 & bliss_score < 0.1) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>% mutate(count = length(well)) 
View(ref_bliss_annot)
nrow(ref_bliss_annot) #26,719

bliss_ref<-ref_bliss_annot %>% filter(count >= 3)
head(bliss_ref)
nrow(bliss_ref) #15,629

bliss_ref_wells<-bliss_ref[,c(1:8,11:12,15:16,19:21,23)]
bliss_ref_wells["outcome"]<-"reference"
View(bliss_ref_wells)
nrow(bliss_ref_wells) #15629
ncol(bliss_ref_wells) #17

##samples

bliss_others<-bliss_annot[,c(1:8,11:12,15:16,19:21,23)]

samples_bliss<-anti_join(bliss_others, bliss_ref, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","Phyla","Species","Sp_short","SFaq","SFa_Drug_name","SFa","SFq_well","SFq","bliss_ex","bliss_score","Product.name"))
nrow(samples_bliss)  #22941
samples_bliss["outcome"]<-"sample"
ncol(samples_bliss) #17
head(samples_bliss)


total_list_bliss<-rbind(data.frame(bliss_ref_wells),data.frame(samples_bliss))
head(total_list_bliss)
nrow(total_list_bliss) #38,570
View(total_list_bliss)
ncol(total_list_bliss) #17

## compute pval

#synergies
syn_bliss<-bliss_annot %>% filter(bliss_score < 0) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>% mutate(count = length(well))
nrow(syn_bliss) #24,069
head(syn_bliss)

syn_bliss_ref<-syn_bliss %>% filter(bliss_score < -0.1 & count >= 3)
head(syn_bliss_ref)
nrow(syn_bliss_ref) #5287

syn_ref_wells<-syn_bliss_ref[,c(1:8,11:12,15:16,19:21,23)]
syn_ref_wells["outcome"]<-"reference"
View(syn_ref_wells)
nrow(syn_ref_wells) #5287
ncol(syn_ref_wells) #17

##samples

syn_others<-syn_bliss[,c(1:8,11:12,15:16,19:21,23)]

syn_samples<-anti_join(syn_others, syn_ref_wells, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","Phyla","Species","Sp_short","SFaq","SFa_Drug_name","SFa","SFq_well","SFq","bliss_ex","bliss_score","Product.name"))
nrow(syn_samples)  #18782
syn_samples["outcome"]<-"sample"
ncol(syn_samples) #17
head(syn_samples)


total_syn<-rbind(data.frame(syn_ref_wells),data.frame(syn_samples))
head(total_syn)
nrow(total_syn) #24,069
View(total_syn)
ncol(total_syn) #17

#Bug wise
syn_hits_rep<-total_syn %>% dplyr::group_by(Bug_ID) %>% do(compute_pval_syn(.))
head(syn_hits_rep)
nrow(syn_hits_rep) #24,069

#Bug and Replicate wise
syn_hits_rep_plate<-total_syn %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
  do(compute_pval_syn(.))
head(syn_hits_rep_plate)
nrow(syn_hits_rep_plate) #24,069

# take max of both pvals (conservative estimate)
syn_max_pvals<-merge(syn_hits_rep,syn_hits_rep_plate,by=c("Bug_ID","Replicate_no","Plate_no",
                                              "Drug_name","well","Phyla",
                                              "Species","Sp_short","bliss_score","outcome"))


syn_max_pvals$pv<-pmax(syn_max_pvals$pv.x,syn_max_pvals$pv.y)
syn_max_pvals$label<-"synergy"
syn_max_pvals<-syn_max_pvals[,c(1:10,25,27:28)]
View(syn_max_pvals)

syn_ref<-syn_hits_rep %>% filter(outcome=='reference') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + ggtitle(label = "synergies ref") + th
syn_samp<-syn_hits_rep %>% filter(outcome=='sample') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + ggtitle(label = "synergies hits") + th


##antagonisms
ant_bliss<-total_list_bliss %>% filter(bliss_score > 0)
nrow(ant_bliss) #14,499

#Bug wise
ant_hits_rep<-ant_bliss %>% dplyr::group_by(Bug_ID) %>% do(compute_pval_bliss_ant(.))
head(ant_hits_rep)
nrow(ant_hits_rep) #14,499

#Bug and Replicate wise
ant_hits_rep_plate<-ant_bliss %>% dplyr::group_by(Bug_ID,Replicate_no) %>%
  do(compute_pval_bliss_ant(.))
View(ant_hits_rep_plate)
nrow(ant_hits_rep_plate) #14,499

# take max of both pvals (conservative estimate)
ant_max_pvals<-merge(hits_rep,hits_rep_plate,by=c("Bug_ID","Replicate_no","Plate_no",
                                              "Drug_name","well","Phyla",
                                              "Species","Sp_short","bliss_score","outcome"))


ant_max_pvals$pv<-pmax(ant_max_pvals$pv.x,ant_max_pvals$pv.y)
ant_max_pvals$label<-"antagonism"
ant_max_pvals<-ant_max_pvals[,c(1:10,27,29:30)]
View(ant_max_pvals)

ant_ref<-ant_max_pvals %>% filter(outcome=='reference') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + ggtitle(label = "antagonisms ref") + th
ant_samp<-ant_max_pvals %>% filter(outcome=='sample') %>% ggplot(aes(x=pv)) + geom_histogram(color="white",fill="#4d4d4d",bins = 30) + ggtitle(label = "antagonisms hits") + th


### plot syn and ant pval distribution

CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/syn_ant_pvals.svg", width = 5, height = 4, bg = "white")
grid.arrange(syn_ref, ant_ref, syn_samp, ant_samp, nrow = 2)
dev.off()

### combining p values across replicates using fisher's method

interactions<-rbind(data.frame(syn_max_pvals),data.frame(ant_max_pvals))
  
View(interactions)
nrow(interactions) #38,568

combined_pval<-interactions %>% dplyr::group_by(Bug_ID,Plate_no,Drug_name,well,Phyla,Species,Sp_short) %>%
    summarise(combined_pv = sumlog(pv)[["p"]])

head(combined_pval)
nrow(combined_pval) #2,496

# significant wells
hits<-combined_pval %>% filter(combined_pv < 0.05)
head(hits)
nrow(hits) #250

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_combined.svg", sep = ""), width = 4, height = 2, bg = "white")
hits %>% ggplot(aes(combined_pv)) + geom_histogram(aes(fill=Plate_no)) + 
  th + scale_fill_jama()
dev.off()

### multiple hypotheses testing: error correction

p_bh<-combined_pval %>% dplyr::group_by(Bug_ID) %>%
  mutate(p_bh = p.adjust(combined_pv,method = "BH"))
View(p_bh)

hits_bh<-p_bh %>% filter(p_bh < 0.05)
nrow(hits_bh)
View(hits_bh)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pval_bh.svg", sep = ""), width = 4, height = 2, bg = "white")
hits_bh %>% ggplot(aes(p_bh)) + geom_histogram(aes(fill=Plate_no)) + 
  th + scale_fill_jama()
dev.off()

## plot hit wells
# read food description plate wise (Misc excel sheet)
annot<-as_tibble(read.table("/Users/vinitaperiwal/GrowthCurver/Figures/food_desc", header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(annot)
nrow(annot) #1,330

#merge with hit wells
hits_bh_annot<-merge(hits_bh, annot, by=c("Drug_name","Plate_no","well"))
head(hits_bh_annot)
nrow(hits_bh_annot) #76

write.table(hits_bh_annot, file = "/Users/vinitaperiwal/GrowthCurver/Figures/hits_bh_annot", sep = "\t", quote = FALSE, row.names = FALSE)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/p_bh_hits_heatmap.svg", sep = ""), width = 6, height = 7, bg = "white")
hits_bh_annot %>% ggplot(aes(x=Sp_short,y=Product.name)) + geom_tile(aes(fill=p_bh)) + theme_bw() +
  th + theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_viridis_b() + scale_y_discrete(name = "Food compound") 
dev.off()



##################### plotting
bliss_annot<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/bliss_annot", sep = "\t", header = TRUE)
head(bliss_annot)

#bliss_annot_heat<-bliss_annot[,c(1,4,5,24,25,28,30)]
#head(bliss_annot_heat)

CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/bliss_heatmap.svg", width = 10, height = 6, bg = "white")
bliss_annot %>% filter(Drug_name != "duloxetine") %>% ggplot(aes(x=Drug_name,y=Product.name)) + geom_tile(aes(fill=bliss_score)) +
  facet_wrap("Sp_short", nrow = 3) + th + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90), axis.ticks = element_blank()) + 
  scale_fill_viridis_b() + scale_y_discrete(name = "Food compound")
dev.off()

#duloxetine
CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/bliss_duloxetine.svg", width = 5, height = 4, bg = "white")
bliss_annot %>% filter(Drug_name == "duloxetine") %>% ggplot(aes(x=Drug_name,y=Product.name)) + geom_tile(aes(fill=bliss_score)) +
  facet_wrap("Sp_short") + th + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank()) + 
  scale_fill_viridis_b()
dev.off()

# defining synergies and antagonisms by threshold

### synergies
View(bliss_desc)
synergies<-bliss_annot %>% filter(bliss_score < -0.1) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>%
  mutate(wells_count = length(well)) %>% filter(wells_count > 2)
View(synergies)
nrow(synergies) #577

syn_num<-unique(synergies[,c(1:4,6:8,22,27)])
nrow(syn_num) #186
head(syn_num)

write.table(syn_num, file = "/Users/vinitaperiwal/GrowthCurver/Figures/synergies.txt", row.names = FALSE, quote = FALSE, sep = '\t')

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/synergies.svg", sep = ""), width = 12, height = 7, bg = "white")
unique(synergies[,c(1:4,6:8,21,22,24,27)]) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>%
  group_by(Bug_ID) %>%
  ggplot(aes(x=Drug_name,y=bliss_score)) + geom_point(aes(color=Drug_name), size=1) + facet_wrap("Sp_short", scales = "free", nrow = 3) + 
  scale_color_nejm() + scale_fill_nejm() + th + theme(axis.text.x = element_blank()) + scale_x_discrete(name = "") +
  scale_y_continuous(name = "synergistic interactions") + geom_text(aes(label=Product.name), size = 2, check_overlap=TRUE)
dev.off()

### antagonisms
antagonisms<-bliss_annot %>% filter(bliss_score > 0.1) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>%
  mutate(wells_count = length(well)) %>% filter(wells_count > 2)
head(antagonisms)
nrow(antagonisms) #301

ant_num<-unique(antagonisms[,c(1:4,6:8,22,27)])
nrow(ant_num) #100
head(ant_num)

write.table(ant_num, file = "/Users/vinitaperiwal/GrowthCurver/Figures/antagonisms.txt", row.names = FALSE, quote = FALSE, sep = '\t')

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/antagonisms.svg", sep = ""), width = 9, height = 6, bg = "white")
unique(antagonisms[,c(1:4,6:8,21,22,24,27)]) %>% dplyr::group_by(Bug_ID,Plate_no,well) %>%
  group_by(Bug_ID) %>%
  ggplot(aes(x=Drug_name,y=bliss_score)) + geom_point(aes(color=Drug_name), size=1.5) + facet_wrap("Sp_short", scales = "free", nrow = 3) + 
  scale_color_nejm() + th + scale_fill_nejm() + theme(axis.text.x = element_blank()) +
  scale_y_continuous(name = "antagonistic interactions") + geom_text(aes(label=Product.name), size = 2.5, check_overlap=TRUE)
dev.off()

#############
syn<-syn_num
syn$label<-"synergy"
head(syn)

ant<-ant_num
ant$label<-"antagonism"
head(ant)

syn_ant<-rbind(syn,ant)
head(syn_ant)
nrow(syn_ant)

syn_ant<-syn_ant %>% dplyr::group_by(Bug_ID,Sp_short,label) %>% summarise(count_label = length(label)) %>% 
  dplyr::group_by(Bug_ID) %>% mutate(sum = sum(count_label))
View(syn_ant)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/num_syn_ant.svg", sep = ""), width = 5, height = 3, bg = "white")
syn_ant %>%  
  ggplot(aes(x=Sp_short,y=count_label, fill=label)) + geom_bar(stat = "identity") + theme_bw() + th + 
  scale_fill_manual(values = c("#ce1256","#00441b")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(name = "count")
dev.off()  

#CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/total_syn_ant.svg", sep = ""), width = 5, height = 3)
syn_ant %>%  
  ggplot(aes(x=label,y=count_label, fill=label)) + geom_bar(stat = "identity") + theme_bw() + th + 
  scale_fill_manual(values = c("#ce1256","#00441b")) +
  scale_y_continuous(name = "count") + theme(legend.position = "none")
#dev.off()



# ######### boxplot of effect of food compounds alone: 2 - using cook's distance
# #### using cook's distance
# library(olsrr)
# 
# food_cmpds<-h12_normAUCl_bugs_desc %>% filter(Drug_name == 'control') %>% melt() %>% filter(variable == 'h12_normAUCl')
# head(food_cmpds)
# nrow(food_cmpds) #7,680
# 
# bugs_cook<-data.frame(stringsAsFactors = FALSE)
# bugs<-unique(food_cmpds$Bug_ID)
# bugs
# 
# for(i in 1:length(bugs)){
#   
#   one_bug<-data.frame(food_cmpds %>% filter(Bug_ID == bugs[i]))
#   head(one_bug)
#   
#   cooks<-glm(value ~ Plate_no + Replicate_no, data = one_bug)
#   one_bug$cooksd <- cooks.distance(cooks)
#   head(one_bug)
#   nrow(one_bug)
#   
#   #define outlier
#   one_bug$outlier<-ifelse(one_bug$cooksd >= 2*mean(one_bug$cooksd), "delete","keep")
#   bugs_cook<-rbind(bugs_cook, data.frame(one_bug))
# 
# }
# 
# head(bugs_cook)
# nrow(bugs_cook)
# 
# #label outliers
# bugs_cook_outliers<-bugs_cook %>% dplyr::group_by(Bug_ID,Plate_no,well,outlier) %>% 
#   mutate(count_of_wells = length(which(outlier=='delete')), label = as.character(lapply(count_of_wells, function(x){
#     if(x %in% c(3,4)){
#       'High-conf'
#     }else if(x == 2){
#       'Med-conf'
#     }else if(x == 1){
#       'Low-conf'
#     }
#   })))
# 
# head(bugs_cook_outliers)
# nrow(bugs_cook_outliers) #7,680
# 
# out<-bugs_cook_outliers %>% filter(label %in% c('High-conf','Med-conf'))
# nrow(out)
# 
# write.table(bugs_cook_outliers, file = "/Users/vinitaperiwal/GrowthCurver/Figures/cooks_outliers", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/cooks_outliers.svg", width = 10, height = 6, bg = "white")
# bugs_cook_outliers %>% mutate(wells = ifelse(outlier=='delete', as.character(well), "")) %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% filter(!label %in% c('NA','','Low-conf')) %>%
#   ggplot(aes(Replicate_no,value)) + geom_boxplot(aes(fill=Plate_no), outlier.size = 0.7, lwd=0.3) + th +
#   facet_wrap("Species", scales = "free") + theme(axis.text.x=element_blank(), axis.title.x = element_blank()) + scale_y_continuous(name = "normalized AUC") +
#   scale_color_nejm() + geom_text(aes(label=wells, color=label), size = 2.5) + scale_fill_jama(alpha = 0.8) 
# dev.off()

##########################


########################################## Reference wells (just to check)


################################################### data filtering
head(params_out) 
nrow(params_out) #54,912




#filtering high sigma values (bug-wise) based on density distribution

#plot of sigma values (gives quality of fit)
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/sigma.svg", sep = ""), width = 10, height = 8, bg = "white") 
# params_doses %>% melt() %>% filter(variable == 'sigma') %>%
#   dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Drug_name,well) %>%
#   ggplot(aes(value)) + geom_density(aes(color=Replicate_no)) + th +
#   scale_color_nejm() + facet_wrap("Bug_ID", scales = "free")
# dev.off()

#identify sigma outliers (bug-wise)
# bugs<-unique(params_doses$Bug_ID)
# bugs
# 
# #detect outliers (auc_l and r)
# group_outliers<-params_doses %>% dplyr::group_by(Bug_ID,Plate_no,well) %>% 
#   mutate(outlier_a = ifelse(is_outlier(auc_l), as.character(well), ""), outlier_r = ifelse(is_outlier(sigma), as.character(well), ""))
# head(group_outliers)
# nrow(group_outliers) #46,080
# 
# #save outliers
# write.table(group_outliers, file = "/Users/vinitaperiwal/GrowthCurver/Figures/outliers", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# #filter outliers (auc_l and r)
# fit_r<-data.frame(group_outliers %>% dplyr::group_by(Bug_ID) %>% filter(outlier_a == '' & outlier_r == ''))
# nrow(fit_r) #43,902
# head(fit_r)
# 
# # auc_l (all vs filtered)
# params_doses %>% melt() %>% filter(variable == 'auc_l') %>%
#   dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Drug_name,well) %>%
#   ggplot(aes(value)) + geom_density(aes(color=Replicate_no)) + th +
#   scale_color_nejm() + facet_wrap("Bug_ID", scales = "free")
# 
# 
# fit_r %>% melt() %>% filter(variable == 'auc_l') %>%
#   dplyr::group_by(Bug_ID,Replicate_no,Plate_no,Drug_name,well) %>%
#   ggplot(aes(value)) + geom_density(aes(color=Replicate_no)) + th +
#   scale_color_nejm() + facet_wrap("Bug_ID", scales = "free")

#applied clustering on each well's parameters for all replicates
# used different script "hclust_gc.R"

# View(params)
# nrow(params) #54,912
# 
# #remove outlier wells obtained from clustering (use all thresholds)
# #read all files
# file.names <- list.files(path = '/Users/vinitaperiwal/GrowthCurver/outliers', recursive = TRUE, pattern = "^outliers_dist_") #recursive reads through all subfolders and files
# file.names
# 
# #loop for each .tab file, fits logistic curve, annotates
# for(i in 1:length(file.names)){
#   
#   #read file/table
#   input_file<-tools::file_path_sans_ext(file.names[i]) #read filename w/o extension
#   input_file
#   
#   outlier_from_hclust<-read.table(paste0("outliers/",file.names[i]), sep = "\t", header = TRUE)
#   head(outlier_from_hclust) 
#   colnames(outlier_from_hclust)<-c("Bug_ID","Plate_no","well","num_wells","count_dup","Replicate_no","outwell_count")
#   nrow(outlier_from_hclust) #7,086
#   
#   minus_outliers<-anti_join(params, outlier_from_hclust, by=c("Bug_ID","Plate_no","well","Replicate_no"))
#   View(minus_outliers)
#   nrow(minus_outliers) 
#   
#   head(minus_outliers)
#   
#   mat_melt<-minus_outliers[,1:20] %>% melt() %>% reshape(idvar = c("Bug_ID","Plate_no","Drug_name","well","variable"), timevar = c("Replicate_no"), direction = "wide")
#   head(mat_melt)  
#   nrow(mat_melt)
#   colnames(mat_melt)<-c("Bug_ID","Plate_no","Drug_name","well","variable","Rep1","Rep2","Rep3","Rep4")
#   
#   #write.table(mat_melt, file = "/Users/vinitaperiwal/GrowthCurver/Figures/mat", sep = "\t", quote = FALSE, row.names = FALSE)
#   
#   #Select variable to compute correlation
#   mat<-mat_melt %>% filter(variable == 'auc_l')
#   head(mat)
#   nrow(mat) 
#   
#   bugs<-unique(mat$Bug_ID)
#   bugs
#   
#   correl_df<-data.frame(Bug_ID=character(),
#                         Plate_no=character(),
#                         comparisons=character(),
#                         value=double(),
#                         Parameter=character(),
#                         stringsAsFactors=FALSE)
#   
#   for(i in 1:length(bugs)){
#     
#     f<-mat %>% filter(Bug_ID == bugs[i])
#     head(f)
#     
#     g<-data.frame(f[colSums(!is.na(f)) > 0])
#     head(g)
#     nrow(g)
#     
#     #g<-g %>% drop_na()
#     
#     plates<-unique(g$Plate_no)
#     plates
#     
#     for(j in 1:length(plates)){
#       
#       c<-g %>% filter(Plate_no == plates[j])
#       head(c)
#       c<-na.omit(c)
#       
#       if(nrow(c) > 4){
#         
#         #print(nrow(c))   
#         Bug_ID<-c[1,1]
#         Bug_ID
#         Plate_no<-plates[j]
#         Plate_no
#         
#         select_numeric_columns<-c %>% select_if(., is.numeric) 
#         
#         correlations<-correlate(select_numeric_columns, method = "pearson")
#         correlations
#         corr<-melt(correlations)
#         corr<-na.omit(corr)
#         corr
#         corr$comp<-paste(corr$rowname,corr$variable,sep = "-")
#         corr<-corr[,c(4,3)]
#         
#         combine_df<-cbind(data.frame(Bug_ID,Plate_no), corr)
#         print(combine_df)
#         correl_df<-rbind(correl_df, data.frame(combine_df))
#         
#       }
#     }
#     
#   }
#   
#   head(correl_df)
#   
#   corr_values<-correl_df
#   head(corr_values)
#   nrow(corr_values) #1,076
#   corr_values<-corr_values %>% filter(!comp %in% c("Rep3-Rep2","Rep4-Rep2","Rep4-Rep3","Rep2-Rep1","Rep3-Rep1","Rep4-Rep1"))
#   nrow(corr_values) #538
#   
#   write.table(corr_values, file = paste0("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_",input_file), sep = "\t", quote = FALSE, row.names = FALSE)
#   
#   median_corr_bug_wise<-corr_values %>% dplyr::group_by(Bug_ID) %>%
#     dplyr::summarise(median_corr = median(value))
#   
#   median_corr_bug_wise
#   
#   write.table(median_corr_bug_wise, file = paste0("/Users/vinitaperiwal/GrowthCurver/outliers/median_correl_pearson_",input_file), sep = "\t", quote = FALSE, row.names = FALSE)
#   
# }
# 
# # read median correlation values and plot them
# #read all files
# file_names <- list.files(path = '/Users/vinitaperiwal/GrowthCurver/outliers/', recursive = TRUE, pattern = "median_correl_pearson_") #recursive reads through all subfolders and files
# file_names
# 
# all_bugs<-data.frame(stringsAsFactors = FALSE)
# 
# for(i in 1:length(file_names)){
#   
#   #read file/table
#   file<-file_names[i]
#   file
#   
#   sp<-strsplit(file, split = "_", perl = TRUE)
#   sp
#   dist<-sp[[1]][6]
#   dist<-paste0("thresh_",dist)
#   dist
#   
#   median_correl<-read.table(paste0("outliers/",file), sep = "\t", header = TRUE)
#   
#   median_correl<-cbind(median_correl, dist)
#   head(median_correl)
#   
#   all_bugs<-rbind(all_bugs, data.frame(median_correl))
#   
# }
# 
# head(all_bugs)
# nrow(all_bugs) #120
# 
# #bug names
# bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
# head(bug_desc)
# nrow(bug_desc)
# 
# all_bugs_desc<-merge(all_bugs,bug_desc, by="Bug_ID")
# head(all_bugs_desc)
# 
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/outliers/median_corr_thresh.svg", sep = ""), width = 5, height = 4.1, bg = "white")
# ggplot(all_bugs_desc, aes(x=Sp_short,y=median_corr)) + geom_point(aes(fill=dist), size=2.5, shape=21) + 
#   scale_fill_d3() + theme_minimal() + th + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.1)) +
#   scale_x_discrete(name="Species") + scale_y_continuous(name = "median replicate correlation")
# dev.off()
# 
# # mix and merge median threshold for individual bugs for final thresholds
# file.names <- list.files(path = '/Users/vinitaperiwal/GrowthCurver/outliers/', recursive = TRUE, pattern = "^outliers_dist_") #recursive reads through all subfolders and files
# file.names
# 
# total_hclust_outliers<-data.frame(stringsAsFactors = FALSE)
# #loop for each .tab file, fits logistic curve, annotates
# for(i in 1:length(file.names)){
#   
#   #read file/table
#   input_file<-tools::file_path_sans_ext(file.names[i]) #read filename w/o extension
#   input_file
#   
#   sp<-strsplit(input_file, "_", perl = TRUE)
#   dist<-paste0("thresh_",sp[[1]][3])
#   
#   out_from_hclust<-read.table(paste0("outliers/",file.names[i]), sep = "\t", header = TRUE)
#   head(out_from_hclust) 
#   colnames(out_from_hclust)<-c("Bug_ID","Plate_no","well","num_wells","count_dup","Replicate_no","outwell_count")
#   nrow(out_from_hclust)
#   
#   out_from_hclust<-cbind(out_from_hclust,dist)
#   total_hclust_outliers<-rbind(total_hclust_outliers,out_from_hclust)
#   
# }
# 
# head(total_hclust_outliers)
# nrow(total_hclust_outliers) #65,405
# 
# total_hclust_outliers<-total_hclust_outliers %>% dplyr::group_by(Bug_ID) %>%
#   mutate(percent_out = (outwell_count/num_wells)*100)
# 
# #bug names
# bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
# head(bug_desc)
# nrow(bug_desc)
# 
# outliers_bugs_desc<-merge(total_hclust_outliers,bug_desc, by="Bug_ID")
# outliers_bugs_desc<-unique(outliers_bugs_desc[,c(1,8:12)])
# head(outliers_bugs_desc)
# 
# #merge with median_corr_values
# corr_percent<-merge(all_bugs_desc,outliers_bugs_desc,by=c("Bug_ID","dist","Phyla","Species","Sp_short"))
# head(corr_percent)
# 
# for_thresh<-merge(all_bugs_desc,for_thresh,by=c("Bug_ID","dist","Phyla","Species","Sp_short"))
# 
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/outliers/corr_and_percent_outliers.svg", sep = ""), width = 8, height = 5, bg = "white")
# corr_percent %>% 
#   ggplot(aes(x=percent_out,y=median_corr)) + geom_point(aes(fill=dist), shape=21, size=2.5) + 
#   scale_fill_jco() + theme_bw() + th +
#   scale_x_continuous(name="% outlier wells", limits = c(0,30,5)) + scale_y_continuous(name = "median correlation") +
#   facet_wrap("Bug_ID", scales = "free")
# dev.off()
# 
# #select thresholds
# bugs<-unique(corr_percent$Bug_ID)
# bugs
# thresholds<-data.frame(stringsAsFactors = FALSE)
# 
# for(i in 1:length(bugs)){
#   
#   bug<-bugs[i]
#   bug
#   
#   if(bug %in% c('NT5021','NT5022','NT5025','NT5026','NT5054')){
#     
#     thresh<-read.table("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_outliers_dist_3", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#     head(thresh)
#     thresh<-thresh %>% filter(Bug_ID == bug)
#     
#   }else if(bug %in% c('NT5004','NT5009','NT5028','NT5078')){
#     
#     thresh<-read.table("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_outliers_dist_4", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#     head(thresh)
#     thresh<-thresh %>% filter(Bug_ID == bug)
#     
#   }else if(bug %in% c('NT5001','NT5011')){
#     
#     thresh<-read.table("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_outliers_dist_5", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#     head(thresh)
#     thresh<-thresh %>% filter(Bug_ID == bug)
#     
#   }else if(bug %in% c('NT5019')){
#     
#     thresh<-read.table("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_outliers_dist_7", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#     head(thresh)
#     thresh<-thresh %>% filter(Bug_ID == bug)
#     
#   }else if(bug %in% c('NT5083')){
#     
#     thresh<-read.table("/Users/vinitaperiwal/GrowthCurver/outliers/auc_l_correl_pearson_outliers_dist_9", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#     head(thresh)
#     thresh<-thresh %>% filter(Bug_ID == bug)
#     
#   }
#   
#   thresholds<-rbind(thresholds,data.frame(thresh))
#   
# }
# 
# head(thresholds)
# thresholds<-merge(thresholds, bug_desc, by=c("Bug_ID"))
# 
# matrix_thresh<-thresholds[,c(7,4)]
# rownames(matrix_thresh)<-unlist(unique(matrix_thresh$Sp_short))
# matrix_thresh<-data.matrix(matrix_thresh[,2], rownames.force = TRUE)
# head(matrix_thresh)
# 
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/outliers/thresh_correlation.svg", sep = ""), width = 6, height = 3.5, bg = "white")
# ggplot(thresholds, aes(value)) + geom_histogram(aes(fill=Sp_short), color="white") + theme_bw() + th + scale_fill_d3(palette = "category20")
# dev.off()
# 
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/outliers/median_correlation.svg", sep = ""), width = 4, height = 4, bg = "white")
# unique(thresholds[,c(1,4,7)]) %>% dplyr::group_by(Bug_ID,Sp_short) %>% summarise(median_corr = median(value)) %>%
#   ggplot(aes(x=Sp_short,y=median_corr)) + geom_bar(stat = "identity") + theme_bw() + th + 
#   theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) + scale_x_discrete(name = "Species")
# dev.off()
# 
# ###############################################################################
# # discard outliers for further analysis
# 
# View(total_hclust_outliers)
# 
# NT5001<-total_hclust_outliers %>% filter(Bug_ID == 'NT5001' & dist == "thresh_5")
# head(NT5001)
# nrow(NT5001) #657
# NT5004<-total_hclust_outliers %>% filter(Bug_ID == 'NT5004' & dist == "thresh_4")
# head(NT5004)
# nrow(NT5004) #752
# NT5009<-total_hclust_outliers %>% filter(Bug_ID == 'NT5009' & dist == "thresh_4")
# head(NT5009)
# nrow(NT5009) #1,064
# NT5011<-total_hclust_outliers %>% filter(Bug_ID == 'NT5011' & dist == "thresh_5")
# head(NT5011)
# nrow(NT5011) #469
# NT5019<-total_hclust_outliers %>% filter(Bug_ID == 'NT5019' & dist == "thresh_7")
# head(NT5019)
# nrow(NT5019) #283
# NT5021<-total_hclust_outliers %>% filter(Bug_ID == 'NT5021' & dist == "thresh_3")
# head(NT5021)
# nrow(NT5021) #748
# NT5022<-total_hclust_outliers %>% filter(Bug_ID == 'NT5022' & dist == "thresh_3")
# head(NT5022)
# nrow(NT5022) #854
# NT5025<-total_hclust_outliers %>% filter(Bug_ID == 'NT5025' & dist == "thresh_3")
# head(NT5025)
# nrow(NT5025) #938
# NT5026<-total_hclust_outliers %>% filter(Bug_ID == 'NT5026' & dist == "thresh_3")
# head(NT5026)
# nrow(NT5026) #801
# NT5028<-total_hclust_outliers %>% filter(Bug_ID == 'NT5028' & dist == "thresh_4")
# head(NT5028)
# nrow(NT5028) #594
# NT5054<-total_hclust_outliers %>% filter(Bug_ID == 'NT5054' & dist == "thresh_3")
# head(NT5054)
# nrow(NT5054) #784
# NT5078<-total_hclust_outliers %>% filter(Bug_ID == 'NT5078' & dist == "thresh_4")
# head(NT5078)
# nrow(NT5078) #628
# NT5083<-total_hclust_outliers %>% filter(Bug_ID == 'NT5083' & dist == "thresh_9")
# head(NT5083)
# nrow(NT5083) #277
# 
# hclust_outliers_total<-do.call(rbind, list(NT5001,NT5004,NT5009,NT5011,NT5019,NT5021,NT5022,NT5025,NT5026,NT5028,NT5054,NT5078,NT5083))
# head(hclust_outliers_total)
# 
# hclust_outliers<-hclust_outliers_total %>% dplyr::group_by(Bug_ID,Plate_no) %>% mutate(plate_wells = length(well))
# hclust_outliers<-unique(hclust_outliers[,c(1,7:8)])
# hclust_outliers<-merge(hclust_outliers,bug_desc,by="Bug_ID")
# head(hclust_outliers)
# 
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/outliers/outlier_wells.svg", sep = ""), width = 6, height = 4, bg = "white")
# hclust_outliers %>%
#   ggplot(aes(x=Sp_short,y=outwell_count)) + geom_bar(aes(fill=dist), stat = "identity") + theme_bw() + th + 
#   scale_fill_jama() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(name = " # outlier wells") +
#   geom_text(aes(label=outwell_count), color="black", size = 2.5, nudge_y = 50)
# dev.off()  
# 
# # subtract outlier wells
# head(params)
# nrow(params) #52,608
# View(hclust_outliers_total)
# nrow(hclust_outliers_total) #8,849
# 
# #filter out all except H12 control wells
# hclust_outliers_total<-hclust_outliers_total %>% filter(well != 'H12')
# nrow(hclust_outliers_total) #8,760
# 
# params_out<-anti_join(params,hclust_outliers_total, by=c("Bug_ID","Plate_no","well","Replicate_no"))
# head(params_out)
# nrow(params_out) #43,848
# 
# write.table(params_out, file = "/Users/vinitaperiwal/GrowthCurver/Figures/params_out", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# ##########################################
# #pair-wise wilcoxon test (compares distribution) for identifying differential replicates
# head(SFall)
# nrow(SFall) #52,608
# 
# bugs<-unique(SFall$Bug_ID)
# bugs
# 
# df1<-data.frame(stringsAsFactors = FALSE)
# 
# for(i in 1:length(bugs)){
#   
#   f<-SFall %>% filter(Bug_ID == bugs[i])
#   head(f)
#   
#   f_plates<-unique(f$Plate_no)
#   f_plates
#   
#   for(j in 1:length(f_plates)){
#     
#     plate<-f_plates[j]
#     #possible pairs
#     rep1_rep2<-f %>% filter(Replicate_no %in% c('Replicate1','Replicate2') & Plate_no == plate)
#     rep1_rep3<-f %>% filter(Replicate_no %in% c('Replicate1','Replicate3') & Plate_no == plate)
#     rep2_rep3<-f %>% filter(Replicate_no %in% c('Replicate2','Replicate3') & Plate_no == plate)
#     rep2_rep4<-f %>% filter(Replicate_no %in% c('Replicate2','Replicate4') & Plate_no == plate)
#     rep3_rep4<-f %>% filter(Replicate_no %in% c('Replicate3','Replicate4') & Plate_no == plate)
#     
#     if(nrow(rep1_rep2) > 96){ #each replicate can have max 1,344 rows (14 plates*96 wells)
#       
#       test<-wilcox.test(as.vector(rep1_rep2$normPlAUC) ~ rep1_rep2$Replicate_no)  
#       t.val<-test$statistic
#       p.val<-test$p.value
#       
#       fr<-data.frame(unique(rep1_rep2$Bug_ID),plate,"r1_r2",t.val,p.val)
#       colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
#       df1<-rbind(df1,fr)
#       rownames(df1)<-NULL
#       df1
#     }
#     if(nrow(rep1_rep3) > 96){
#       
#       test<-wilcox.test(as.vector(rep1_rep3$normPlAUC) ~ rep1_rep3$Replicate_no)  
#       t.val<-test$statistic
#       p.val<-test$p.value
#       
#       fr<-data.frame(unique(rep1_rep3$Bug_ID),plate,"r1_r3",t.val,p.val)
#       colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
#       df1<-rbind(df1,fr)
#       rownames(df1)<-NULL
#       df1
#       
#     }
#     if(nrow(rep2_rep3) > 96){
#       
#       test<-wilcox.test(as.vector(rep2_rep3$normPlAUC) ~ rep2_rep3$Replicate_no)  
#       t.val<-test$statistic
#       p.val<-test$p.value
#       
#       fr<-data.frame(unique(rep2_rep3$Bug_ID),plate,"r2_r3",t.val,p.val)
#       colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
#       df1<-rbind(df1,fr)
#       rownames(df1)<-NULL
#       df1
#       
#     }
#     if(nrow(rep2_rep4) > 96){
#       
#       test<-wilcox.test(as.vector(rep2_rep4$normPlAUC) ~ rep2_rep4$Replicate_no)  
#       t.val<-test$statistic
#       p.val<-test$p.value
#       
#       fr<-data.frame(unique(rep2_rep4$Bug_ID),plate,"r2_r4",t.val,p.val)
#       colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
#       df1<-rbind(df1,fr)
#       rownames(df1)<-NULL
#       df1
#       
#     }
#     if(nrow(rep3_rep4) > 96){
#       
#       test<-wilcox.test(as.vector(rep3_rep4$normPlAUC) ~ rep3_rep4$Replicate_no) 
#       test
#       t.val<-test$statistic
#       p.val<-test$p.value
#       
#       fr<-data.frame(unique(rep3_rep4$Bug_ID),plate,"r3_r4",t.val,p.val)
#       colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
#       df1<-rbind(df1,fr)
#       rownames(df1)<-NULL
#       df1
#       
#     }
#   }
#   
# }
# 
# View(df1)
# #bug names
# bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
# head(bug_desc)
# nrow(bug_desc)
# 
# all_bugs_df1<-merge(df1,bug_desc, by="Bug_ID")
# head(all_bugs_df1)
# nrow(all_bugs_df1) #538
# 
# #list significant p-values (problem plates)
# problem_plates<-all_bugs_df1 %>% filter(pval < 0.05)
# head(problem_plates)
# nrow(problem_plates) #6
# 
# #plot significant p-values
# CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/wilcoxon_ranktest_platewise.svg", sep = ""), width = 8, height = 5, bg = "white")
# all_bugs_df1 %>% filter(pval < 0.05) %>% ggplot(aes(x=comp,y=pval)) + geom_point(aes(fill=plate), shape=21, size = 3, position = position_dodge2(width = 0.3)) +
#   theme_bw() + th + scale_fill_d3(palette = "category20") + facet_wrap("Sp_short", scales = "free") +
#   theme(axis.text.x = element_text()) + scale_y_continuous(trans = "exp")
# dev.off()  

###################### boxplot of effect of food compounds alone: 1 - outlier method
# head(SFa)
# 
# food_outliers<-SFa %>% melt() %>% 
#   filter(variable == 'SFa') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
#   mutate(outlier = ifelse(is_outlier(value), as.character(well), "")) %>% dplyr::group_by(Bug_ID,Plate_no,outlier) %>% 
#   mutate(count_of_wells = length(outlier), label = as.character(lapply(count_of_wells, function(x){
#     if(x == 3 | x == 4){
#       'High-conf'
#     }else if(x == 2){
#       'Med-conf'
#     }else if(x == 1){
#       'Low-conf'
#     }
#   })))
# 
# head(food_outliers)
# food_outliers<-merge(bug_desc, food_outliers)
# head(food_outliers)
# nrow(food_outliers) #8,455
# 
# write.table(food_outliers, file = "/Users/vinitaperiwal/GrowthCurver/Figures/food_outliers", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/food_normAUCl.svg", width = 12, height = 6, bg = "white")
# food_outliers %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% filter(!label %in% c('NA','','Low-conf')) %>%
#   ggplot(aes(Replicate_no,value)) + geom_boxplot(aes(fill=Plate_no), outlier.size = 0.7, lwd=0.3) + th +
#   facet_wrap("Species", scales = "free", nrow = 3) + theme(axis.text.x=element_blank(), axis.title.x = element_blank()) + scale_y_continuous(name = "normalized AUC") +
#   scale_color_nejm() + geom_text(aes(label=outlier, color=label), size = 2.5) + scale_fill_jama(alpha = 0.8) 
# dev.off()

######################################################

#pair-wise welch two sample t-test for identifying differential replicates
View(SFall) 
nrow(SFall) #52,128

bugs<-unique(SFall$Bug_ID)
bugs

total_plates<-SFall %>% dplyr::group_by(Bug_ID,Replicate_no) %>% summarise(plate_count = length(unique(Plate_no)))
nrow(total_plates)
View(total_plates)
sum(total_plates$plate_count) #543

df<-data.frame(stringsAsFactors = FALSE)

for(i in 1:length(bugs)){
  
  f<-SFall %>% filter(Bug_ID == bugs[i])
  head(f)
  
  f_plates<-unique(f$Plate_no)
  f_plates
  
  for(j in 1:length(f_plates)){
    
    plate<-f_plates[j]
    #possible pairs
    rep1_rep2<-f %>% filter(Replicate_no %in% c('Replicate1','Replicate2') & Plate_no == plate)
    rep1_rep3<-f %>% filter(Replicate_no %in% c('Replicate1','Replicate3') & Plate_no == plate)
    rep1_rep4<-f %>% filter(Replicate_no %in% c('Replicate1','Replicate4') & Plate_no == plate)
    rep2_rep3<-f %>% filter(Replicate_no %in% c('Replicate2','Replicate3') & Plate_no == plate)
    rep2_rep4<-f %>% filter(Replicate_no %in% c('Replicate2','Replicate4') & Plate_no == plate)
    rep3_rep4<-f %>% filter(Replicate_no %in% c('Replicate3','Replicate4') & Plate_no == plate)
    
    if(nrow(rep1_rep2) > 96){ #each replicate can have max 1,344 rows (14 plates*96 wells)
      
      test<-t.test(as.vector(rep1_rep2$normAUC) ~ rep1_rep2$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep1_rep2$Bug_ID),plate,"r1_r2",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
    }
    if(nrow(rep1_rep3) > 96){
      
      test<-t.test(as.vector(rep1_rep3$normAUC) ~ rep1_rep3$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep1_rep3$Bug_ID),plate,"r1_r3",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
      
    }
    if(nrow(rep1_rep4) > 96){
      
      test<-t.test(as.vector(rep1_rep4$normAUC) ~ rep1_rep4$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep1_rep4$Bug_ID),plate,"r1_r4",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
      
    }
    if(nrow(rep2_rep3) > 96){
      
      test<-t.test(as.vector(rep2_rep3$normAUC) ~ rep2_rep3$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep2_rep3$Bug_ID),plate,"r2_r3",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
      
    }
    if(nrow(rep2_rep4) > 96){
      
      test<-t.test(as.vector(rep2_rep4$normAUC) ~ rep2_rep4$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep2_rep4$Bug_ID),plate,"r2_r4",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
      
    }
    if(nrow(rep3_rep4) > 96){
      
      test<-t.test(as.vector(rep3_rep4$normAUC) ~ rep3_rep4$Replicate_no)  
      t.val<-test$statistic
      p.val<-test$p.value
      
      fr<-data.frame(unique(rep3_rep4$Bug_ID),plate,"r3_r4",t.val,p.val)
      colnames(fr)<-c("Bug_ID","plate","comp","t.val","pval")
      df<-rbind(df,fr)
      rownames(df)<-NULL
      df
      
    }
  }
  
}

nrow(df) #541
View(df)
#bug names
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc)

all_bugs_df<-merge(df,bug_desc, by="Bug_ID")
head(all_bugs_df)
nrow(all_bugs_df) #541

#plot distribution of p values
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/t_test_pvalues.svg", sep = ""), width = 5, height = 3, bg = "white")
ggplot(all_bugs_df, aes(pval)) + geom_histogram(aes(fill=""), color="white") + 
  theme_bw() + th + scale_fill_jama()
dev.off()

#list significant p-values (problem plates)
problem_plates<-all_bugs_df %>% filter(pval < 0.01)
head(problem_plates)
nrow(problem_plates) #66

write.table(problem_plates, file = paste0("/Users/vinitaperiwal/GrowthCurver/problem_plates"), sep = "\t", quote = FALSE, row.names = FALSE)

#plot significant p-values
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/t_test_platewise.svg", sep = ""), width = 8, height = 5, bg = "white")
all_bugs_df %>% filter(pval < 0.01) %>% ggplot(aes(x=comp,y=pval)) + geom_point(aes(fill=plate), shape=21, size = 3, position = position_dodge2(width = 0.3)) +
  theme_bw() + th + scale_fill_d3(palette = "category20") + facet_wrap("Sp_short", scales = "free") +
  theme(axis.text.x = element_text()) + scale_y_continuous(trans = "exp")
dev.off()  

#plot count of problem plates
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/t_test_problemplates.svg", sep = ""), width = 4, height = 3, bg = "white")
problem_plates %>% dplyr::group_by(Bug_ID,Sp_short) %>% summarise(count = length(plate)) %>%
  ggplot(aes(x=Sp_short,y=count)) + geom_bar(aes(fill=""), stat = "identity") +
  theme_bw() + th + scale_fill_jama() + geom_text(aes(label=count), check_overlap = TRUE, nudge_y = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) + scale_y_discrete(name = "# problem plates") +
  scale_x_discrete(name = "Species")
dev.off()


##############################################################################
# correlation after removing problem plates
problem_plates<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/problem_plates", header = TRUE, sep = '\t', stringsAsFactors = FALSE )
head(SFall)
nrow(SFall) #52128

#remove problem plates
head(problem_plates)
nrow(problem_plates) #66
p_plates<-problem_plates[,1:3]
head(p_plates)

p_plates<-separate(p_plates, comp, into = c("repA","repB"), sep = "_")
head(p_plates)

#rename replicates
p_plates[p_plates=="r1"]<-"Replicate1"
p_plates[p_plates=="r2"]<-"Replicate2"
p_plates[p_plates=="r3"]<-"Replicate3"
p_plates[p_plates=="r4"]<-"Replicate4"

nrow(p_plates) #66
head(p_plates) 

#melt
melt_p_plates<-p_plates %>% dplyr::group_by(Bug_ID,plate) %>% melt(measure.vars = c("repA","repB"))
head(melt_p_plates)
nrow(melt_p_plates) #132

colnames(melt_p_plates)[4]<-"Replicate_no"
melt_p_plates<-melt_p_plates[,c(1,2,4)]
head(melt_p_plates)
colnames(melt_p_plates)[2]<-"Plate_no"

control_plates<-melt_p_plates %>% filter(Plate_no %in% c('plate12','plate13'))
View(unique(control_plates)) #14

#distinct plates
View(melt_p_plates %>% distinct()) #109, equals to 10,272 wells

# now subtract melt_p_plates from SFall
SFall_filtered<-anti_join(SFall, melt_p_plates, by=c("Bug_ID","Plate_no","Replicate_no"))
head(SFall_filtered)
nrow(SFall_filtered) #41664

# proceed for correlation with filtered data
mat_melt<-SFall_filtered %>% melt() %>% reshape(idvar = c("Bug_ID","Plate_no","Drug_name","well","Phyla","Species","Sp_short","variable"), timevar = c("Replicate_no"), direction = "wide")
head(mat_melt)  
nrow(mat_melt) #30528
colnames(mat_melt)<-c("Bug_ID","Plate_no","Drug_name","well","Phyla","Species","Sp_short","variable","Rep1","Rep2","Rep3","Rep4")

#Select variable to compute correlation
mat<-mat_melt %>% filter(variable == 'normAUC')
head(mat)
nrow(mat) 

bugs<-unique(mat$Bug_ID)
bugs

correl_df<-data.frame(Bug_ID=character(),
                      Plate_no=character(),
                      comparisons=character(),
                      value=double(),
                      Parameter=character(),
                      stringsAsFactors=FALSE)

for(i in 1:length(bugs)){
  
  f<-mat %>% filter(Bug_ID == bugs[i])
  head(f)
  
  g<-data.frame(f[colSums(!is.na(f)) > 0])
  head(g)
  nrow(g)
  
  #g<-g %>% drop_na()
  
  plates<-unique(g$Plate_no)
  plates
  
  for(j in 1:length(plates)){
    
    c<-g %>% filter(Plate_no == plates[j])
    head(c)
    c<-na.omit(c)
    
    if(nrow(c) > 4){
      
      #print(nrow(c))   
      Bug_ID<-c[1,1]
      Bug_ID
      Plate_no<-plates[j]
      Plate_no
      
      select_numeric_columns<-c %>% select_if(., is.numeric) 
      
      correlations<-correlate(select_numeric_columns, method = "pearson")
      correlations
      corr<-melt(correlations)
      corr<-na.omit(corr)
      corr
      corr$comp<-paste(corr$rowname,corr$variable,sep = "-")
      corr<-corr[,c(4,3)]
      
      combine_df<-cbind(data.frame(Bug_ID,Plate_no), corr)
      print(combine_df)
      correl_df<-rbind(correl_df, data.frame(combine_df))
      
    }
  }
  
}

head(correl_df)

corr_values<-correl_df
head(corr_values)
nrow(corr_values) #784
corr_values<-corr_values %>% filter(!comp %in% c("Rep3-Rep2","Rep4-Rep2","Rep4-Rep3","Rep2-Rep1","Rep3-Rep1","Rep4-Rep1"))
nrow(corr_values) #392

write.table(corr_values, file = paste0("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC_correl_pearson_filter"), sep = "\t", quote = FALSE, row.names = FALSE)

median_corr_bug_wise<-corr_values %>% dplyr::group_by(Bug_ID) %>%
  dplyr::summarise(median_corr = median(value))

median_corr_bug_wise

write.table(median_corr_bug_wise, file = paste0("/Users/vinitaperiwal/GrowthCurver/median_correl_pearson_filtered"), sep = "\t", quote = FALSE, row.names = FALSE)

#bug names
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc)

all_bugs_desc<-merge(median_corr_bug_wise,bug_desc, by="Bug_ID")
head(all_bugs_desc)

all_bugs_corr_desc<-merge(corr_values,bug_desc, by="Bug_ID")
head(all_bugs_corr_desc)

#plots colored histogram of all correlation values
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC_correlation_filtered.svg", sep = ""), width = 6, height = 4, bg = "white")
ggplot(all_bugs_corr_desc, aes(value)) + geom_histogram(aes(fill=Sp_short), color="white") + theme_bw() + th + scale_fill_d3(palette = "category20")
dev.off()

#plots median correlation values of plates of each bug
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/median_correlation_filtered.svg", sep = ""), width = 4, height = 4, bg = "white")
all_bugs_desc %>% dplyr::group_by(Bug_ID,Sp_short) %>% 
  ggplot(aes(x=Sp_short,y=median_corr)) + geom_bar(stat = "identity") + theme_bw() + th + 
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1)) + scale_x_discrete(name = "Species")
dev.off()

