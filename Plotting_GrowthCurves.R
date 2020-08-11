setwd("/Users/vinitaperiwal/GrowthCurver/")
library(dplyr)
library(tibble)
library(reshape2)
library(ggfortify)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(lattice)
library(gridExtra)
library(ggsci) 
library(Cairo) 
library(ggforce) #facet wrapping
library(ggpubr) #reg equation

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

##################################### plot bug only wells ########################################################
control_h12<-fits %>% filter(Drug_name == 'control' & (Plate_no %in% c("plate12","plate13")) & col_name == 'H12')
head(control_h12)
nrow(control_h12) #2,162

#bug names
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc)

control_h12_desc<-merge(control_h12,bug_desc, by="Bug_ID")
head(control_h12_desc)

CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/H12_control_raw_fitted_plate12.svg", width = 12, height = 4, bg = "white")
facet1<-control_h12_desc %>% filter(Plate_no == "plate12") %>% dplyr::group_by(Bug_ID,Replicate_no)  %>%
  ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=Replicate_no)) + 
  geom_point(aes(y=OD,color=Replicate_no), size=0.01) + scale_color_nejm() +
  scale_x_continuous(name = "Time (Hours)", limits = c(0,24)) + scale_y_continuous(name="bg corrected OD",) +
  th + theme_bw(base_rect_size = 0.1) + theme(legend.position="bottom") + facet_wrap("Sp_short", scales = "free", ncol = 8)
print(facet1)
dev.off()
CairoSVG(file="/Users/vinitaperiwal/GrowthCurver/Figures/H12_control_raw_fitted_plate13.svg", width = 12, height = 4, bg = "white")
facet2<-control_h12_desc %>% filter(Plate_no == "plate13") %>% dplyr::group_by(Bug_ID,Replicate_no)  %>%
  ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=Replicate_no)) + 
  geom_point(aes(y=OD,color=Replicate_no), size=0.01) + scale_color_nejm() +
  scale_x_continuous(name = "Time (Hours)", limits = c(0,24)) + scale_y_continuous(name="bg corrected OD",) +
  th + theme_bw(base_rect_size = 0.1) + theme(legend.position="bottom") + facet_wrap("Sp_short", scales = "free", ncol = 8)
print(facet2)
dev.off()


#################################### plot all replicates
bugs<-unique(fits$Bug_ID)
bugs

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "control")
  head(f)
  
  CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i], "_control_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "promazine")
  head(f)
  
  CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i], "_promazine_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "loxapine")
  head(f)
  
  CairoSVG( file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i],"_loxapine_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "omeprazole")
  head(f)
  
  CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i], "_omeprazole_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "felodipine")
  head(f)
  
  CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i], "_felodipine_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "5fu")
  head(f)
  
  CairoSVG( file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i],"_5fu_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

for(i in 1:length(bugs)){
  
  if(bugs[i] %in% c('NT5003','NT5022','NT5048')){
    
    f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "duloxetine")
    head(f)
    
    CairoSVG( file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i],"_duloxetine_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
    
    p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
      ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
      geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
      ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
      th + theme_minimal()
    
    print(p)
    dev.off()
  }
}

# for(i in 1:length(bugs)){
#   
#   if(bugs[i] %in% c('NT5003','NT5022','NT5048')){
#   
#   f<-fits %>% filter(Bug_ID == bugs[i], Drug_name == "drug")
#   head(f)
#   
#   CairoSVG( file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i],"_drug_all_replicates.svg", sep = ""), width = 20, height = 20, bg = "white")  
#   
#   p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
#     ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
#     geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
#     ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD", limits = c(0.0,1.0)) + 
#     th + theme_minimal()
#   
#   print(p)
#   dev.off()
#   }
# }

#Only H12 wells of bug

for(i in 1:length(bugs)){
  
  f<-fits %>% filter(Bug_ID == bugs[i], col_name == "H12", Drug_name == 'control')
  head(f)
  
  CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/", bugs[i], "_controlH12_replicates.svg", sep = ""), width = 6, height = 4, bg = "white")
  
  p<-f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
    ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
    geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(col_name ~ Plate_no) + scale_color_nejm() + 
    ggtitle(bugs[i]) + scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
    th + theme_minimal()
  
  print(p)
  dev.off()
  
}

#individual well

f<-fits %>% filter(Bug_ID == 'NT5019', col_name %in% c('H12','A5'), Plate_no %in% c('plate12','plate2','plate8'), Drug_name %in% c('control','promazine','felodipine'))
head(f)

CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/pcop_A5_syn.svg", sep = ""), width = 6, height = 4, bg = "white")
f %>% dplyr::group_by(Bug_ID,Plate_no,col_name) %>%
  ggplot(aes(x=time, y=fitted,group=Replicate_no)) + geom_line(aes(color=Replicate_no)) + 
  geom_point(aes(y=OD,color=Replicate_no), size=0.001) + facet_wrap(c("col_name","Drug_name")) + scale_color_nejm() + 
  scale_x_continuous(name = "Time (Hours)", limits = c(0,25)) + scale_y_continuous(name="bg corrected OD") + 
  th + theme_minimal()
dev.off()

