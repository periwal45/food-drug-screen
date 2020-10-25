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

#bug annotations
bug_desc<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/Bug_desc", sep = "\t", header = TRUE)
head(bug_desc)
nrow(bug_desc) #23


########## This script contains main analysis code
#' 1. Reads fitted parameters
#' 2. Data preprocessing (removing missing and variable dose data)
#' 3. Normalization
#' 4. Correlation analysis of replicates
#' 5. Hit calling - food compounds (fit t-distribution, computes p - value, combine p-values, correct for FDR)
#' 6. Bliss interactions for combination plates

################ Data preprocessing (filter out incomplete and variable data)
#read parameters
params<-data.frame(read.table("/Users/vinitaperiwal/GrowthCurver/GC_fit_params", header = TRUE, sep = '\t', stringsAsFactors = FALSE))

params<-params %>% filter(Plate_no != 'plate16' & !(Bug_ID %in% c('NT5003','NT5048')))
#View(params)
nrow(params) #47,232

#select specific columns
params<-params[,c(1:5,19)]
head(params)

#removed varied dose plates
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

################ Data normalization
#' normalizing aucs of food only plates by plate median
#' defining reference wells in food compound plates (bcos there are only two control wells)
#' normalize other plates by mean auc of reference wells
#' 
SFall<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/SFall", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
head(SFall)
nrow(SFall) #46,752

#normalizing by plate median only for food plates
normAUCl_food<-SFall %>% filter(Drug_name == 'control') %>% dplyr::group_by(Bug_ID,Replicate_no,Plate_no) %>% 
  mutate(normAUC = round((auc_l/median(auc_l)), digits = 3))

View(normAUCl_food)
nrow(normAUCl_food) #7776

#plot median normalized AUCs
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
write.table(ref, file = "/Users/vinitaperiwal/GrowthCurver/ref_all_wells", sep = "\t", quote = FALSE, row.names = FALSE)

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

## plot of all normalized plates
CairoSVG(file=paste("/Users/vinitaperiwal/GrowthCurver/Figures/normAUC.svg", sep = ""), width = 8, height = 4.6, bg = "white")
all_norm %>% dplyr::group_by(Bug_ID,Replicate_no) %>% filter(normAUC < 2) %>%
  ggplot(aes(x=Sp_short,y=normAUC)) + geom_boxplot(aes(fill = Plate_no), outlier.size = 0.01, lwd=0.15) +
  th + scale_fill_d3(palette = "category20") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#################### correlation analysis (for wells having normAUC < 1.1)
# present in Correlation_analysis.R

#################### Hit calling (food compounds)

###### fit t-distribution and compute p-value per bug plate wise (all replicates pooled)
ref_all_wells<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/ref_all_wells", sep = "\t", header = TRUE)

# label reference wells
ref_all_wells<-ref_allreps %>% filter(count >= 3)
head(ref_all_wells)
nrow(ref_all_wells) #7372

ref_all_wells<-ref_all_wells[,c(1:10)]
ref_all_wells["label"]<-"reference"
head(ref_all_wells)
nrow(ref_all_wells) #7372
ncol(ref_all_wells) #11

#read all normalized aucs
final_normAUCs<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/final_normAUCs", sep = "\t", header = TRUE)
head(final_normAUCs)
nrow(final_normAUCs) #46,752

data_foodplates<-final_normAUCs %>% filter(Drug_name == 'control')
nrow(data_foodplates) #7776
head(data_foodplates)

# label sample wells
samples_data_foodplates<-anti_join(data_foodplates, ref_all_wells, by=c("Bug_ID","Replicate_no","Plate_no","Drug_name","well","auc_l","normAUC","Phyla","Species","Sp_short"))
nrow(samples_data_foodplates) #404
head(samples_data_foodplates)

samples_data_foodplates["label"]<-"sample"
head(samples_data_foodplates)
ncol(samples_data_foodplates)

# combine all reference and sample wells
total_list<-rbind(data.frame(ref_all_wells),data.frame(samples_data_foodplates))
head(total_list)
nrow(total_list) #7776
View(total_list)
ncol(total_list) #11

#### computing p-values from student distribution

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

# plotting p-value distribtuion

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

############## Bliss interactions

final_normAUCs<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/final_normAUCs", sep = "\t", header = TRUE)
head(final_normAUCs)
nrow(final_normAUCs) #46,752

# separate normAUCs of food compounds alone
SFa<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well != 'H12' & Plate_no %in% c("plate12","plate13")) %>%
  mutate(SFa = normAUC)
head(SFa)  
nrow(SFa) #7695

# separate normAUCs of drugs alone
SFq<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well == 'H12' & !(Plate_no %in% c("plate12","plate13","plate16"))) %>%
  mutate(SFq = normAUC)
head(SFq)
nrow(SFq) #406

###################### boxplot of effect of HTD drugs alone (compared with bug alone)
SF_drugs_bugs<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well == 'H12' & Plate_no != "plate16") %>%
  mutate(SFq = normAUC)
head(SF_drugs_bugs)
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

# separate normAUCs of combination wells
SFaq<-final_normAUCs %>% dplyr::group_by(Bug_ID) %>% filter(well != 'H12' & Drug_name != 'control') %>%
  mutate(SFaq = normAUC)
head(SFaq)
nrow(SFaq) #38,570

##### start joining SFaq and SFa and SFq

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

# calculate bliss scores
bliss<-merge_SFaq_SFa_SFq %>% mutate(bliss_ex = SFa*SFq, bliss_score = SFaq-bliss_ex, label = ifelse(SFaq < bliss_ex, "synergy","antagonism"))
#bliss<-merge_SFaq_SFa_SFq %>% mutate(bliss_ex = SFa*SFq, bliss_score = SFaq-bliss_ex)
head(bliss)

write.table(bliss, file = "/Users/vinitaperiwal/GrowthCurver/Figures/bliss_scores", sep = "\t", quote = FALSE, row.names = FALSE)

############# reading food annotations to compound wells #########################
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

short_blissannot<-bliss_annot[,c(1:8,11,15,19,20,21,22,24)]
View(short_blissannot)
nrow(short_blissannot) #38,570

A<-short_blissannot %>% dplyr::group_by(Bug_ID,Plate_no,Drug_name,well) %>% 
  mutate(lSFa = log(SFa),lSFq = log(SFq), lSFaq = log(SFaq))

nrow(A) #38,570

A_filtered<-A %>% filter(lSFa != "-Inf" & lSFq != "-Inf" & lSFaq != "-Inf") %>% 
  dplyr::group_by(Bug_ID,Plate_no,Drug_name,well) %>% 
  mutate(n1 = length(SFa), y1=mean(lSFa), s1=var(lSFa)*(n1-1),
         n2 = length(SFq), y2=mean(lSFq), s2=var(lSFq)*(n2-1),
         n3 = length(SFaq), y3=mean(lSFaq), s3=var(lSFaq)*(n3-1),
         sy=s1+s2+s3, dft=n1+n2+n3-3, denf=1/n1+1/n2+1/n3,
         lbliss=y1+y2-y3,
         expbliss=exp(y1+y2-y3),
         syn_percent=(expbliss-1)*100,
         tss=(y1+y2-y3)/sqrt(sum(sy)/dft)/denf,
         pv=2*(1-pt(abs(tss),df=dft)),
         pvP=1-pt(tss,df=dft)) #n - no of observations (replicates), sy - total sum of squares, dft,denf - df, tss - t-statistic, pv - p-value for Bliss independence hypothesis, pvP - One-sided p-value

nrow(A_filtered) #38,473
View(A_filtered)

B<-A_filtered %>% filter(pv<0.05)
nrow(B) #223

View(B)


B %>% 
  ggplot(aes(x=bliss_score,y=Product.name)) + geom_boxplot(aes(fill=Sp_short), lwd=0.3) + theme_bw() +
  th+ theme(axis.text.x = element_text(angle = 90,hjust = 1), axis.ticks = element_blank(), panel.grid = element_blank()) + 
  scale_fill_d3() + scale_y_discrete(name = "Food compound") + facet_grid(~Drug_name, scales = "free")


