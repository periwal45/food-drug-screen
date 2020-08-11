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

########################### correlation analysis (for wells having normAUC < 1.1)
all_norm<-read.table(file = "/Users/vinitaperiwal/GrowthCurver/Figures/final_normAUCs", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
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
