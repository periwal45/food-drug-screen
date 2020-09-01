############################################### Functions ##########################################################

# set global theme for all plots

th<-theme(plot.title = element_text(size = 12, face = "bold"),axis.title=element_text(size=12,color = "black"),
          axis.text.x = element_text(size=10, color = "black"),axis.text.y = element_text(size=10, color = "black"))

#melt dataframe of time and ODs (Used in Preprocess_Raw_Screens.R)
myfun <- function(x){
  melt(x,id.vars = "Time_h", variable.name = "wells", value.name= "OD")
}

# function for creating annotation from file name (Used in GC_processing.R)

create_annot <- function(input_file){
  
  elements<-strsplit(input_file, "_", perl = TRUE)
  ID<-elements[[1]]
  Replicate_no<-elements[[1]][2]
  Plate_no<-elements[[1]][3]
  source<-elements[[1]][4]
  Drug_name<-elements[[1]][5]
  sp<-strsplit(ID, "/", perl = TRUE)
  Bug_ID<-sp[[1]][1]
  
  annotation<-rbind(c(Bug_ID, Replicate_no, Plate_no, source, Drug_name))
  
  return(annotation)
  
}

## fit and plot different distributions (read paper fitdistrplus, 2018)
### function to compute p-value from t-distribution

compute_pval_syn<-function(total_list_bliss){
  
  head(total_list_bliss)
  nrow(total_list_bliss)
  scores<-NULL
  scores<-(total_list_bliss %>% filter(outcome == 'reference'))$bliss_score

  f <- fitStudent(abs(scores))
  f
  
  total_list_bliss$pv <- pstudent(abs(total_list_bliss$bliss_score), f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  
  df<-data.frame(total_list_bliss)
  return(df)

}

compute_pval_bliss_ant<-function(total_list_bliss){
  
  head(total_list_bliss)
  nrow(total_list_bliss)
  scores<-NULL
  scores<-(total_list_bliss %>% filter(outcome == 'reference'))$bliss_score
  
  f <- fitStudent_ant(scores)
  f
  
  total_list_bliss$pv <- 1-pstudent(total_list_bliss$bliss_score, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  
  df<-data.frame(total_list_bliss)
  return(df)
  
}

fitStudent <- function(y) {
  
  f <- NULL
  f <- fitdistrplus::fitdist(y, "student", start=list(nu=30, mean=mean(y), sigma=sd(y)))
  f
}

compute_pval<-function(bug_aucs){
  
  head(bug_aucs)
  nrow(bug_aucs)
  
  normAUCs<-(bug_aucs %>% filter(label == 'reference'))$normAUC
  head(normAUCs)
  
  #f<- fitdistrplus::fitdist(normAUCs, "weibull", method="mge", gof="AD")
  f <- fitStudentOrNormal(normAUCs)
  f
  
  #if (is.na(f$estimate["nu"])) {
  
  #bug_aucs$pv <- pnorm(bug_aucs$normAUC, f$estimate["mean"], f$estimate["sd"])
  #bug_aucs$pv_upper <- pnorm(bug_aucs$normAUC, f$estimate["mean"], f$estimate["sd"], lower.tail = F)
  
  #} else {
  
  bug_aucs$pv <- pstudent(bug_aucs$normAUC, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"])
  #bug_aucs$pv_upper <- pstudent(bug_aucs$normAUC, f$estimate["nu"], f$estimate["mean"], f$estimate["sigma"], lower.tail = F)
  #}
  
  
  df<-data.frame(bug_aucs)
  return(df)
  
}

########


dstudent <- function(x, nu, mean, sigma) dt((x-mean)/sigma, nu, log=F)/sigma
pstudent <- function(q, nu, mean, sigma, lower.tail = TRUE, log.p = FALSE) pt((q-mean)/sigma, nu, lower.tail = lower.tail, log.p=log.p)
qstudent <- function(p, nu, mean, sigma) (qt(p, nu)*sigma+mean)

fitStudentOrNormal <- function(y) {
  f <- NULL
  #try({
   f <- fitdistrplus::fitdist(y, "student", start=list(nu=30, mean=mean(y), sigma=sd(y)))
  #})
  #if (is.null(f) || f$estimate[1] > 100) {
    #f <- fitdistrplus::fitdist(y, "norm")
  #}
  f
}


# n<-length(sample_wells$well)
# n
# 
# total_wells<-rbind(ref_wells,sample_wells)
# nrow(total_wells)
# head(total_wells)
# 
# out_df<-data.frame(stringsAsFactors = FALSE)
# 
# for(i in 1:length(total_wells$well)){
#   
#   Bug_ID<-total_wells[i,1]
#   Replicate_no<-total_wells[i,2]
#   Plate_no<-total_wells[i,3]
#   Drug_name<-total_wells[i,4]
#   well<-total_wells[i,5]
#   auc_l<-total_wells[i,6]
#   normAUC<-total_wells[i,7]
#   Phyla<-total_wells[i,8]
#   Species<-total_wells[i,9]
#   Sp_short<-total_wells[i,10]
#   label<-total_wells[i,11]
#   
#   t.val<-(as.numeric(normAUC)-mean_ref_normAUC)/(sd_sample_normAUC/sqrt(n))
#   t.val
#   p.val<-pt(-abs(t.val),df=n-1)
#   p.val
#   
#   df<-data.frame(Bug_ID,Replicate_no,Plate_no,Drug_name,well,auc_l,normAUC,t.val,p.val,Phyla,Species,Sp_short,label)
#   out_df<-rbind(out_df,df)
#   
# }
# 
# return(out_df)
# }