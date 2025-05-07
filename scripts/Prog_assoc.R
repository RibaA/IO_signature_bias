## library
# devtools::install_github("bhklab/PredictioR")
library(PredictioR)
library(survival)
library(GSVA)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
#library(Hmisc)

###########################################################
## set up directory
###########################################################

dir <- "data/female"
dir_sig <- "result/female/sig"
dir_output <- "result/female/assoc"

#########################################################
## Association with OS
#########################################################

for(i in 1:length(list.files(dir))){
  
  load(file.path(dir, list.files(dir)[i]))
  expr <- dat_icb
  study.icb <- substr(list.files(dir)[i], 1, nchar(list.files(dir)[i])-4) # change to  -7 for sex=specific analysis
  cancer_type <- strsplit(study.icb, "__")[[1]][2]
  treatment_type <- strsplit(study.icb, "__")[[1]][3]
    
  load(file.path(dir_sig, list.files(dir_sig)[i]))
  
  res.os <- lapply(1:nrow(geneSig.score), function(k){
    
    res <- geneSigSurvCont(dat.icb = expr,
                           geneSig = geneSig.score[k,],
                           time.censor = 36,
                           n.cutoff = 25,
                           study =  paste(study.icb, cancer_type, treatment_type, sep="__"),
                           surv.outcome = "OS",
                           sig.name = rownames(geneSig.score)[k],
                           cancer.type = cancer_type,
                           treatment = treatment_type)
    
    res
    
  })
  
  res.os <- do.call(rbind, res.os)
  res.os$FDR <- p.adjust(res.os$Pval, method="BH")
  
  
  write.csv(res.os, file = file.path(dir_output, "os", paste(study.icb , "sig_os.csv", sep = "_")), row.names = FALSE)
  
}

#########################################################
## Association with PFS
#########################################################

for(i in 1:length(list.files(dir))){
  
  load(file.path(dir, list.files(dir)[i]))
  expr <- dat_icb
  study.icb <- substr(list.files(dir)[i], 1, nchar(list.files(dir)[i])-4) # change to 7 for sex-specific analysis
  cancer_type <- strsplit(study.icb, "__")[[1]][2]
  treatment_type <- strsplit(study.icb, "__")[[1]][3]
  
  load(file.path(dir_sig, list.files(dir_sig)[i]))
    
  res.pfs <- lapply(1:nrow(geneSig.score), function(k){
    
    res <- geneSigSurvCont(dat.icb = expr,
                           geneSig = geneSig.score[k,],
                           time.censor = 24,
                           n.cutoff = 25,
                           study =  paste(study.icb, cancer_type, treatment_type, sep="__"),
                           surv.outcome = "PFS",
                           sig.name = rownames(geneSig.score)[k],
                           cancer.type = cancer_type,
                           treatment = treatment_type)
    
    res
    
  })
  
  res.pfs <- do.call(rbind, res.pfs)
  res.pfs$FDR <- p.adjust(res.pfs$Pval, method="BH")
  
  
  write.csv(res.pfs, file = file.path(dir_output, "pfs", paste(study.icb , "sig_pfs.csv", sep = "_")), row.names = FALSE)
  
}
#########################################################
## Association with response
#########################################################

for(i in 1:length(list.files(dir))){
  
  load(file.path(dir, list.files(dir)[i]))
  expr <- dat_icb
  study.icb <- substr(list.files(dir)[i], 1, nchar(list.files(dir)[i])-4) # change to 7 for sex-specific analysis
  cancer_type <- strsplit(study.icb, "__")[[1]][2]
  treatment_type <- strsplit(study.icb, "__")[[1]][3]
  
  load(file.path(dir_sig, list.files(dir_sig)[i]))
  
  res.logreg <- lapply(1:nrow(geneSig.score), function(k){
  
    res <- geneSigLogReg(dat.icb = expr,
                       geneSig = geneSig.score[k,],
                       n.cutoff = 25,
                       study =  paste(study.icb, cancer_type, treatment_type, sep="__"),
                       sig.name = rownames(geneSig.score)[k],
                       n0.cutoff = 3, 
                       n1.cutoff = 3,
                       cancer.type = cancer_type,
                       treatment = treatment_type)
  
  res
  
})

res.logreg <- do.call(rbind, res.logreg)
res.logreg$FDR <- p.adjust(res.logreg$Pval, method="BH")


write.csv(res.logreg, file = file.path(dir_output, "response", paste(study.icb , "sig_logreg.csv", sep = "_")), row.names = FALSE)

}
