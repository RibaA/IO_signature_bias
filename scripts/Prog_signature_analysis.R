## library
# devtools::install_github("bhklab/PredictioR")
library(PredictioR)
library(survival)
library(GSVA)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
#library(Hmisc)

##################################################################
## set-up directory/path
##################################################################

url <- "https://raw.githubusercontent.com/bhklab/SignatureSets/main/data-raw/signature_information.csv"
dir <- "data/pan-sex"
dir_signature <- "data/SignatureSets"
dir_output <- "result"

############################################
## load signature data
############################################
signature_info <- read.csv(url)
files <- list.files(dir_signature)
signature <- lapply(1:length(files), function(k){
  
  load(file.path(dir_signature, files[k]))
  sig
})

names(signature) <- substr(files, 1, nchar(files)-4)

############################################
## load pan-sex data
############################################
data_files <- list.files(dir)

dat <- lapply(1:length(data_files), function(k){
  
  load(file.path(dir, data_files[k]))
  assay(dat_icb)
  
})

names(dat) <- substr(data_files, 1, nchar(data_files)-4)

##################################################################
## Compute signature score
##################################################################

for(j in 1:length(dat)){

  expr <- dat[[j]]
  study.icb <- names(dat)[j]
  print(study.icb)
  
  geneSig.score <- lapply(1:length(signature), function(i){ 
    
    print(paste(i , names(signature)[i], sep="/"))
    sig_name <- names(signature)[i]
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "GSVA"){
      
      geneSig <- geneSigGSVA(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study.icb)
      
      
      if(sum(!is.na(geneSig)) > 0){
        geneSig <- geneSig[1,]
      }     
      
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Weighted Mean"){
      
      geneSig <- geneSigMean(dat.icb = expr,
                             sig = signature[[i]],
                             sig.name = sig_name,
                             missing.perc = 0.5,
                             const.int = 0.001,
                             n.cutoff = 15,
                             sig.perc = 0.8,
                             study = study.icb)
      
    }
    
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "ssGSEA"){
      
      geneSig <- geneSigssGSEA(dat.icb = expr,
                               sig = signature[[i]],
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = study.icb)
      
      if(sum(!is.na(geneSig)) > 0){
        geneSig <- geneSig[1,]
      }     
      
      
    }
    
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "COX-IS_Bonavita"){
      
      geneSig <- geneSigCOX_IS(dat.icb = expr,
                               sig = signature[[i]],
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = study.icb)
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPS_Charoentong"){
      
      geneSig <- geneSigIPS(dat.icb = expr,
                            sig = signature[[i]],
                            sig.name = sig_name,
                            missing.perc = 0.5,
                            const.int = 0.001,
                            n.cutoff = 15,
                            study = study.icb)
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PredictIO_Bareche"){
      
      geneSig <- geneSigPredictIO(dat.icb = expr,
                                  sig = signature[[i]],
                                  sig.name = sig_name,
                                  missing.perc = 0.5,
                                  const.int = 0.001,
                                  n.cutoff = 15,
                                  sig.perc = 0.8,
                                  study = study.icb)
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPRES_Hugo"){
      
      geneSig <- geneSigIPRES(dat.icb = expr,
                              sig = signature[[i]],
                              sig.name = sig_name,
                              missing.perc = 0.5,
                              const.int = 0.001,
                              n.cutoff = 15,
                              sig.perc = 0.8,
                              study = study.icb)
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "PassON_Du"){
      
      geneSig <- geneSigPassON(dat.icb = expr,
                               sig = signature[[i]],
                               sig.name = sig_name,
                               missing.perc = 0.5,
                               const.int = 0.001,
                               n.cutoff = 15,
                               sig.perc = 0.8,
                               study = study.icb)
      
    }
    
    if(signature_info[signature_info$signature == sig_name, "method"] == "Specific Algorithm" & sig_name == "IPSOV_Shen"){
      
      geneSig <- geneSigIPSOV(dat.icb = expr,
                              sig = signature[[i]],
                              sig.name = sig_name,
                              missing.perc = 0.5,
                              const.int = 0.001,
                              n.cutoff = 15,
                              sig.perc = 0.8,
                              study = study.icb)
      
    }
    
    
    if(sum(!is.na(geneSig)) > 0){
      
      geneSig <- geneSig
      
    }     
    
    if(sum(!is.na(geneSig)) == 0){
      
      geneSig <- rep(NA, ncol(expr))
      
    }
    
    geneSig
    
  })
  
  geneSig.score <- do.call(rbind, geneSig.score)
  rownames(geneSig.score) <- names(signature)
  remove <- which(is.na(rowSums(geneSig.score)))
  if(length(remove) > 0){
    
    geneSig.score <- geneSig.score[-remove, ]
    
  }
  
  save(geneSig.score, file=file.path(dir_output, "pan-sex/sig", paste(study.icb , "sig_score.rda", sep = "_")))
  
}

############################################################
## Pearson Correlation analysis
############################################################

#fit <- rcorr(t(geneSig.score))
#cor.val <- list ("r" = fit$r,
#                 "n" = fit$n,
#                 "p" = fit$P)

#names(cor.val) <- rep(study.icb, 3)
#save(cor.val, file = file.path(dir_output, "pan-sex", paste(study.icb , "sig_pcor.rda", sep = "_")))

