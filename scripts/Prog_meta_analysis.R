## library
# devtools::install_github("bhklab/PredictioR")
library(PredictioR)
library(meta)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
#library(Hmisc)

###########################################################
## set up directory
###########################################################

dir <- "result/female/assoc"
dir_output <- "result/female/meta"

############################################
## load association data
############################################
data.files <- list.files(file.path(dir, "os"))
res.os <- lapply(1:length(data.files), function(k){
  
  res <- read.csv(file.path(dir, "os", data.files[k]))
  res
})

res.os <- do.call(rbind, res.os)
res.os <- res.os[!is.na(res.os$Coef), ]
write.csv(res.os, file = file.path(dir_output, "res.os.csv"), row.names= FALSE)

## extract signature score association (PFS) data
data.files <- list.files(file.path(dir, "pfs"))
res.pfs <- lapply(1:length(data.files), function(k){
  
  res <- read.csv(file.path(dir, "pfs", data.files[k]))
  res
  
})

res.pfs <- do.call(rbind, res.pfs)
res.pfs <- res.pfs[!is.na(res.pfs$Coef), ]
write.csv(res.pfs, file = file.path(dir_output, "res.pfs.csv"), row.names= FALSE)

## extract signature score association (response) data
data.files <- list.files(file.path(dir, "response"))
res.logreg <- lapply(1:length(data.files), function(k){
  
  res <- read.csv(file.path(dir, "response", data.files[k]))
  res
  
})

res.logreg <- do.call(rbind, res.logreg)
res.logreg <- res.logreg[!is.na(res.logreg$Coef), ]
write.csv(res.logreg, file = file.path(dir_output, "res.logreg.csv"), row.names= FALSE)

################################################################################
################################################################################
## Association meta-analysis: pan-cancer and response
################################################################################
################################################################################
df <- res.logreg
signature <- unique(df$Gene)

AllGeneSig_meta <- lapply(1:length(signature), function(j){
  
  print(j)
  res <- metafun(coef = df[df$Gene == signature[j], "Coef"],
                 se = df[df$Gene == signature[j], "SE"],
                 study  = df[df$Gene == signature[j], "Study"],
                 pval = df[df$Gene == signature[j], "Pval"],
                 n = df[df$Gene == signature[j], "N"],
                 cancer.type = df[df$Gene == signature[j], "Cancer_type"],
                 treatment = df[df$Gene == signature[j], "Treatment"],
                 cancer.spec = FALSE,
                 treatment.spec = FALSE,
                 feature = unique(df[df$Gene == signature[j], "Gene"]))
  
  res$meta_summery
})

AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]
AllGeneSig_meta$FDR <- p.adjust(AllGeneSig_meta$Pval, method = "BH")
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

write.csv(AllGeneSig_meta, file=file.path(dir_output, "meta_pan_logreg.csv"), row.names = FALSE)

################################################################################
################################################################################
## Association meta-analysis: pan-cancer and OS
################################################################################
################################################################################
df <- res.os
signature <- unique(df$Gene)

AllGeneSig_meta <- lapply(1:length(signature), function(j){
  
  print(j)
  res <- metafun(coef = df[df$Gene == signature[j], "Coef"],
                 se = df[df$Gene == signature[j], "SE"],
                 study  = df[df$Gene == signature[j], "Study"],
                 pval = df[df$Gene == signature[j], "Pval"],
                 n = df[df$Gene == signature[j], "N"],
                 cancer.type = df[df$Gene == signature[j], "Cancer_type"],
                 treatment = df[df$Gene == signature[j], "Treatment"],
                 cancer.spec = FALSE,
                 treatment.spec = FALSE,
                 feature = unique(df[df$Gene == signature[j], "Gene"]))
  
  res$meta_summery
})

AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]
AllGeneSig_meta$FDR <- p.adjust(AllGeneSig_meta$Pval, method = "BH")
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

write.csv(AllGeneSig_meta, file=file.path(dir_output, "meta_pan_os.csv"), row.names = FALSE)

################################################################################
################################################################################
## Association meta-analysis: pan-cancer and PFS
################################################################################
################################################################################
df <- res.pfs
signature <- unique(df$Gene)

AllGeneSig_meta <- lapply(1:length(signature), function(j){
  
  print(j)
  res <- metafun(coef = df[df$Gene == signature[j], "Coef"],
                 se = df[df$Gene == signature[j], "SE"],
                 study  = df[df$Gene == signature[j], "Study"],
                 pval = df[df$Gene == signature[j], "Pval"],
                 n = df[df$Gene == signature[j], "N"],
                 cancer.type = df[df$Gene == signature[j], "Cancer_type"],
                 treatment = df[df$Gene == signature[j], "Treatment"],
                 cancer.spec = FALSE,
                 treatment.spec = FALSE,
                 feature = unique(df[df$Gene == signature[j], "Gene"]))
  
  res$meta_summery
})

AllGeneSig_meta <- do.call(rbind, AllGeneSig_meta)
AllGeneSig_meta <- AllGeneSig_meta[!is.na(AllGeneSig_meta$Coef), ]
AllGeneSig_meta$FDR <- p.adjust(AllGeneSig_meta$Pval, method = "BH")
AllGeneSig_meta <- AllGeneSig_meta[order(AllGeneSig_meta$FDR), ]

write.csv(AllGeneSig_meta, file=file.path(dir_output, "meta_pan_pfs.csv"), row.names = FALSE)

