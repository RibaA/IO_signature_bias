## library
# devtools::install_github("bhklab/PredictioR")
library(PredictioR)
library(ggplot2)
library(ComplexHeatmap)
library(forestplot)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)

###########################################################
## set up directory
###########################################################
url <- "https://raw.githubusercontent.com/bhklab/SignatureSets/main/data-raw/signature_information.csv"
dir <-  "result/pan-sex/meta"
dir_male <-  "result/male/meta"
dir_female <-  "result/female/meta"
dir_output <- "result"

sig.info <- read.csv(url)

############################################################
## number of samples per sex
############################################################
dir0 <-  "C:/IO_signature_bias/data/pan-sex" 
files <- list.files(dir0)

res <- lapply(1:length(files), function(k){
 
   load(file.path(dir0, files[k]))
  
  id <- strsplit(substr(files[k], 5, nchar(files[k])-4), "__")
  data.frame(study = paste(id[[1]][1], id[[1]][2], id[[1]][3], sep = "_"),
             Female = length(dat_icb$sex[which(dat_icb$sex == "F")]),
             Male = length(dat_icb$sex[which(dat_icb$sex == "M")]))
  
})

res <- do.call(rbind, res)
res$pan <- sapply(1:nrow(res), function(k){
  
  sum(c(res$Female[k], res$Male[k]))
  
})
  
res_sample <- res  
#########################################################
## Pan sex: OS
#########################################################
os <- read.csv(file.path(dir, "meta_pan_os.csv"))
os$FDR <- p.adjust(os$Pval, method = "BH")
os_pan <- os

os <- read.csv(file.path(dir_female, "meta_pan_os.csv"))
os$FDR <- p.adjust(os$Pval, method = "BH")
os_female <- os

os <- read.csv(file.path(dir_male, "meta_pan_os.csv"))
os$FDR <- p.adjust(os$Pval, method = "BH")
os_male <- os

################################################################
## load association pan-sex (OS) 
################################################################
## pan-sex
dir <- "result/pan-sex/assoc/os"
data.files <- list.files(dir)

res.os <- lapply(1:length(data.files), function(k){
  
  res <- read.csv(file.path(dir, data.files[k]))
  res
})

res.os <- do.call(rbind, res.os)
res.os <- res.os[!is.na(res.os$Coef), ]

## CYT_Rooney and TLS_Cabrita signatures
# sig <- c("TLS_Cabrita", "CYT_Rooney")

df <- res.os[res.os$Gene == "CYT_Rooney", ]
df$study <- sapply(1:nrow(df), function(k){
  
  id <-  substr(df$Study[k], 5, nchar(df$Study[k])-4)
  paste(strsplit(id, "__")[[1]][1], strsplit(id, "__")[[1]][2], strsplit(id, "__")[[1]][3], sep="_")

})

df_pan <- df

## female
dir <- "result/female/assoc/os"
data.files <- list.files(dir)

res.os <- lapply(1:length(data.files), function(k){
  res <- read.csv(file.path(dir, data.files[k]))
  res
})

res.os <- do.call(rbind, res.os)
res.os <- res.os[!is.na(res.os$Coef), ]

df <- res.os[res.os$Gene == "CYT_Rooney", ]
df$study <- sapply(1:nrow(df), function(k){
  
  id <-  substr(df$Study[k], 5, nchar(df$Study[k]))
  paste(strsplit(id, "__")[[1]][1], strsplit(id, "__")[[1]][2], strsplit(id, "__")[[1]][3], sep="_")
})

df_female <- df

## male
dir <- "result/male/assoc/os"
data.files <- list.files(dir)

res.os <- lapply(1:length(data.files), function(k){
  
  res <- read.csv(file.path(dir, data.files[k]))
  res
})

res.os <- do.call(rbind, res.os)
res.os <- res.os[!is.na(res.os$Coef), ]

df <- res.os[res.os$Gene == "CYT_Rooney", ]
df$study <- sapply(1:nrow(df), function(k){
  
  id <-  substr(df$Study[k], 5, nchar(df$Study[k]))
  paste(strsplit(id, "__")[[1]][1], strsplit(id, "__")[[1]][2], strsplit(id, "__")[[1]][3], sep="_")
})

df_male <- df

######################################################
## prepare data for forestplot
######################################################
int <- union(union(df_pan$study, df_male$study), df_female$study)

df_pan_upd <- lapply(1:length(int), function(k){
  
  df_pan[df_pan$study == int[k], ]
  
})

df_pan_upd <- do.call(rbind, df_pan_upd)
df_pan_upd$Lower <- df_pan_upd$Coef - (1.96 * df_pan_upd$SE)
df_pan_upd$Upper <- df_pan_upd$Coef + (1.96 * df_pan_upd$SE)
df_pan_upd$N <- sapply(1:nrow(df_pan_upd), function(k){
  
  res_sample[res_sample$study == df_pan_upd$study[k], "pan"]
  
})

df_female_upd <- lapply(1:length(int), function(k){
  
  res <- df_female[df_female$study == int[k], ]
  if(nrow(res) == 0){
    
    res <- data.frame(Outcome = "OS",
                      Gene = "CYT_Rooney",
                      Study = NA,
                      Coef = NA,
                      SE = NA,
                      N = NA, 
                      Pval = NA,
                      Cancer_type = NA,
                      Treatment = NA,
                      FDR = NA,
                      study = int[k])
    
  }
  
  res
  
})

df_female_upd <- do.call(rbind, df_female_upd)
df_female_upd$Lower <- df_female_upd$Coef - (1.96 * df_female_upd$SE)
df_female_upd$Upper <- df_female_upd$Coef + (1.96 * df_female_upd$SE)
df_female_upd$N <- sapply(1:nrow(df_female_upd), function(k){
  
  res_sample[res_sample$study == df_female_upd$study[k], "Female"]
  
})

df_male_upd <- lapply(1:length(int), function(k){
  
  res <- df_male[df_male$study == int[k], ]
  if(nrow(res) == 0){
    
    res <- data.frame(Outcome = "OS",
                      Gene = "CYT_Rooney",
                      Study = NA,
                      Coef = NA,
                      SE = NA,
                      N = NA, 
                      Pval = NA,
                      Cancer_type = NA,
                      Treatment = NA,
                      FDR = NA,
                      study = int[k])
    
  }
  
  res
  
})

df_male_upd <- do.call(rbind, df_male_upd)
df_male_upd$Lower <- df_male_upd$Coef - (1.96 * df_male_upd$SE)
df_male_upd$Upper <- df_male_upd$Coef + (1.96 * df_male_upd$SE)
df_male_upd$N <- sapply(1:nrow(df_male_upd), function(k){
  
  res_sample[res_sample$study == df_male_upd$study[k], "Male"]
  
})

dat <- lapply(1:length(int), function(k){
  
  data.frame(study = int[k],
             N = paste(df_pan_upd$N[k], df_male_upd$N[k], df_female_upd$N[k], sep="/"),
             pansex_Effect = df_pan_upd$Coef[k],
             pansex_Lower = df_pan_upd$Lower[k],
             pansex_Upper = df_pan_upd$Upper[k],
             male_Effect = df_male_upd$Coef[k],
             male_Lower = df_male_upd$Lower[k],
             male_Upper = df_male_upd$Upper[k],
             female_Effect = df_female_upd$Coef[k],
             female_Lower = df_female_upd$Lower[k],
             female_Upper = df_female_upd$Upper[k])
})

dat <- do.call(rbind, dat)

## create format of forestplot
table_text <- rbind(
  c("Study", "N (Pan/M/F)", "Pan", "Male", "Female"),
  cbind(dat$study, 
        formatC(dat$N, format="d"), 
        ifelse(is.na(dat$pansex_Effect), " ", formatC(dat$pansex_Effect, format="f", digits=2)),
        ifelse(is.na(dat$male_Effect), " ", formatC(dat$male_Effect, format="f", digits=2)),
        ifelse(is.na(dat$female_Effect), " ", formatC(dat$female_Effect, format="f", digits=2)))
)

meta_pansex_vals <- as.numeric(os_pan[os_pan$Gene == "CYT_Rooney", c("Coef", "CI_lower", "CI_upper")])
meta_male_vals <- as.numeric(os_male[os_male$Gene == "CYT_Rooney", c("Coef", "CI_lower", "CI_upper")])
meta_female_vals <- as.numeric(os_female[os_female$Gene == "CYT_Rooney", c("Coef", "CI_lower", "CI_upper")])

# prepare effect size matrices
mean_values <- rbind(
  c(NA, NA, NA),  # Header row
  cbind(dat$pansex_Effect, dat$male_Effect, dat$female_Effect)
)
lower_values <- rbind(
  c(NA, NA, NA), 
  cbind(dat$pansex_Lower, dat$male_Lower, dat$female_Lower)
)
upper_values <- rbind(
  c(NA, NA, NA), 
  cbind(dat$pansex_Upper, dat$male_Upper, dat$female_Upper)
)

# Meta-estimate row
mean_values <- rbind(
  mean_values,
  c(meta_pansex_vals[1], meta_male_vals[1], meta_female_vals[1])  
)

# 95% CI: lower 
lower_values <- rbind(
  lower_values,
  c(meta_pansex_vals[2], meta_male_vals[2], meta_female_vals[2])
)

# 95% CI: upper
upper_values <- rbind(
  upper_values,
  c(meta_pansex_vals[3], meta_male_vals[3], meta_female_vals[3])
)

table_text <- rbind(
  table_text,
  c("Meta-Estimate (RE)", " ",  # Empty space for sample size column
    formatC(as.numeric(meta_pansex_vals)[1], format="f", digits=2),
    formatC(as.numeric(meta_male_vals)[1], format="f", digits=2),
    formatC(as.numeric(meta_female_vals)[1], format="f", digits=2))
)

# generate the forest plot
pdf(file=file.path(dir_output, "CYT_Rooney.pdf"), 
    width = 7, height = 5.7, useDingbats = FALSE)

forestplot(
  labeltext = table_text,
  mean = mean_values, 
  lower = lower_values, 
  upper = upper_values, 
  is.summary = c(TRUE, rep(FALSE, nrow(dat))),  
  xlab = "Estimated logHR",
  col = fpColors(box = c("#647D96FF", "#A8554EFF", "#A89E5EFF"), lines = "#181830FF",
                 summary = c("#647D96FF", "#A8554EFF", "#A89E5EFF")),
  txt_gp = fpTxtGp(label = gpar(fontsize = 9), 
                   xlab = gpar(fontsize = 9)),
  boxsize = 0.2,  
  zero = 0,  
  lwd.zero = 2)

dev.off()






