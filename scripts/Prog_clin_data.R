## library
# devtools::install_github("bhklab/PredictioR")
library(stringr)
library(forestplot)
library(data.table)
library(PredictioR)
library(ComplexHeatmap)
library(MultiAssayExperiment)

##################################################################
## set-up directory/path
##################################################################

dirPublic <- "data/public"
dirPrivate <- "data/private"
dir_output <- "result"

############################################
## load all data
############################################
## public/discovery cohorts
publicExpr <- lapply(1:length(list.files(dirPublic)), function(k){
  load(file.path(dirPublic, list.files(dirPublic)[k]))
  dat_icb
})

study.icb <- substr(list.files(dirPublic), 5, nchar(list.files(dirPublic)) - 4)
names(publicExpr) <- study.icb
study.icb.public <- names(publicExpr)

## private/validation cohorts
privateExpr <- lapply(1:length(list.files(dirPrivate)), function(k){
  load(file.path(dirPrivate, list.files(dirPrivate)[k]))
  dat_icb
})

study.icb <- substr(list.files(dirPrivate), 5, nchar(list.files(dirPrivate)) - 4)
names(privateExpr) <- study.icb
study.icb.private <- names(privateExpr)

## merge expression data
expr <- c(publicExpr, privateExpr)
study.icb <- c(study.icb.public, study.icb.private)

######################################################
## Table S1: Descriptive analysis all clinical data
######################################################
# expr <- publicExpr
# expr <- privateExpr

clin <- lapply(1:length(expr), function(k){

  dat.expr <- assay(expr[[k]])
  dat.clin <- data.frame(colData(expr[[k]]))

  sub.dat.clin <- dat.clin[, c("age", "sex", "treatment", "treatmentid", "cancer_type",
               "response", "survival_time_pfs", "event_occurred_pfs",
               "survival_time_os", "event_occurred_os")]
 
 data.frame(study = study.icb[k] , sub.dat.clin)
 
 })

clin <- do.call(rbind, clin)

# cancer types with limitted number of patients were categorized as 'Other'
clin$cancer_type <- ifelse(clin$cancer_type %in% c("Bladder Bone or Soft tissue", "Colorectum", "HNC", "Head and neck",
                                                   "Liver", "Ovary", "Mesothelium", "Prostate", "Small intestine",
                                                   "Vulva"), "Other", clin$cancer_type) # other (21), Prostate (7), Mesothelium (3), HNC (15), OV (13), and Breast (11)

write.csv(clin, file = file.path(dir_output, "data/clin_summary.csv"))
##################################################################
## total samples, cancer types, treatment
##################################################################
df <- lapply(1:length(expr), function(k){
  
  clin <- colData(expr[[k]])
  cancer <- unique(clin$cancer_type)
  
  if(length(cancer) > 1){
    
  df0 <- data.frame(study = study.icb[k],
                    cancer_type = names(table(clin$cancer_type)),
                    n = as.numeric(table(clin$cancer_type)),
                    treatment = unique(clin$treatment))
    
  }else{
    
    df0 <- data.frame(study = study.icb[k],
                      cancer_type = cancer,
                      n = ncol(expr[[k]]),
                      treatment = unique(clin$treatment))
  }
  
  df0
  
})

df_cancer_treatment <- do.call(rbind, df)
df_cancer_treatment <- df_cancer_treatment[!(df_cancer_treatment$cancer_type %in% c("Unknown", "Lymph_node")), ]
df_cancer_treatment$cancer_type <- ifelse(df_cancer_treatment$cancer_type %in% c("HNC", "Head and neck", "Ovary", "Breast","Mesothelium", "Prostate"), "Other", df_cancer_treatment$cancer_type) # other: other (21), Prostate (7) and Mesothelium (3)
df_cancer_treatment_1 <- df_cancer_treatment[df_cancer_treatment$n <= 10, ]
df_cancer_treatment_1$cancer_type <- "Other"
df_cancer_treatment_2 <- df_cancer_treatment[df_cancer_treatment$n > 10, ]
df_cancer_treatment <- rbind(df_cancer_treatment_1, df_cancer_treatment_2)

cancer <- unique(df_cancer_treatment$cancer_type)
df_cancer <- sapply(1:length(cancer), function(k){
  
  sum(df_cancer_treatment[df_cancer_treatment$cancer_type == cancer[k], "n"])
  
})
names(df_cancer) <- cancer

treatment <- unique(df_cancer_treatment$treatment)
df_treatment <- sapply(1:length(treatment), function(k){
  
  sum(df_cancer_treatment[df_cancer_treatment$treatment == treatment[k], "n"])
  
})
names(df_treatment) <- treatment

##########################################
## Heatmap 
##########################################

files <- list.files(file.path(dirPublic))
files.name <- substr(files, 1, nchar(files) - 4) 

df_public <- lapply(1:length(files), function(k){
  
  load(file.path(dirPublic, files[k]))
  
  clin <- data.frame(colData(dat_icb))
  data.frame(study = strsplit(files.name[k], "__")[[1]][1],
             n = ncol(dat_icb),
             sex = ifelse(sum(!is.na(clin$sex)) > 0, "Yes", "No" ),
             age = ifelse(sum(!is.na(clin$age)) > 0, "Yes", "No" ),
             cancer_type = unique(clin$cancer_type),
             treatment = unique(clin$treatment),
             response = ifelse(sum(!is.na(clin$response)) > 0, "Yes", "No" ),
             os = ifelse(sum(!is.na(clin$survival_time_os)) > 0, "Yes", "No" ),
             pfs = ifelse(sum(!is.na(clin$survival_time_pfs)) > 0, "Yes", "No" ),
             access = "public"
  )
  
})

df_public <- do.call(rbind, df_public)

## private data

files <- list.files(file.path(dirPrivate))
files.name <- substr(files, 1, nchar(files) - 4) 

df_private <- lapply(1:length(files), function(k){
  
  load(file.path(dirPrivate, files[k]))
  
  clin <- data.frame(colData(dat_icb))
  cancer <- unique(clin$cancer_type)
  if(length(cancer) == 1 ){
    
    res <- data.frame(study = strsplit(files.name[k], "__")[[1]][1],
                      n = ncol(dat_icb),
                      sex = ifelse(sum(!is.na(clin$sex)) > 0, "Yes", "No" ),
                      age = ifelse(sum(!is.na(clin$age)) > 0, "Yes", "No" ),
                      cancer_type = unique(clin$cancer_type),
                      treatment = unique(clin$treatment),
                      response = ifelse(sum(!is.na(clin$response)) > 0, "Yes", "No" ),
                      os = ifelse(sum(!is.na(clin$survival_time_os)) > 0, "Yes", "No" ),
                      pfs = ifelse(sum(!is.na(clin$survival_time_pfs)) > 0, "Yes", "No" ),
                      access = "private")
    
  }else{
    
    res <- lapply(1:length(cancer), function(i){
      
      clin_cancer <- clin[clin$cancer_type == cancer[i], ]
      res <- data.frame(study = strsplit(files.name[k], "__")[[1]][1],
                        n = nrow(clin_cancer),
                        sex = ifelse(sum(!is.na(clin_cancer$sex)) > 0, "Yes", "No" ),
                        age = ifelse(sum(!is.na(clin_cancer$age)) > 0, "Yes", "No" ),
                        cancer_type = unique(clin_cancer$cancer_type),
                        treatment = unique(clin_cancer$treatment),
                        response = ifelse(sum(!is.na(clin_cancer$response)) > 0, "Yes", "No" ),
                        os = ifelse(sum(!is.na(clin_cancer$survival_time_os)) > 0, "Yes", "No" ),
                        pfs = ifelse(sum(!is.na(clin_cancer$survival_time_pfs)) > 0, "Yes", "No" ),
                        access = "private")
      
      
    })
    
    res <- do.call(rbind, res)
    
  }
  
  res
  
})

df_private <- do.call(rbind, df_private)

df <- rbind(df_public, df_private)
df <- df[!(df$cancer_type %in% c("Lymph_node", "Unknown")), ]
df$cancer_type <- ifelse(df$cancer_type %in% c("Bladder Bone or Soft tissue", "Colorectum", "HNC", "Head and neck",
                                                   "Liver", "Ovary", "Mesothelium", "Prostate", "Small intestine", 
                                                   "Vulva"), "Other", df$cancer_type) # other (21), Prostate (7), Mesothelium (3), HNC (15), OV (13), and Breast (11)
#df.dat <- df[df$n > 15, ]
sum(df$n) 

## Pie chart: cancer types
df.dat <- df
cancer <- unique(df.dat$cancer_type)
df_cancer <- lapply(1:length(cancer), function(k){
  
  dat <- df.dat[df.dat$cancer_type == cancer[k], ]
  data.frame(cancer = cancer[k],
             freq = sum(dat$n))
  
})

df_cancer <- do.call(rbind, df_cancer)
df_cancer$cancer <- ifelse(df_cancer$freq < 100, "Other", df_cancer$cancer)
df_cancer <- df_cancer[order(df_cancer$freq, decreasing = TRUE), ]
df0 <- df_cancer[1:6, ]
df0[df0$cancer == "Other", "freq"] <- sum(df_cancer[df_cancer$cancer == "Other", "freq"])
df_cancer <- df0
df_cancer$perc <- round((df_cancer$freq/ sum(df_cancer$freq)) * 100)
#df <- df[df$freq >= 50, ]

pdf(file = file.path (dir_output, "data/pie_io_cancer.pdf"), width = 4, height = 4)

bp <- ggplot(df_cancer, aes(x="", y=freq, fill=cancer))+
  #geom_bar(width = 1, stat = "identity") +
  geom_col(color = "#4d4d4d") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#2C6AA5FF", "#BFAB25FF", "#4EA699FF", "#8491BEFF", "#64894DFF",
                              "#D65B5AFF")) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position="none",
    legend.text = element_text(size = 6, face="bold"),
    legend.title = element_blank()
  ) + geom_text(aes(label = paste(perc, "%", sep="")),
           position = position_stack(vjust = 0.5),
            color = "#000000",
            size = 5)

bp

dev.off()



## Pie chart: treatment
df.dat <- df
treatment <- unique(df.dat$treatment)
df_treatment <- lapply(1:length(treatment), function(k){
  
  dat <- df.dat[df.dat$treatment == treatment[k], ]
  data.frame(treatment = unique(dat$treatment),
             freq = sum(dat$n))
  
})

df_treatment <- do.call(rbind, df_treatment)
df_treatment <- df_treatment[order(df_treatment$freq, decreasing = TRUE), ]
df_treatment$perc <- round((df_treatment$freq/sum(df_treatment$freq)) * 100)
#df <- df[df$freq >= 50, ]

pdf(file = file.path(dir_output, "data/pie_io_treatment.pdf"), 
     width = 4, height = 4)

bp <- ggplot(df_treatment, aes(x="", y=perc, fill=treatment))+
  #geom_bar(width = 1, stat = "identity") +
  geom_col(color = "#4d4d4d") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c( "#f1e4df","#dab7aa", "#D6A285FF", "#b17863", "#7b5345", "#5b3d33")) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position="none",
    legend.text = element_text(size = 6, face="bold"),
    legend.title = element_blank()
  )  + 
  geom_text(aes(label = paste(perc, "%", sep="")),
           position = position_stack(vjust = 0.5),
            color = "#000000",
           size = 5)

bp

dev.off()

## heatmap 
df.dat <- df
#df.dat <- df[df$n > 15, ]
study <- unique(df.dat[duplicated(df.dat$study), "study"])
df0 <- df.dat[!(df.dat$study %in% study), ]

df1 <- lapply(1:length(study), function(k){
  print(k)
  df.merged <- df.dat[df.dat$study == study[k], ]
  data.frame(study = study[k],
             n = sum(as.numeric(df.merged$n)),
             sex = df.merged$sex[1],
             age = df.merged$age[1],
             cancer_type = df.merged$cancer_type[1],
             treatment = df.merged$treatment[1],
             response = df.merged$response[1],
             os = df.merged$os[1],
             pfs = df.merged$pfs[1],
             access = df.merged$access[1])
})

df1 <- do.call(rbind, df1)
df.dat <- rbind(df1, df0)
df.dat <- df.dat[order(df.dat$access, decreasing = TRUE), ]
df.dat$study <- substr(df.dat$study, 5, nchar(df.dat$study))
df.dat$study[which(df.dat$study == "Powles")] <- "ABACUS"
df.dat$study[which(df.dat$study == "Wolf")] <- "I-SPY"
df.dat$study[which(df.dat$study == "Wolf2")] <- "I-SPY2"
df.dat$study[which(df.dat$study == "Fehrenbacher")] <- "POPLAR"
df.dat$study[which(df.dat$study == "Rittmeyer")] <- "OAK"

df.heatmap <- df.dat[, colnames(df.dat) %in% c("study", "age", "sex", "response", "os", "pfs", "access" )]
df.heatmap$access <- ifelse(df.heatmap$access == "public", 1, 0)
df.heatmap$sex <- ifelse(df.heatmap$sex == "Yes", 1, 0)
df.heatmap$age <- ifelse(df.heatmap$age == "Yes", 1, 0)
df.heatmap$response <- ifelse(df.heatmap$response == "Yes", 1, 0)
df.heatmap$os <- ifelse(df.heatmap$os == "Yes", 1, 0)
df.heatmap$pfs <- ifelse(df.heatmap$pfs == "Yes", 1, 0)
rownames(df.heatmap) <- df.heatmap$study
df.heatmap <- df.heatmap[, -1]


## heatmap
df.heatmap <- t(df.heatmap)

pdf(file=file.path(dir_output, "data/heatmap_data.pdf"),
     width = 8, height = 2)

ht_data <- Heatmap(df.heatmap, name = "status",
                   #top_annotation = ha,
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 8),
                   col =c("#E1E1E1FF", "#7D96AFFF"),
                   border = TRUE,
                   rect_gp = gpar(col = "#647184FF", lwd = 1))

ht_data
#draw(ht_data,
#     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     #merge_legends = TRUE,
#     heatmap_legend_side = "right",
#     annotation_legend_side="bottom")

dev.off()

