## library
# devtools::install_github("bhklab/PredictioR")
library(stringr)
library(paletteer)
library(forestplot)
library(data.table)
library(PredictioR)
library(ggplot2)
library(MultiAssayExperiment)

##################################################################
## set-up directory/path
##################################################################

url <- "https://raw.githubusercontent.com/bhklab/SignatureSets/main/data-raw/signature_information.csv"
dir_output <- "result"

##################################################################
## signature data: association type
##################################################################
sig <- read.csv(url)
assoc <- unique(sig$association)
df_assoc <- lapply(1:length(assoc), function(k){ 
  data.frame(assoc = assoc[k],
             freq = nrow(sig[sig$association == assoc[k], ]))
})

df_assoc <- do.call(rbind, df_assoc)

## Pie chart: association
df_assoc$perc <- round((df_assoc$freq / sum(df_assoc$freq)) * 100)

pdf(file = file.path(dir_output, "signature/pie_sig_assoc.pdf"), width = 4, height = 4)

bp <- ggplot(df_assoc, aes(x="", y=freq, fill=assoc))+
  #geom_bar(width = 1, stat = "identity") +
  geom_col(color = "#4d4d4d") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#4B5A4BFF", "#96873CFF")) +
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
            size = 7)

bp

dev.off()


## signature data: method
group <- unique(sig$method)

df_method <- lapply(1:length(group), function(k){
  print(k)
  data.frame(method = group[k],
             freq = nrow(sig[sig$method == group[k], ]))
})

df_method <- do.call(rbind, df_method)

## Pie chart: association
df_method$perc <- round((df_method$freq / sum(df_method$freq)) * 100)

pdf(file = file.path (dir_output, "signature/pie_sig_method.pdf"), width = 4, height = 4)

bp <- ggplot(df_method, aes(x="", y=freq, fill=method))+
  #geom_bar(width = 1, stat = "identity") +
  geom_col(color = "#4d4d4d") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#B49696FF", "#E3CA97FF", "#D29678FF", "#96A5A5FF")) +
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
                size = 7)

bp

dev.off()
