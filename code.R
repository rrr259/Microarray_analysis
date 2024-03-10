#changing the directory by choosing working directory in session (directory with all tables needed)
# making defensive R
#function 1 to check if files are there 
check_file_existence <- function(file_paths) {
  if (!all(sapply(file_paths, file.exists))) {
    stop("One or more required files are missing.")
  }
}
#function 2 to check if LongName column in gene annotation
check_long_name_column <- function(gene_annotation) {
  stopifnot(
    is.data.frame(gene_annotation),
    "LongName" %in% colnames(gene_annotation)
  )
}
#function 3 to check if all samples are in data_all
check_samples_columns <- function(data_all) {
  stopifnot(
    is.data.frame(data_all),
    "A" %in% colnames(data_all),
    "B" %in% colnames(data_all),
    "C" %in% colnames(data_all),
    "D" %in% colnames(data_all),
    "E" %in% colnames(data_all),
    "F" %in% colnames(data_all),
    "G" %in% colnames(data_all),
    "H" %in% colnames(data_all),
    "I" %in% colnames(data_all),
    "J" %in% colnames(data_all),
    "K" %in% colnames(data_all),
    "L" %in% colnames(data_all)
  )
}

#loading libraries 
library(readr)
library(dplyr)
library(pheatmap)
library(viridis)

#doing check if all files are present 
files_to_check <- c("data_all.csv", "gene_annotation.csv", "sample_annotation.csv", "genelist_76.txt")
check_file_existence(files_to_check)
#getting tables 
data_all <- read_csv("data_all.csv")
gene_annotation <- read_csv("gene_annotation.csv")
sample_annotation <- read_csv("sample_annotation.csv")
#checking if all information needed is in tables 
check_samples_columns(data_all)
check_long_name_column(gene_annotation)
#view tables 
View(data_all)
View(gene_annotation)
View(sample_annotation)

#get genelist.txt convert it to the csv file 
genes <- read.table("genelist_76.txt", header=TRUE)
View(genes)
#rename the column in the data_all, so it's a gene number 
colnames(data_all)[1] <- "gene_no"
#as gene_no corresponds to x in genes, make genes x vector 
genes_new <- genes$x
#take the data from data_all that corresponds to the genelist given by using %in%
df_new <- data_all[data_all$gene_no %in% genes_new, ] 
#convert numbers to log
df_new_log <- log2(df_new[2:13]+1)
#merge tables 
df_final <- cbind(df_new[,1, drop=FALSE], df_new_log)
#merging according two data frames as correspond to each other 
gene_ann <- merge(df_final, gene_annotation, by.x = "gene_no", by.y = "Gene")
#getting rid of column that contains same information
df_ann <- gene_ann[, -14]

#moving type and longname in the beginning of the table
df_ann_new <- df_ann %>% relocate("Type", .after = "gene_no")
df_ann_final <- df_ann_new %>% relocate("LongName", .after = "gene_no")
#removing unwanted column 
df_ann_final_1 <- df_ann_final[, -1]
#renaming column 
colnames(df_ann_final_1)[1] <- "Gene_name"
#getting new column with merging names 
df_ann_final_1$d <- paste(df_ann_final_1$Gene_name, df_ann_final_1$Type, sep="_")
#removing those columns (as have one that contains information)
df_ann_final_1 <- df_ann_final_1[, -c(1,2)]
#making this column a rownames for the  heatmaps 
row.names(df_ann_final_1) <- df_ann_final_1$d
#removing unwanted column 
df_ann_final_1 <- df_ann_final_1[,-13]

#creating new data frame for annotations using sample_annotation column 
sample_annotation <- sample_annotation[, -1]
data_annotation <- data.frame(row.names = sample_annotation$SampleName, sample_annotation$TreatmentGroup)
colnames(data_annotation)[1] <- "Treatment"

#getting unique values 
treatments <- unique(data_annotation$Treatment)
#calculating lenght 
num_treatments <- length(treatments)
#using viridis to assign specific colours
treatment_colors <- viridis(num_treatments, option = "rocket")
#making a list 
annotation_colors <- list(Treatment = treatment_colors[data_annotation$Treatment])

#creating heatmaps 
#heatmap1
pheatmap(
  df_ann_final_1,
  scale='row',
  annotation_col= data_annotation,
  annotation_colors = annotation_colors,
  cluster_col=FALSE,
  main = 'Heatmap #1',
  legend = FALSE,
  fontsize_row = 6,
  fontsize_col = 12,
  color = viridis(100, option = "plasma"),
  annotation_legend=TRUE
)


#heatmap2 
pheatmap(
  df_ann_final_1,
  scale='row',
  annotation_col= data_annotation,
  annotation_colors = annotation_colors,
  cluster_col=TRUE,
  main = 'Heatmap #2',
  legend = FALSE,
  fontsize_row = 6,
  fontsize_col = 12, 
  color = viridis(100, option= "viridis")
)








