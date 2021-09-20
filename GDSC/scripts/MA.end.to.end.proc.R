# Sharvari Gujja
# Script to download and process microarray data for downstream modeling
# Usage: Rscript <script> <drugname> <outputDir>
# Example: Rscript MA.end.to.end.proc.R Gemcitabine GDSC_Gemcitabine_100cuts/
# Output: Generates 100 cuts of train and test data matrices with normalized gene counts for sig. DE genes
# Outputs:
# Output dir: 1. GDSC_<drugname>_100cuts
  # 1. *affy_limma_train*.csv - DE output on each training set
  # 2. *Affy_DE_TrainingSet_new* -normalized count matrix with id,label,DE sig genes as rows and samples as columns
  # 3. *Affy_boruta_TrainingSet_new* - normalized count matrix with id,label,sig genes from Boruta feature selection as rows and samples as columns
# Output PCA plot of the calibrated, summarized data: GDSC_<drugname>e_PCA.pdf

library(Biobase)
library(oligoClasses)
# 
# #Annotation and data import packages
library(ArrayExpress)
library(hgu219.db)
# 
# #Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
# 
# #Analysis and statistics packages
library(limma)
# library(topGO)
# library(ReactomePA)
# library(clusterProfiler)
# 
library(ggplot2)
library(openxlsx)
library(dplyr)
library(caret)
library(plyr)
library(Boruta)

rm(list = ls())
gc()

args <- commandArgs(TRUE)

# select drug name
prefix <- "Gemcitabine" #args[1] 
# output dir
out <- "test" #args[2] 
if(!dir.exists(out)) { dir.create(out) }

## list for boruta_signif
boruta_signif_list = list()
## list for DE sig genes
tlist_list = list()

#### Note: Download RAW CEL files from GDSC to location below after cloning the repo
# read the raw data
raw_data_dir <- "../data/raw/cel_files/"

SDRF <- read.delim("../data/E-MTAB-3610.sdrf.txt")

# Binarized IC50 file
clin <- as.data.frame(readxl::read_excel("../data/GDSC_binarizedIC50.xlsx"))

# subset for specific drug
clin.sub <- clin[c("Screened Compounds:",prefix)]
names(clin.sub)[2] <- "drug.response"
# merge w/ SDRF
SDRF.clin <- merge(SDRF, clin.sub,by.x="Characteristics.cell.line.", by.y="Screened Compounds:")
# remove missing values for drug response
SDRF.clin <- SDRF.clin[!is.na(SDRF.clin$drug.response), ]
print(dim(SDRF.clin))
print(table(SDRF.clin$drug.response))

rownames(SDRF.clin) <- SDRF.clin$Array.Data.File
SDRF.clin <- AnnotatedDataFrame(SDRF.clin)

# import the files and chip annotation package
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                       SDRF.clin$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF.clin)
stopifnot(validObject(raw_data))

print(dim(raw_data))

print(head(Biobase::pData(raw_data)))
print(Biobase::exprs(raw_data)[1:5, 1:5])


#RMA calibration of the data
data_eset_norm <- oligo::rma(raw_data) #, target = "core")

# PCA analysis
exp_data <- Biobase::exprs(data_eset_norm)
PCA <- prcomp(t(exp_data), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Phenotype = 
                       Biobase::pData(data_eset_norm)$drug.response)

#, colour = Phenotype
p <- ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Phenotype)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) # 
ggsave( paste0("GDSC_",prefix,"_PCA.pdf"), plot = p, device = 'pdf', width = 8, height = 7)

# filter out lowly expressed probes - “soft” intensity based filtering
data_medians <- rowMedians(Biobase::exprs(data_eset_norm))

hist_res <- hist(data_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")


man_threshold <- 3  #change
abline(v = man_threshold, col = "coral4", lwd = 2)

# Transcripts that do not have intensities larger than the threshold in at least as many arrays as the smallest experimental group are excluded.
no_of_samples <- 
  table(pData(data_eset_norm)$drug.response)

samples_cutoff <- min(no_of_samples)

#summarize the results and get an overview over how many genes are filtered out
idx_man_threshold <- apply(Biobase::exprs(data_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
print(table(idx_man_threshold))

data_manfiltered <- subset(data_eset_norm, idx_man_threshold)


# Annotation of the transcript clusters
anno_data <- AnnotationDbi::select(hgu219.db,
                                       keys = (featureNames(data_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_data <- subset(anno_data, !is.na(SYMBOL))


# Removing multiple mappings
anno_grouped <- group_by(anno_data, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

print(head(anno_summarized))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)

head(anno_filtered)

probe_stats <- anno_filtered 

print(nrow(probe_stats))

ids_to_exlude <- (featureNames(data_manfiltered) %in% probe_stats$PROBEID)

print(table(ids_to_exlude))

# an expression set without the ids_to_exclude.
data_final <- subset(data_manfiltered, !ids_to_exlude)

print(validObject(data_final))

# exclude them from the feature data 
print(head(anno_data))
fData(data_final)$PROBEID <- rownames(fData(data_final))
# left join
fData(data_final) <- left_join(fData(data_final), anno_data)
rownames(fData(data_final)) <- fData(data_final)$PROBEID 

print("Final data is valid:")
print(validObject(data_final))
print("Final data dimensions:")
print(dim(exprs(data_final)))

# assign probe ids with gene symbols
norm_mat <- as.data.frame(exprs(data_final))
fdata_mat <- as.data.frame(fData(data_final)[,-3])

fdata_norm_mat <- merge(fdata_mat, norm_mat, by= "row.names")[,-c(1,2)]

fdata_norm_mat_gene <- as.data.frame(limma::avereps(fdata_norm_mat[,2:ncol(fdata_norm_mat)], ID = fdata_norm_mat$SYMBOL))
sg <- as.data.frame(pData(data_final))
  
print(dim(fdata_norm_mat_gene))
print(dim(sg))

# limma package for differential expression analysis
sg$drug.response <- factor(sg$drug.response)
design <- model.matrix(~0 + sg$drug.response)
colnames(design) <- levels(sg$drug.response)

fit <- lmFit(fdata_norm_mat_gene, design)  # fit each probeset to model

contrast.matrix<-makeContrasts(R-S,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
output <- topTable(fit2, number=Inf, sort.by = "P", adjust.method="BH")

print(length(which(output$adj.P.Val < 0.05)))

# Summary of results (number of differentially expressed genes)
result <- decideTests(fit2)
print(summary(result))


#### read annotation file for chip
#### Note download file "HG-U219.na36.annot" into the data folder 
#### http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus

anno_affy <- read.csv("../data/HG-U219.na36.annot.csv", header=T, comment.char = '#')
anno_affy_sub <- anno_affy[c("Probe.Set.ID","Transcript.Assignments")]
anno_affy_sub$type <- ifelse(grepl('protein_coding', anno_affy_sub$Transcript.Assignments),'protein_coding','non_coding' )

anno_affy_sub <- merge(anno_affy_sub,fdata_mat,by.x="Probe.Set.ID", by.y= "PROBEID")
anno_affy_sub <- unique(anno_affy_sub [c("type","SYMBOL")])


###########################################################
##########################################################
### DE for 100 cuts of the data

# format the data to have id, group, gene1,..,genen as rows and cell lines as columns
input.class <- function(df.sub,infile.sub.m2)
{
  df.sub.t <- as.data.frame(t(df.sub))
  df.sub.t.grp <- merge(infile.sub.m2,df.sub.t,by.x="id", by.y="row.names")
  print(dim(df.sub.t.grp))
  names(df.sub.t.grp)[1:2] <- c("id","group")
  df.sub.t.grp$group <- as.factor(df.sub.t.grp$group)
  return(df.sub.t.grp)
}

sg.sub <- sg["drug.response"]
sg.sub$id <- row.names(sg.sub)
df.sub.t.grp <- input.class(fdata_norm_mat_gene,sg.sub)

# split data into train and test
set.seed(123)
inTrain <- createDataPartition(df.sub.t.grp$group, p = .8, list = FALSE, times=100)

for( cut in c(1:100)){
  
  # train data
  train <- df.sub.t.grp[ inTrain[,cut], ]
  print(dim(train))
  trainClass <- df.sub.t.grp$group[ inTrain[,cut]]
  # test data
  test  <- df.sub.t.grp[-inTrain[,cut], ]
  print(dim(test))
  testClass  <- df.sub.t.grp$group[-inTrain[,cut]]
  
  # transpose
  ctmg0 <- as.data.frame(t(train[3:ncol(train)]))
  names(ctmg0) <- train$id
  sg <- as.data.frame(train[c(1:2)])
  
  print(dim(ctmg0))
  print(dim(sg))
  
  #DE on train
  design <- model.matrix(~0 + sg$group)
  colnames(design) <- levels(sg$group)
  
  fit <- lmFit(ctmg0, design)  # fit each probeset to model
  
  
  contrast.matrix<-makeContrasts(R-S,levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  output <- topTable(fit2, number=Inf, sort.by = "P", adjust.method="BH")
  output <- merge(output,anno_affy_sub, by.x="row.names", by.y="SYMBOL")
  
  print(length(which(output$adj.P.Val < 0.05)))
 
  write.csv(output, file = paste0(out, "/",prefix,"_affy_limma_train_",cut,".csv"),quote=F,row.names=T)
  
  # sig protein coding genes with FDR < 0.05
  tlist <- as.data.frame(output[output$adj.P.Val<0.05 & output$type == "protein_coding" ,c("Row.names")])
  colnames(tlist) <- "gene"
  print(dim(tlist))
  tlist_list[[cut]] <- tlist$gene
  
  if(nrow(tlist) != 0){
    # get subset for train data sig genes
    nct1.train.sig.t.grp <- train[,colnames(train) %in% c("id", "group", tlist$gene)]
    print(dim(nct1.train.sig.t.grp))
    
    # get subset for test data sig genes
    nct1.test.sig.t.grp <- test[,colnames(test) %in% c("id", "group", tlist$gene)]
    print(dim(nct1.test.sig.t.grp))
    
    
    # save
    write.csv(nct1.train.sig.t.grp,paste0(out, prefix, "_Affy_DE_TrainingSet_new_",cut,".csv"),row.names = T,quote = F)
    write.csv(nct1.test.sig.t.grp,paste0(out,prefix,"_Affy_DE_TestingSet_new_",cut,".csv"),row.names = T, quote = F)
    
    
    
    ############################################################
    ## feature selection using boruta
    x1<- as.data.frame.matrix(nct1.train.sig.t.grp[3:ncol(nct1.train.sig.t.grp)])
    row.names(x1) <- nct1.train.sig.t.grp$id
    normalization <- preProcess(x1)
    x1 <- predict(normalization, x1)
    x1 <- as.data.frame(x1)
    y1 <- sg$group
    y1 <- mapvalues(sg$group , from = c("R","S"  ), to = c(-1,1))
    names(y1) <- sg$id
    
    
    # remove features with 0 variance or v.few uniq values
    zero.var <- nearZeroVar(x1)
    if(length(zero.var) != 0){x1 <- x1[,-zero.var]}
    # remove attributes with an absolute correlation of 0.8 or higher.
    #Remove Redundant Features
    if(nrow(tlist) >=2 ){
      cor.fe <- findCorrelation(cor(x1, method="spearman"), .9)
      if(length(cor.fe) != 0){x1 <- x1[, -cor.fe]}
      print(dim(x1))
    }
    
    ### Boruta feature ranking and selection algorithm based on random forests algorithm
    boruta_output <- Boruta(x1,y1, doTrace=0)  
    names(boruta_output)
    boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
    print(boruta_signif)  
    
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_signif <- getSelectedAttributes(roughFixMod)
    print(boruta_signif)
    print(length(boruta_signif))
    
    # Variable Importance Scores
    imps <- attStats(roughFixMod)
    imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
    print(head(imps2[order(-imps2$meanImp), ]))  # descending sort
    
    boruta_signif_list[[cut]] <- boruta_signif
    
    # get subset train norm matrix for sig genes
    # get subset test norm matrix for sig genes
    # get subset for train data sig genes
    boruta.train.sig.t.grp <- train[,colnames(train) %in% c("id", "group", boruta_signif)]
    print(dim(boruta.train.sig.t.grp))
    print(table(boruta.train.sig.t.grp$group))
    
    # get subset for test data sig genes
    boruta.test.sig.t.grp <- test[,colnames(test) %in% c("id", "group", boruta_signif)]
    print(dim(boruta.test.sig.t.grp))
    print(table(boruta.test.sig.t.grp$group))
    
    # save
    write.csv(boruta.train.sig.t.grp,paste0(out, prefix, "_Affy_boruta_TrainingSet_new_",cut,".csv"),row.names = T,quote = F)
    write.csv(boruta.test.sig.t.grp,paste0(out,prefix,"_Affy_boruta_TestingSet_new_",cut,".csv"),row.names = T, quote = F)
  }
  
}

saveRDS(boruta_signif_list, file=paste0(out, prefix, "_Affy_boruta_signif_list100.RDS"))
saveRDS(tlist_list, file=paste0(out, prefix, "_Affy_DE_signif_list100.RDS"))

common.ids=tlist_list[[1]]
for(i in 2:100)
{
  common.ids <- Reduce(intersect,list(tlist_list[[i]] , common.ids))
  print(length(common.ids))
}

print(length(common.ids))
print((common.ids))
saveRDS(common.ids, file=paste0(out, prefix, "_Affy_DE_signif_common_ids.RDS"))


#### tlist list length
list.length <- c()
for(index in 1:100)
{
  list.length <- append(list.length, length(tlist_list[[index]]))
}

print(summary(list.length))
paste("Time to run the whole workflow:",intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))
