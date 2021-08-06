# Sharvari Gujja
# script to generate normalized input matrix for R vs. NR classification for a specific drug and
# generate ML models for 100 cuts of train and test
# data columns formatted as ID, group, feature1,.., featureN
# Usage: Rscript <script> <drugname> <outputDir>
# Example: Rscript get.R.NR.input.for.classifier.100cuts.R Gemcitabine TCGA_Gemcitabine_100cuts/
# Output: Summary table of classification metrics across 100 runs

library(dplyr)
library(edgeR)
library(ggplot2)
library(umap)
library(matrixStats)
library(plyr)
library(Boruta)
library(caret)
library(sigFeature)

library(humanid)
library(data.table)
library(tidyr)

rm(list = ls())
gc()

args <- commandArgs(TRUE)

# select drug name
prefix <- args[1] 
# output dir
out <- args[2] 
if(!dir.exists(out)) { dir.create(out) }

logFC_cutoff <- as.numeric(args[3])

## list for boruta_signif
boruta_signif_list = list()
## list for DE sig genes
tlist_list = list()

#### format clinical group info
########################
# read R and NR data for a drug

infile <- read.csv("../data/raw/drug_response_ding_etal_2016.csv", header=T, sep=",")
print(dim(infile))

# remove duplicated rows
infile <- infile %>% distinct()

# create new variable group for NR and R
infile<- infile %>% mutate(group = case_when(response == "Clinical Progressive Disease"|response == "Stable Disease"  ~ "non.responders", 
                                             response == "Complete Response" | response == "Partial Response" ~ "responders", TRUE ~ "NA"))

print(dim(infile))

# subset for drug name
infile.sub <- infile[infile$drug.name %in% prefix,]
print(dim(infile.sub))
print(table(infile.sub$cancers, infile.sub$group))

# remove cancers with < 10 samples in each group
cancers.sel <- as.data.frame.matrix(table(infile.sub$cancers, infile.sub$group)) 
cancers.sel <- cancers.sel[(cancers.sel$non.responders >= 10 & cancers.sel$responders >= 10) ,]
infile.sub <- infile.sub[infile.sub$cancers %in% row.names(cancers.sel), ]
print(table(infile.sub$cancers, infile.sub$group))

# select one cancer type:
#infile.sub <- infile.sub[infile.sub$cancers %in% "PAAD",]
infile.sub <- infile.sub[c("cancers", "patient.arr" ,"group","drug.name")] %>% distinct()
print(dim(infile.sub))
print(table(infile.sub$group))

#####
# > head(infile.sub)
# cancers  patient.arr          group drug.name
# 1    BLCA TCGA-G2-A2EC     responders Cisplatin
# 2    BLCA TCGA-G2-A2EJ non.responders Cisplatin
# 3    BLCA TCGA-G2-A2EF     responders Cisplatin

#####

#### read TCGA reference annotation file for gene symbol and gene type
gtf.file <- "../data/annot//gencode.v22.annotation.gtf.gz"
gtf.df <- fread(gtf.file, header=F,sep="\t", skip=1)
gtf.df.gene <- gtf.df[gtf.df$V3 == "gene", "V9"]

find.string <- paste(unlist(list("gene_id ", " gene_type ", " gene_status ", " gene_name ", "\\\"")), collapse = "|")
gtf.df.gene$V9 <- gsub(find.string, replacement = "", gtf.df.gene$V9)

gtf.df.gene.split <- separate(gtf.df.gene, col = V9, into = c("gene_id", "gene_type", "gene_status","gene_name"), sep = "\\;")
gtf.df.gene.split <- gtf.df.gene.split  %>% 
  mutate(across(where(is.character), str_trim))

# > head(gtf.df.gene.split)
# gene_id                          gene_type gene_status    gene_name
# 1: ENSG00000223972.5 transcribed_unprocessed_pseudogene       KNOWN      DDX11L1
# 2: ENSG00000227232.5             unprocessed_pseudogene       KNOWN       WASH7P
# 3: ENSG00000278267.1                              miRNA       KNOWN    MIR6859-3

###### TCGA raw count and sample group info
df <- readRDS("../data/raw/tcga_rna_raw_counts_33cancers.rds")
print(dim(df))
#[1] 60483 11093

# sample group
files <- list.files(pattern = "htseq_manifest.csv", recursive = TRUE, path="../data/manifest_files//", full.names = TRUE)
tables <- lapply(files, read.csv, header = TRUE)
mf <- do.call(rbind , tables)

# replace df colnames with patient barcodes
mf$file_name <- gsub(".htseq.counts.gz","",mf$file_name )
# subset for Primary Tumors
mf.sub <- mf[mf$sample_type == "Primary Tumor" ,]
print(dim(mf.sub))

# merge
infile.sub.m <- merge(infile.sub, mf.sub, by.x="patient.arr", by.y="cases.submitter_id")
print(dim(infile.sub.m))
print(table(infile.sub.m$drug.name, infile.sub.m$group))

# subset count matrix for the patients on the candidate drug
df.sub <- df[, colnames(df) %in% infile.sub.m$file_name]
print(dim(df.sub))

# sort sample id in sg in the same order as row names in count matrix
infile.sub.m <-infile.sub.m[order(match(infile.sub.m$file_name,colnames(df.sub))),]
print(dim(infile.sub.m))
print(table(infile.sub.m$project,infile.sub.m$group))

# > str(infile.sub.m)
# 'data.frame':	135 obs. of  32 variables:
#   $ patient.arr              : chr  "TCGA-G2-A2EC" "TCGA-FJ-A3Z9" "TCGA-FD-A3SO" "TCGA-DK-A3WW" ...
# $ cancers                  : chr  "BLCA" "BLCA" "BLCA" "BLCA" ...
# $ group                    : chr  "responders" "responders" "non.responders" "responders" ...
# $ drug.name                : chr  "Cisplatin" "Cisplatin" "Cisplatin" "Cisplatin" ...
# $ id                       : chr  "01eff87c-9e9b-4cc8-964d-9c0aa137424d" "026112a8-334d-4c82-90b0-006d48499069" "0a3135cb-f854-4944-898e-1985754c754f" "1274623c-3d90-4a0e-8f62-967c08a8c736" ...


# create DGE object for all samples
y <- DGEList(counts=df.sub, group=infile.sub.m$group)
# design matrix
design <- model.matrix(~0+group, data=infile.sub.m)
rownames(design) <- colnames(y)

#e filtering low expressed genes
keep <- filterByExpr(y, design)
print(table(keep))
df.sub <- df.sub[keep,]
print(dim(df.sub))

infile.sub.m <- infile.sub.m[,c("file_name","group")]
names(infile.sub.m) <- c("id","group")

# generate input for classification
# transponse to get rows as samples and columns as genes to merge with sample group
input.class <- function(df.sub,infile.sub.m2)
{
  df.sub.t <- as.data.frame(t(df.sub))
  #infile.sub.m2 <- infile.sub.m[,c("file_name","group")]
  df.sub.t.grp <- merge(infile.sub.m2,df.sub.t,by.x="id",by.y="row.names")
  print(dim(df.sub.t.grp))
  #names(df.sub.t.grp)[1:2] <- c("id","group")
  df.sub.t.grp$group <- as.factor(df.sub.t.grp$group)
  return(df.sub.t.grp)
}

df.sub.t.grp <- input.class(df.sub,infile.sub.m)

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
  
  ctmg0 <- as.data.frame(t(train[3:ncol(train)]))
  names(ctmg0) <- train$id
  sg <- as.data.frame(train[c(1:2)])
  
  print(dim(ctmg0))
  print(dim(sg))
  
  # TMM normalization of gene counts
  
  # create DGE object for all samples
  y <- DGEList(counts=ctmg0, group=sg$group)
  
  # design matrix
  design <- model.matrix(~0+group, data=sg)
  print(colnames(design))
  rownames(design) <- colnames(y)
  
  # normalization
  y <- calcNormFactors(y,method='TMM')
  nct1 <- cpm(y,normmalized.lib.sizes=TRUE,log=T)
  nct1 <- data.frame(round(nct1,3), check.names=F)
  print(dim(nct1))
  write.table(nct1, file = paste0(out,"/", prefix,"_res_nores_TMM_count_matrix_train_",cut,".txt"),quote=F,row.names=T, sep="\t")
  
  # create GSEA input
  #createGSEAinput(prefix,nct1,sg$group)
  
  # estimate dispersion
  y <- estimateDisp(y,design,robust=TRUE)
  
  # contrasts
  my.contrasts <- makeContrasts(
    NR.vs.R = groupnon.responders-groupresponders,
    levels=design
  )
  
  #DE
  fit.ql <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit.ql, contrast=my.contrasts)
  t.qlf <- as.data.frame(topTags(qlf,n=Inf))
  
  t.qlf <- merge(t.qlf, gtf.df.gene.split, by.x="row.names", by.y="gene_id")
  t.qlf <- t.qlf[order(t.qlf$FDR),]
  row.names(t.qlf) <- t.qlf$Row.names
  write.csv(t.qlf, file = paste0(out, "/",prefix,"_edgerOut_QL_train_",cut,".csv"),quote=F,row.names=T)
  
  # sig protein coding genes with FDR < 0.05
  tlist <- as.data.frame(row.names(t.qlf[t.qlf$FDR<0.05 & t.qlf$gene_type == "protein_coding" & abs(t.qlf$logFC) >= logFC_cutoff ,]))
  colnames(tlist) <- "gene"
  print(dim(tlist))
  tlist_list[[cut]] <- tlist$gene
  
  # get subset TMM matrix for sig genes
  nct1.sig <- nct1[row.names(nct1) %in% tlist$gene,]
  print(dim(nct1.sig))
  
  nct1.sig.t.grp <- input.class(nct1.sig,sg)
  
  ########
  ### normalize test
  ctmg1 <- as.data.frame(t(test[3:ncol(test)]))
  names(ctmg1) <- test$id
  sg1 <- as.data.frame(test[c(1:2)])
  
  print(dim(ctmg1))
  print(dim(sg1))
  
  # TMM normalization of gene counts
  
  # create DGE object for all samples
  y <- DGEList(counts=ctmg1, group=sg1$group)
  
  # design matrix
  design <- model.matrix(~0+group, data=sg1)
  print(colnames(design))
  rownames(design) <- colnames(y)
  
  # normalization
  y <- calcNormFactors(y,method='TMM')
  nct1.test <- cpm(y,normmalized.lib.sizes=TRUE,log=T)
  nct1.test <- data.frame(round(nct1.test,3), check.names=F)
  print(dim(nct1.test))
  
  # get subset TMM matrix for sig genes
  nct1.test.sig <- nct1.test[row.names(nct1.test) %in% tlist$gene,]
  print(dim(nct1.test.sig))
  
  nct1.sig.test.t.grp <- input.class(nct1.test.sig,sg1)
  
  # save
  write.csv(nct1.sig.t.grp,paste0(out, "/",prefix, "_DE_TrainingSet_new_",cut,".csv"),row.names = T,quote = F)
  write.csv(nct1.sig.test.t.grp,paste0(out,"/",prefix,"_DE_TestingSet_new_",cut,".csv"),row.names = T, quote = F)
  
 
  
  ############################################################
  ## feature selection using boruta
  x1<- as.matrix(t(nct1.sig))
  normalization <- preProcess(x1)
  x1 <- predict(normalization, x1)
  x1 <- as.data.frame(x1)
  y1 <- sg$group
  y1 <- mapvalues(sg$group , from = c("non.responders","responders"  ), to = c(-1,1))
  names(y1) <- sg$id
  
  
  # remove features with 0 variance or v.few uniq values
  zero.var <- nearZeroVar(x1)
  if(length(zero.var) != 0){x1 <- x1[,-zero.var]}
  # remove attributes with an absolute correlation of 0.8 or higher.
  #Remove Redundant Features
  cor.fe <- findCorrelation(cor(x1, method="spearman"), .9)
  if(length(cor.fe) != 0){x1 <- x1[, -cor.fe]}
  print(dim(x1))
  
  
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
  
  # get subset train TMM matrix for sig genes
  boruta.sig <- nct1[row.names(nct1) %in% boruta_signif,]
  print(dim(boruta.sig))
  
  boruta.sig.t.grp <- input.class(boruta.sig,sg)
  
  # get subset test TMM matrix for sig genes
  boruta.test.sig <- nct1.test[row.names(nct1.test) %in% boruta_signif,]
  print(dim(boruta.test.sig))
  
  boruta.sig.test.t.grp <- input.class(boruta.test.sig,sg1)
  # save
  write.csv(boruta.sig.t.grp,paste0(out, "/",prefix, "_boruta_TrainingSet_new_",cut,".csv"),row.names = T,quote = F)
  write.csv(boruta.sig.test.t.grp,paste0(out,"/",prefix,"_boruta_TestingSet_new_",cut,".csv"),row.names = T, quote = F)
  

}

saveRDS(boruta_signif_list, file=paste0(out,"/", prefix, "_boruta_signif_list100.RDS"))
saveRDS(tlist_list, file=paste0(out, "/",prefix, "_DE_signif_list100.RDS"))

common.ids=tlist_list[[1]]
for(i in 2:100)
{
  common.ids <- Reduce(intersect,list(tlist_list[[i]] , common.ids))
  print(length(common.ids))
}

print(length(common.ids))


##### Boruta feature selection
#length(Reduce(intersect,list(tlist_list[[1]] , tlist_list[[2]], tlist_list[[3]],tlist_list[[4]], tlist_list[[5]], tlist_list[[6]],tlist_list[[7]],tlist_list[[8]],tlist_list[[9]],tlist_list[[10]],tlist_list[[11]],tlist_list[[12]],tlist_list[[13]],tlist_list[[14]],tlist_list[[15]],tlist_list[[16]],tlist_list[[17]],tlist_list[[18]])))
x1 <- as.data.frame(t(nct1.sig[row.names(nct1.sig) %in% common.ids,]))
normalization <- preProcess(x1)
x1 <- predict(normalization, x1)
x1 <- as.data.frame(x1)
y1 <- sg$group
y1 <- mapvalues(sg$group , from = c("non.responders","responders"  ), to = c(-1,1))
names(y1) <- sg$id


# remove features with 0 variance or v.few uniq values
zero.var <- nearZeroVar(x1)
if(length(zero.var) != 0){x1 <- x1[,-zero.var]}
# remove attributes with an absolute correlation of 0.8 or higher.
#Remove Redundant Features
cor.fe <- findCorrelation(cor(x1, method="spearman"), .9)
if(length(cor.fe) != 0){x1 <- x1[, -cor.fe]}
print(dim(x1))


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

boruta_signif.df <- subset(gtf.df.gene.split ,gtf.df.gene.split$gene_id %in% row.names(imps2) )
print(boruta_signif.df)

write.csv(boruta_signif.df,paste0(out,"/",prefix,"_DEcommon_100cuts_boruta_genes.csv"),row.names = F, quote = F)



#### tlist list length
list.length <- c()
for(index in 1:100)
{
  list.length <- append(list.length, length(tlist_list[[index]]))
}

print(summary(list.length))
