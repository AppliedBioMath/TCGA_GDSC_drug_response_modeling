# Sharvari Gujja
# script to generate normalized input matrix for R vs. NR classification for a specific drug and
# generate ML models for 100 cuts of train and test
# data columns formatted as ID, group, feature1,.., featureN
# Usage: Rscript <script> <drugname> <outputDir>
# Example: Rscript get.R.NR.input.for.classifier.100cuts.R Gemcitabine TCGA_Gemcitabine_100cuts/
# Output: Summary table of classification metrics across 100 runs
# FE of DE genes that are sig in atleast 50% of the training cuts...

library(dplyr)
library(edgeR)
library(ggplot2)
library(matrixStats)

library(caret)
library(sigFeature)

library(humanid)
library(data.table)
library(tidyr)

library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)
library(PRROC)
library(plyr)


library(parallel)
library(doParallel)

library(org.Hs.eg.db)
library(clusterProfiler)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

Sys.setenv(TZ="US/Eastern")
intervalStart <- Sys.time()

rm(list = ls())
gc()

args <- commandArgs(TRUE)


##################################################################
# input arguments
# select drug name
prefix <- "Gemcitabine" #args[1] 
# output dir
out <- "test_output" #args[2] 

# create output directory if it doesn't exist
if(!dir.exists(out)) { dir.create(out) }

# log FC cut-off to use for filtering sig DE genes
#logFC_cutoff <- as.numeric(2)

## list for DE sig genes for each cut of train 
tlist_list = list()

##########################################################################
# Arguments for ML algorithms
pos_class = "responders"
### input ML algorithms
methods = c('glmnet','glmnet','glmnet','svmLinear2','svmRadial','rf','nb')
method_name = c('lasso','ridge','enet','svmLinear','svmRadial','rf','nb')

lassoGrid=expand.grid(.alpha=1, .lambda=10^seq(10, -2, length = 100)) 
ridgeGrid=expand.grid(.alpha=0,.lambda=10^seq(10, -2, length = 100)) 
enetGrid=expand.grid(.alpha=0.5,.lambda=10^seq(10, -2, length = 100))
svmRGrid = expand.grid(sigma = c(.01, .015, 0.2), C = c(0.75, 0.9, 1, 1.1, 1.25))

cvCtrl = trainControl(method = "cv", number = 10,classProbs = TRUE,allowParallel = TRUE)

# Stats Mat
summary_df=data.frame(method=character(),train_acc=double(),test_acc=double(),
                        train_bacc=double(),test_bacc=double(),train_auroc=double(), test_auroc=double(),
                        train_prec=double(), test_prec=double(), train_recall=double(), test_recall=double(),
                        train_F1=double(),test_F1=double(),train_Sensitivity=double(),test_Sensitivity=double(),
                        train_Specificity=double(),test_Specificity=double(),
                        train_PPV=double(),test_PPV=double(),
                        train_NPV=double(),test_NPV=double(),stringsAsFactors=FALSE)

n=1
###############################################################################3

# DATA PROCESSING
########################
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
print(paste0("Responders and Non-responders for the drug ",prefix," by cancer type:"))
print(table(infile.sub$cancers, infile.sub$group))

# select one cancer type:
#infile.sub <- infile.sub[infile.sub$cancers %in% "PAAD",]
infile.sub <- infile.sub[c("cancers", "patient.arr" ,"group","drug.name")] %>% distinct()
print(dim(infile.sub))
print(paste0("Total Responders and Non-responders for the drug ",prefix,":"))
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
print("Genes and samples for raw expression counts from TCGA:")
print(dim(df))

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
print(paste0("Responders and Non-responders for the drug ",prefix," with raw counts:"))
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
print("Filtering low expressed genes:")
print(table(keep))
df.sub <- df.sub[keep,]
print(dim(df.sub))

infile.sub.m <- infile.sub.m[,c("file_name","group")]
names(infile.sub.m) <- c("id","group")

# generate input for classification
# transpose to get rows as samples and columns as genes to merge with sample group
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

# set seed for reproducibility
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

  # normalization
  y <- calcNormFactors(y,method='TMM')
  nct1 <- cpm(y,normmalized.lib.sizes=TRUE,log=T)
  nct1 <- data.frame(round(nct1,3), check.names=F)
  print(dim(nct1))
  print(paste0("Saving normalized count matrix for training set ", cut, ":"))
  write.table(nct1, file = paste0(out,"/", prefix,"_res_nores_TMM_count_matrix_train_",cut,".txt"),quote=F,row.names=T, sep="\t")
  
  # remove features with 0 variance or v.few uniq values
  print("Remove features with near zero variance:")
  # transpose normalized count matrix to have columns as genes
  nct1.t <- as.data.frame(t(nct1))
  zero.var <- nearZeroVar(nct1.t)
  if(length(zero.var) != 0){nct1.t <- nct1.t[,-zero.var]}
  # remove attributes with an absolute correlation of 0.9 or higher.
  #Remove Redundant Features
  cor.fe <- findCorrelation(cor(nct1.t, method="spearman"), .9)
  if(length(cor.fe) != 0){nct1.t <- nct1.t[, -cor.fe]}
  print("Dimensions of normalized count matrix after removing features with near zero variance and correlation > 0.9:")
  print(dim(nct1.t))
  
  # subset raw count matrix with the same features retained after filtering
  ctmg0.sub <- ctmg0[row.names(ctmg0) %in% colnames(nct1.t),]
  
  print(paste0("Run DGE analysis with Responders as reference for training set ", cut, ":"))
  # create DGE object for all samples
  y <- DGEList(counts=ctmg0.sub, group=sg$group)
  
  # design matrix
  design <- model.matrix(~0+group, data=sg)
  print(colnames(design))
  rownames(design) <- colnames(y)
  
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
  print(paste0("Saving DGE output for training set ", cut, ":"))
  write.csv(t.qlf, file = paste0(out, "/",prefix,"_edgerOut_QL_train_",cut,".csv"),quote=F,row.names=T)
  
  print("Filtering for sig DE genes with FDR < 0.05") #that are protein coding with FDR < 0.05 and abs(logFC) >= 2")
  # sig protein coding genes with FDR < 0.05 and logFC >= |2|
  #tlist <- as.data.frame(row.names(t.qlf[t.qlf$FDR < 0.05 & t.qlf$gene_type == "protein_coding" & abs(t.qlf$logFC) >= logFC_cutoff ,]))
  tlist <- as.data.frame(row.names(t.qlf[t.qlf$FDR < 0.05 & t.qlf$gene_type == "protein_coding" ,]))
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
  print(paste0("Saving train and test data set for cut ", cut))
  write.csv(nct1.sig.t.grp,paste0(out, "/",prefix, "_DE_TrainingSet_new_",cut,".csv"),row.names = T,quote = F)
  write.csv(nct1.sig.test.t.grp,paste0(out,"/",prefix,"_DE_TestingSet_new_",cut,".csv"),row.names = T, quote = F)
  
  
  ### Code to run multiple machine learning algorithms
  train_data = nct1.sig.t.grp
  test_data = nct1.sig.test.t.grp
  
  y_train = train_data[[2]]
  y_train = as.factor(y_train)
  x_train = train_data[,-c(1,2)]
  
  y_test = test_data[[2]]
  y_test = as.factor(y_test)
  x_test = test_data[,-c(1,2)]

  
  # Scale Data : Apply mean and SD of train to test
  cat("\n---------------\nScale Data...\n---------------\n")
  train_mean = apply(x_train, 2, mean)
  train_sd = apply(x_train, 2, sd)
  # Train
  z_x_train = sweep(x_train, 2, train_mean, "-")
  z_x_train = sweep(z_x_train, 2, train_sd, "/")
  # Test
  z_x_test = sweep(x_test, 2, train_mean, "-")
  z_x_test = sweep(z_x_test, 2, train_sd, "/")
  
  
  # make group/label as categorical
  y_train = as.factor(y_train)
  y_test = factor(y_test, levels = levels(y_train))
  
  
  # RF
  rfGrid=expand.grid(.mtry=sqrt(ncol(x_train)))
  tuneGrids=list(lassoGrid,ridgeGrid,enetGrid,NULL,svmRGrid,rfGrid,NULL)
  
  print("Running ML models:")
  for (m in 1:length(methods)) {
    method=methods[m]
    print(method_name[[m]])
    
    ######### Model Runs
    fit = train(z_x_train,y_train,method=method,trControl = cvCtrl,tuneGrid = tuneGrids[[m]])
    
    # Class returns the class of the prediction 
    pred_train = predict(fit,z_x_train)
    pred_test = predict(fit,z_x_test)
    
    # Class Stats
    confmat_train = confusionMatrix(pred_train, y_train, positive = as.character(pos_class))
    confmat_test = confusionMatrix(pred_test, y_test, positive = as.character(pos_class))
    
    # calculate response probabilities
    response_train = predict(fit,z_x_train,type="prob")
    response_test = predict(fit,z_x_test,type="prob")
    
    # AUC train and test
    auc.train = ModelMetrics::auc(as.numeric(y_train), predict(fit, newdata = z_x_train, type = 'prob')[,1])
    auc.test = ModelMetrics::auc(as.numeric(y_test), predict(fit, newdata = z_x_test, type = 'prob')[,1])
    
    summary_df[n,'method']=method_name[[m]]
    summary_df[n,'train_acc']=confmat_train$overall["Accuracy"]
    summary_df[n,'test_acc']=confmat_test$overall["Accuracy"]
    summary_df[n,'train_bacc']=confmat_train$byClass['Balanced Accuracy']
    summary_df[n,'test_bacc']=confmat_test$byClass['Balanced Accuracy']
    summary_df[n,'train_prec']=confmat_train$byClass['Precision']
    summary_df[n,'test_prec']=confmat_test$byClass['Precision']
    summary_df[n,'train_recall']=confmat_train$byClass['Recall']
    summary_df[n,'test_recall']=confmat_test$byClass['Recall']
    summary_df[n,'train_F1']=confmat_train$byClass['F1']
    summary_df[n,'test_F1']=confmat_test$byClass['F1']
    summary_df[n,'train_Sensitivity']=confmat_train$byClass['Sensitivity']
    summary_df[n,'test_Sensitivity']=confmat_test$byClass['Sensitivity']
    summary_df[n,'train_Specificity']=confmat_train$byClass['Specificity']
    summary_df[n,'test_Specificity']=confmat_test$byClass['Specificity']
    summary_df[n,'train_PPV']=confmat_train$byClass['Pos Pred Value']
    summary_df[n,'test_PPV']=confmat_test$byClass['Pos Pred Value']
    summary_df[n,'train_NPV']=confmat_train$byClass['Neg Pred Value']
    summary_df[n,'test_NPV']=confmat_test$byClass['Neg Pred Value']
    summary_df[n,'train_auroc']=auc.train
    summary_df[n,'test_auroc']=auc.test
    
    n = n+1 
  }

}


# aggregate mean and sem for each method and metrics
DF <- as.data.frame(summary_df)
means <- aggregate( .~ method, DF, mean)
rownames(means) <- paste0("mean_", means$method)
sds <- aggregate( .~ method, DF, sd)
rownames(sds) <- paste0("sd_", means$method)
tab <- rbind(means, sds)
tab <- tab[, -1]
tab <- t(tab)
print(tab)

# save output files
saveRDS(tlist_list, file=paste0(out, "/",prefix, "_DE_signif_list100.RDS"))
print(dim(summary_df))
write.csv(summary_df,paste0(out,"/",prefix,"_DE_100runs_ML_metrics_info.csv"),quote=F,row.names=F )
write.csv(as.data.frame(tab),paste0(out,"/",prefix,"_DE_100runs_ML_metrics_info_mean_SD.csv"),quote=F,row.names=F )


#### sig DE genes range: tlist list length
print("Range of sig. DE genes:")
list.length <- c()
for(index in 1:100)
{
  list.length <- append(list.length, length(tlist_list[[index]]))
}

print(summary(list.length))

######### FUNCTIONAL ENRICHMENT ANALYSIS

## common sig genes across all 100 runs
common.ids=tlist_list[[1]]
for(i in 2:100)
{
  common.ids <- Reduce(intersect,list(tlist_list[[i]] , common.ids))
}

print("Common sig DE genes across all 100 runs:")
print(length(common.ids))

gtf.df.gene.split.common <- gtf.df.gene.split[gtf.df.gene.split$gene_id %in% common.ids,]
print(gtf.df.gene.split.common)

## union of sig genes across all 100 runs
union.ids=tlist_list[[1]]
for(i in 2:100)
{
  union.ids <- Reduce(union,list(tlist_list[[i]] , union.ids))
}

union.ids <- unique(union.ids)
print("Union sig DE genes across all 100 runs:")
print(length(union.ids))

gtf.df.gene.split.union <- gtf.df.gene.split[gtf.df.gene.split$gene_id %in% union.ids,]
print(gtf.df.gene.split.union)


### Run GO enrichment analysis
# create tlist and blist

input.list <- gsub("\\..*", "" ,gtf.df.gene.split.union$gene_id)
print(length(input.list))
bg.list <- gsub("\\..*", "" ,row.names(df.sub))
print(length(bg.list))

## Run GO enrichment analysis
ego <- enrichGO(gene = input.list,
                universe = bg.list,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.01,
                minGSSize = 2,
                maxGSSize = 500,
                readable = TRUE,
                pool = TRUE)

## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
print((dim(cluster_summary)))

write.csv(cluster_summary, "brca_matched_pairs_clusterProfiler_go_2.csv")


##### Calculating time to run the whole workflow

intervalEnd <- Sys.time()
paste("Time to run the whole pipeline for 100 iterations of Train and Test:",intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))

