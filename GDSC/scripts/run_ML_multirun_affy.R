# Sharvari Gujja
# Script to run ML algorithms for 100 iterations
# Usage: Rscript <script> <drugname> <outputDir>
# Example: Rscript run_ML_multirun_affy.R Gemcitabine GDSC_Gemcitabine_100cuts/
# Output: Summary table of classification metrics across 100 runs

rm(list = ls())
gc()
library(feather)
library(caret)
library(stringr)
library(HandTill2001)
library(MLmetrics)
library(PRROC)
library(plyr)


library(parallel)
library(doParallel)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

Sys.setenv(TZ="US/Eastern")
intervalStart <- Sys.time()

args <- commandArgs(TRUE)

prefix = args[1] #"Gemcitabine" # 
block_path = args[2] #"test_output" #


pos_class = "S"

### stuff for cv
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


n=1;

# Load Data -- NR vs. R
for(cut in c(1:100))
{
  if (file.exists(paste0(block_path,"/",prefix, '_Affy_DE_TrainingSet_new_',cut,'.csv'))){
    train_data = read.csv(paste0(block_path,"/",prefix, '_Affy_DE_TrainingSet_new_',cut,'.csv'), row.names=1, header=T)
    test_data = read.csv(paste0(block_path,"/",prefix, '_Affy_DE_TestingSet_new_',cut,'.csv'), row.names=1, header=T)
    
    if (ncol(train_data) >= 4) {
      y_train = train_data[[2]]
      y_train = as.factor(y_train)
      x_train = train_data[,-c(1,2)]
      
      y_test = test_data[[2]]
      y_test = as.factor(y_test)
      x_test = test_data[,-c(1,2)]
      
      # RF
      rfGrid=expand.grid(.mtry=sqrt(ncol(x_train)))
      tuneGrids=list(lassoGrid,ridgeGrid,enetGrid,NULL,svmRGrid,rfGrid,NULL)
      
      ### Main Model Loop
      cat(paste('\n', format(Sys.time(), "%H:%M"), '\n', sep=""))
      
      
      # Remove Zero Variance
      train_var = apply(x_train, 2, sd)
      no_var_idx = which(train_var == 0.0)
      if (length(no_var_idx) > 0) {
        x_train = x_train[, -no_var_idx]
        x_test = x_test[, -no_var_idx]
      }
      
      #Scale Data
      print("\n---------------\nScale Data...\n---------------\n")
      train_mean = apply(x_train, 2, mean)
      train_sd = apply(x_train, 2, sd)
      # Train
      z_x_train = sweep(x_train, 2, train_mean, "-")
      z_x_train = sweep(z_x_train, 2, train_sd, "/")
      # Test
      z_x_test = sweep(x_test, 2, train_mean, "-")
      z_x_test = sweep(z_x_test, 2, train_sd, "/")

      
      # Make y categorical
      y_train = as.factor(y_train)
      y_test = factor(y_test, levels = levels(y_train))
      
      # Run ML for multiple methods
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
    print(summary_df)

    saveRDS(summary_df,paste0(prefix, "_Affy_DE_100runs_ML_summary.RDS"))
  }
}


stopCluster(cluster)
intervalEnd <- Sys.time()
paste("100 iterations of data took",intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))

# aggregate mean and sem for each method and metrics
print(dim(summary_df))

DF <- as.data.frame(summary_df)
means <- aggregate( .~ method, DF, mean)
rownames(means) <- paste0("mean_", means$method)
sds <- aggregate( .~ method, DF, sd)
rownames(sds) <- paste0("sd_", means$method)

tab <- rbind(means, sds)
tab <- tab[, -1]
tab <- t(tab)

print(tab)
