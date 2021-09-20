# Sharvari Gujja
# 
# generate DL models for 100 cuts of train and test
# 
# Usage: Rscript <script> <drugname> <outputDir>
# Example: Rscript DNN_h2o_binary_classf.R Gemcitabine TCGA_Gemcitabine_100cuts/
# Output: Summary table of classification metrics across 100 runs


rm(list = ls())
gc()

library(h2o)
library(skimr)
library(data.table)


args <- commandArgs(TRUE)
options(stringsAsFactors = FALSE)
##################################################################
# input arguments
# select drug name
prefix <- args[1] #"Gemcitabine" #
# output dir
out <- args[2] #"/Users/sgujja/repos/drug_R_NR_modeling//TCGA/scripts/test_output/" # 

# Stats Mat
summary_df=data.frame(train_acc=double(),test_acc=double(),
                      train_auroc=double(), test_auroc=double(),
                      train_prec=double(), test_prec=double(), 
                      train_recall=double(), test_recall=double(),
                      train_F1=double(),test_F1=double(),
                      train_Specificity=double(),test_Specificity=double()
                      )



  
n=1
for( cut in c(1:100)){
  
  if (file.exists(paste0(out,"/",prefix, '_Affy_DE_TrainingSet_new_',cut,'.csv'))){
    train = read.csv(paste0(out,"/",prefix, '_Affy_DE_TrainingSet_new_',cut,'.csv'), row.names=1, header=T, check.names = FALSE)
    test = read.csv(paste0(out,"/",prefix, '_Affy_DE_TestingSet_new_',cut,'.csv'), row.names=1, header=T, check.names = FALSE)
  
    
    train$group <- ifelse(train[,"group"]=="S", 0, 1)
    test$group <- ifelse(test[,"group"]=="S", 0, 1)
    
    # set group column as factor
    train$group<- as.factor(train$group)
    test$group <- as.factor(test$group)
    
    setDT(train)
    setDT(test)
    
    # print dimensions
    print(dim(train))
    print(dim(test))
    
    #check target variable
    #binary in nature check if data is imbalanced
    train[,.N/nrow(train),group]
    test[,.N/nrow(test),group]
    
    if (ncol(train) >= 4) {
      y_train = train[[2]]
      y_train = as.factor(y_train)
      x_train = train[,-c(1,2)]
      
      y_test = test[[2]]
      y_test = as.factor(y_test)
      x_test = test[,-c(1,2)]

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
      
      
      
      train <- as.data.frame(cbind(y_train,z_x_train))
      test <- as.data.frame(cbind(y_test,z_x_test))
      
      # reassign column name
      names(train)[names(train) == 'y_train'] <- 'group'
      names(test)[names(test) == 'y_test'] <- 'group'
      
      localH2o <- h2o.init(nthreads = -1, max_mem_size = "20G",port = 54321+(cut*2))
      
      #load data on H2o
      trainh2o <- as.h2o(train)
      testh2o <- as.h2o(test)
      
      
      #set variables
      y <- "group"
      x <- setdiff(colnames(trainh2o),y)
      
      ### Random Grid Search
      
      hyper_params <- list(
        activation=c("Rectifier","Tanh","Maxout","RectifierWithDropout","TanhWithDropout","MaxoutWithDropout"),
        hidden=lapply(1:100, function(x)10+sample(50,sample(4), replace=TRUE)),
        input_dropout_ratio=c(0,0.05),
        l1=seq(0,1e-3,1e-5),
        l2=seq(0,1e-3,1e-3)
      )
      
      #set search criteria
      
      search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 900, max_models = 100, seed=1234567, stopping_rounds=15, stopping_tolerance=0.001)
      #train model
      
      model_grid <- h2o.grid("deeplearning",
                             grid_id = paste0("mygrid_",cut),
                             hyper_params = hyper_params,
                             search_criteria = search_criteria,
                             x = x,
                             y = y,
                             training_frame = trainh2o,
                             validation_frame = testh2o,
                             nfolds = 5,
                             epochs = 1000,
                             balance_classes = TRUE,
                             fold_assignment = "Stratified",
                             score_interval = 2,
                             stopping_metric = "misclassification"
      )
      
      #get best model
      
      d_grid <- h2o.getGrid(paste0("mygrid_",cut),sort_by = "auc",decreasing=TRUE)
      best_dl_model <- h2o.getModel(d_grid@model_ids[[1]])
      
      # save model
      path <- h2o.saveModel(best_dl_model, path=paste0(out,"/",prefix,"./mybest_deeplearning_model_",cut), force=TRUE)
      
      # Retrieve the variable importance
      print(h2o.varimp(best_dl_model))
      #compute variable importance and performance
      h2o.varimp_plot(best_dl_model,num_of_features = 20)
      
      # metrics
      train_perf <- h2o.performance(best_dl_model) # train metrics
      test_perf <- h2o.performance(best_dl_model, testh2o)  # test metrics
      
      #head(h2o.metric(test_perf))
      print("Confusion Matrix: Train:")
      print(h2o.confusionMatrix(train_perf))
      print("Confusion Matrix: Test:")
      print(h2o.confusionMatrix(test_perf))
      
      # Classify the test set (predict class labels)
      # This also returns the probability for each class
      pred <- h2o.predict(best_dl_model, newdata = testh2o)
      # Take a look at the predictions
      print(head(pred))
      
      ##
      
      fpr = best_dl_model@model$training_metrics@metrics$thresholds_and_metric_scores$fpr
      tpr = best_dl_model@model$training_metrics@metrics$thresholds_and_metric_scores$tpr
      fpr_val = best_dl_model@model$validation_metrics@metrics$thresholds_and_metric_scores$fpr
      tpr_val = best_dl_model@model$validation_metrics@metrics$thresholds_and_metric_scores$tpr
      plot(fpr,tpr, type='l')
      title('AUC')
      lines(fpr_val,tpr_val,type='l',col='red')
      legend("bottomright",c("Train", "Test"),col=c("black","red"),lty=c(1,1),lwd=c(3,3))
      
      training_metrics = data.frame(row.names= best_dl_model@model$training_metrics@metrics$max_criteria_and_metric_scores$metric , 
                                    val=best_dl_model@model$training_metrics@metrics$max_criteria_and_metric_scores$value)
      
      
      testing_metrics = data.frame(row.names= best_dl_model@model$validation_metrics@metrics$max_criteria_and_metric_scores$metric , 
                                   val=best_dl_model@model$validation_metrics@metrics$max_criteria_and_metric_scores$value)
      
      
      summary_df[n,'train_auroc']=h2o.auc(best_dl_model, valid=FALSE) # train
      summary_df[n,'test_auroc']=h2o.auc(best_dl_model, valid=TRUE) # test
      summary_df[n,'train_acc']= training_metrics["max accuracy",]
      summary_df[n,'test_acc']=testing_metrics["max accuracy",]
      summary_df[n,'train_prec']=training_metrics["max precision",]
      summary_df[n,'test_prec']=testing_metrics["max precision",]
      summary_df[n,'train_recall']=training_metrics["max recall",]
      summary_df[n,'test_recall']=testing_metrics["max recall",]
      summary_df[n,'train_F1']=training_metrics["max f1",]
      summary_df[n,'test_F1']=testing_metrics["max f1",]
      summary_df[n,'train_Specificity']=training_metrics["max specificity",]
      summary_df[n,'test_Specificity']=testing_metrics["max specificity",]
      
      n = n+1

    }
  }
}

# aggregate mean and sd for each method and metrics
tab <- skim(summary_df)

print(dim(summary_df))
write.csv(summary_df,paste0(out,"/",prefix,"_DL_metrics_info.csv"),quote=F,row.names=F )
write.csv(as.data.frame(tab),paste0(out,"/",prefix,"_DL_metrics_info_mean_SD.csv"),quote=F,row.names=F )

#shutdown H2O
h2o.shutdown(prompt=FALSE)

  
