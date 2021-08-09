library(h2o)
h2o.init()

# Import the insurance dataset into H2O:
train <- read.csv("/Users/sgujja/repos/drug_R_NR_modeling//TCGA/scripts/test_output/Gemcitabine_DE_TrainingSet_new_1.csv",sep=",", header=T,row.names = 1,stringsAsFactors = F, check.names = F )[,-1]
test <- read.csv("/Users/sgujja/repos/drug_R_NR_modeling//TCGA/scripts/test_output/Gemcitabine_DE_TestingSet_new_1.csv",sep=",", header=T, row.names = 1 ,stringsAsFactors = F, check.names = F)[,-1]

# set group column as factor
train$group <- as.factor(train$group)
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

localH2o <- h2o.init(nthreads = -1, max_mem_size = "20G")

#load data on H2o
trainh2o <- as.h2o(train)
testh2o <- as.h2o(test)


#set variables
y <- "group"
x <- setdiff(colnames(trainh2o),y)


### Random Grid Search
activation_opt <- c("Rectifier","RectifierWithDropout", "Maxout","MaxoutWithDropout")
hidden_opt = lapply(1:100, function(x)10+sample(50,sample(4), replace=TRUE))
l1_opt <- c(0,1e-3,1e-5)
l2_opt <- c(0,1e-3,1e-5)

hyper_params <- list( activation=activation_opt,
                      hidden=hidden_opt,
                      l1=l1_opt,
                      l2=l2_opt )

#set search criteria
search_criteria = list(strategy = "RandomDiscrete",
                       max_models = 10, max_runtime_secs = 1000,
                       seed=123456)

#train model
model_grid <- h2o.grid("deeplearning",
                       grid_id = "mygrid",
                       hyper_params = hyper_params,
                       search_criteria = search_criteria,
                       x = x,
                       y = y,
                       training_frame = trainh2o,
                       nfolds = 10,
                       balance_classes = TRUE,
                       score_interval = 2,
                       epochs = 10000,
                       stopping_rounds = 3,
                       stopping_tolerance = 0.05,
                       stopping_metric = "misclassification")


#get best model
d_grid <- h2o.getGrid("mygrid",sort_by = "accuracy")
best_dl_model <- h2o.getModel(d_grid@model_ids[[1]])


# Retrieve the variable importance
h2o.varimp(best_dl_model)
#compute variable importance and performance
h2o.varimp_plot(best_dl_model,num_of_features = 20)

# metrics
train_perf <- h2o.performance(best_dl_model) # train metrics
test_perf <- h2o.performance(best_dl_model, testh2o)  # test metrics

#head(h2o.metric(test_perf))
h2o.confusionMatrix(train_perf)
h2o.confusionMatrix(test_perf)
h2o.accuracy(test_perf)
plot(test_perf, type = "roc")

# Classify the test set (predict class labels)
# This also returns the probability for each class
pred <- h2o.predict(best_dl_model, newdata = testh2o)

# Take a look at the predictions
print(head(pred))


##
h2o.auc(best_dl_model, train = TRUE)
h2o.auc(best_dl_model, xval = TRUE)

#shutdown H2O
h2o.shutdown(prompt=FALSE)


