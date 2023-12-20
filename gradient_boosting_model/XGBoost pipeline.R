#LOAD PACKAGES AND DATA ----
##Packages----
library(ParBayesianOptimization)
library(xgboost)
library(Matrix)
library(tidyverse)
library(doParallel)
library(parallel)
library(pROC)
library(Ckmeans.1d.dp)
##Dataset----
load("./Viral_Receptor_Prediction_Dataset.RData") #Modify the address to the location of the dataset file. 
#Data: contains all variables needed for the boosting model.
#Gene_info: contains gene names and ensembl and uniprot codes. 
label <- Data[,1] # Extraction of the objective variable for future steps.

#BAYESIAN OPTIMIZATION ----
##Model Preparation ----
###Optimization function ----
opt_fn <- function(eta, gamma, max_depth, subsample, colsample_bytree, alpha, lambda) {#The arguments of the function are the hyperparameters to be optimized.
  
  #Transforms our dataset (sparse matrix) to xgbDMatrix data structure
  xgbData <- xgboost::xgb.DMatrix(data = Data[,-1], label=label, missing = NA)
  
  #List of hyperparameters
  pars <- list(
    eta = eta,
    gamma = gamma,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    alpha = alpha,
    lambda = lambda
    )
  
  #XGBoost model with cross-validation 
  xgb_model_CV <- xgb.cv(
    nthread=0,
    params = pars,
    data = xgbData,
    early_stopping_rounds = 50, #Early stopping criterion
    nfold = 10, #Ten-fold cross validation
    nrounds = 10000,#Maximum number of xgb model iteration
    scale_pos_weight = sum(label==0)/sum(label==1), #Class imbalance correction
    objective = "binary:logistic",
    eval_metric = "auc",
    verbosity = 0,
    prediction = T,
    showsd = T,
    maximize = T,
    stratified = T
  )
  
  #Definition of output metrics
  nrounds = xgb_model_CV$best_iteration #Number of rounds to reach optimal AUC in the cross validation test subset.
  conf_mat <- as.vector(table(pred=as.integer(xgb_model_CV$pred>0.5),lab=label))#Confusion matrix based on model predictions
  if(length(conf_mat)!=4){conf_mat <- matrix(NA,2,2)}
  precision <- conf_mat[4]/sum(conf_mat[c(2,4)])#Precision = TP/(TP+FP) What fraction of the predicted receptors are real (or known).
  recall <- conf_mat[4]/sum(conf_mat[c(3,4)])#Recall = TP/(FN+TP) What fraction of the known receptors is recovered by the model prediction.
  Score = (precision+recall)/2 #This is our objetive metric to maximize (Average or combined precision and recall)
  auc = max(xgb_model_CV$evaluation_log$test_auc_mean) #Model AUC.
  
  #Function output
  return(
    list(
      Score = Score,
      nrounds = nrounds,
      Precision=precision,
      Recall=recall,
      TN=conf_mat[1],
      FP=conf_mat[2],
      FN=conf_mat[3],
      TP=conf_mat[4],
      AUC=auc
    )
  )
}

###Boundaries for sampling our hyperparameters ----
boundaries <- list(
  eta = c(0.01, 0.5),
  gamma = c(0, 10),
  max_depth = c(1L, 10L),
  subsample = c(0.1, 1),
  colsample_bytree = c(0.1,1),
  alpha = c(0,20),
  lambda = c(0,20)
)

##Parallel Running of Bayesian Optimization ----
no_cores <- detectCores() #Or select the desired amount of cores to use.  
cl <- makeCluster(no_cores) #Cluster creation
registerDoParallel(cl)
clusterExport(cl, c('Data','label')) #Export data to clusters
clusterEvalQ(cl,expr= { #Load xgboost library to use in our parallelized function
  library(xgboost)
})

#Run our Bayesian Optimization
opt_bayes <- bayesOpt(
      FUN = opt_fn,
      bounds = boundaries,
      initPoints = 100, #Number of initial parameter sets sampled 
      iters.n = 500, #Total number of xgb models run apart from the 100 initial ones. 
      iters.k = 10, #Of these total models, 10 will be run per each refinement epoch.
      parallel = TRUE
    )

stopCluster(cl) #End clusters
registerDoSEQ()

#SELECTION OF THE BEST MODEL ON AVERAGE ----
##Obtaining a list of top-ranked hyperparamenter sets ----
N <- 50 #Set the number of top-ranked hyperparameter sets chosen. 
replicates <- 100 #Set the number of replicates per hyperparameter set
Parms_grid <- opt_bayes$scoreSummary %>% arrange(desc(Score))%>%
  select(eta, gamma, max_depth, subsample, colsample_bytree,alpha, lambda) %>%
  slice(1:N) %>% slice(rep(1:n(), each = replicates))

##Definition of a XGBModel function ----
XGBModel=function(iter){
  #Transforms our dataset (sparse matrix) to xgbDMatrix data structure
  xgbData <- xgboost::xgb.DMatrix(data = Data[,-1], label=label, missing = NA)
  
  #Obtain a list of the hyperparameter of a given set in Parms_grid 
  l <- Parms_grid[iter,]
  l <- lapply(l,function(x)list(x))
  
  #Fit xgb model
  xgb_model_CV <- xgb.cv(data = xgbData, nrounds = 10000, nthread = 0,missing = NA,params = l,scale_pos_weight = sum(label==0)/sum(label==1),
                   metrics = "auc",early_stopping_rounds = 50,nfold = 10,verbose = F,
                   objective = "binary:logistic",prediction = T,showsd = T,stratified = T,maximize = T)
  #Return Output
  return(xgb_model_CV)
}

##Running XGBModel replicates in parallel ----
iter <- 1:(N*replicates) #Number of iterations in parallel
no_cores <- detectCores() #Or select the desired amount of cores to use.  
cl <- makeCluster(no_cores) #Cluster creation
registerDoParallel(cl)
clusterExport(cl,varlist = c("Data","label","Parms_grid")) #Export data to clusters
clusterEvalQ(cl,expr= { #Load xgboost library to use in our parallelized function
  library(xgboost)
})
z <- parLapply(cl = cl,X = iter,fun = XGBModel)#Running our function in parallel. Outputs combined in a list.
stopCluster(cl) #End clusters
registerDoSEQ()

##Obtaining XGBModel logs, predictions and output metrics ----
Logs <- list()
Preds <- list()
Output <- data.frame()

for(i in 1:length(z)){
  Preds[[i]] <- z[[i]]$pred
  Logs[[i]] <- z[[i]]$evaluation_log
  Output <- rbind(Output,cbind(Parms_grid[i,],z[[i]]$evaluation_log[z[[i]]$best_iteration,]))
  print(i)}
#Create a dataframe with metrics (TN, FP, FN, TP, Precision, Recall, F1, Mean_PR, AUC) for each hyperparamenter set
XGBMetrics <- data.frame()
for (i in 1:length(Preds)){
  conf_mat <- as.vector(table(pred=as.integer(Preds[[i]]>0.5),lab=label))#Confusion matrix based on model predictions
  if(length(conf_mat)!=4){next}else{
    precision <- conf_mat[4]/sum(conf_mat[c(2,4)])#Precision = TP/(TP+FP) What fraction of the predicted receptors are real (or known).
    recall <- conf_mat[4]/sum(conf_mat[c(3,4)]) #Recall = TP/(FN+TP) What fraction of the known receptors is recovered by the model prediction.
    f1 <- 2/(1/precision+1/recall) #F1 score: harmonic mean of precision and recall
    mean_precision_recall <- (precision+recall)/2 #Combined precision and recall (average)
    XGBMetrics[i,1:8] <- c(conf_mat,precision,recall,f1,mean_precision_recall)}}
colnames(XGBMetrics) <- c("TN","FP","FN","TP","Precision","Recall","F1","Mean_PR")#Set data frame column names.

ModelID <- rep(1:N,each=replicates) #ModelID to differentiate replicates of each hyperparameter set.
#Add ModelID and cross validation test auc mean to the XGBMetrics dataframe.
XGBMetrics <- cbind(ModelID,Output[,c(1:8)],XGBMetrics,AUC=Output$test_auc_mean)
#View(XGBMetrics) #Run to visualize the output

##Creating a data frame with mean and variance of metrics for each ModelID ----
XGBM1 <- XGBMetrics %>% group_by(ModelID)%>%summarise(across(1:ncol(Parms_grid),mean))
XGBM2 <- XGBMetrics%>%group_by(ModelID)%>%summarise(across((ncol(Parms_grid)+1):(ncol(XGBMetrics)-1),list(avg=mean,var=var)))
XGBMeanMetrics <- cbind(XGBM1,XGBM2[,-1])
#View(XGBMeanMetrics) #Run to visualize the output

##Deciding the best model ----
#Best model based on meanPR (Average precision and recall)
test_mu <- data.frame(mu_mPR=numeric(),model = numeric(),p_value= numeric())
Best_model <- data.frame()
min_mu <- quantile(XGBMetrics$Mean_PR,probs = 0.1)#Probs can be modified if needed
max_mu <- quantile(XGBMetrics$Mean_PR,probs = 0.9)#Probs can be modified if needed
delta=0.0001 #Step size between means
mu_seq <- seq(min_mu,max_mu,delta)#Fixed means (mu) to test in one-sample t-test
for(j in 1:length(mu_seq)){
  for(i in 1:N){
    ttest <- XGBMetrics%>%filter(ModelID==i)%>%pull(Mean_PR)%>%t.test(mu= mu_seq[j],alternative = "greater")
    test_mu[i,] <-c(mu_seq[j],i,ttest$p.value)}
    Best_model <- rbind(Best_model,test_mu)
}

# View(Best_model) #Run to see dataframe

#Plot to visualize models by meanPR significanly greater than threshold (y axis)
plot_meanPR <- Best_model%>%group_by(model)%>%summarise(th_meanPR=mu_seq[which.min(p_value<0.01)])%>%ggplot(aes(y=th_meanPR,x=factor(model)))+geom_bar(stat = "identity",width = 0.8,col="black",fill="green4")+coord_cartesian(ylim=c(min_mu,max_mu))+theme_bw()+
  labs(x="Model",y="Significance threshold meanPR (t-test)")
plot_meanPR

#Find and save the best model ID
Best_model_ID <- Best_model%>%filter(p_value<0.01)%>%filter(mu_mPR==max(mu_mPR))%>%pull(model) #If more than one model appears reduce delta in previous steps (line 188). 
Best_model_ID

##Average best model prediction and metrics ----
preds <- unlist(Preds) #Unlist N*replicates cross validation model predictions.
dim(preds) <- c(nrow(Data),N*replicates) #Give them proper dimensions.
preds <- as.data.frame(cbind(ModelID,t(preds))) #Generate a data frame with model ID and predictions per row.
mean_preds <- preds%>%group_by(ModelID)%>%summarise(across(everything(),mean)) #Calculate averaged predictions for each model.
Avg_Pred <- cbind(Gene_info,Avg_model_pred=t(mean_preds[Best_model_ID,-1]),label) #Save averaged predictions for the best model.

#Calculate metrics based on averaged best model predictions.
conf_mat <- as.vector(table(Avg_Pred$Avg_model_pred>0.5,Avg_Pred$label))#Confusion matrix based on averaged best model predictions
precision <- conf_mat[4]/sum(conf_mat[c(2,4)])#Precision = TP/(TP+FP) What fraction of the predicted receptors are real (or known).
recall <- conf_mat[4]/sum(conf_mat[c(3,4)])#Recall = TP/(FN+TP) What fraction of the known receptors is recovered by the average best model prediction.
f1 <- 2/(1/precision+1/recall) #F1 score: harmonic mean of precision and recall
Mean_PR <-  (precision+recall)/2 #Averaged precision and recall
AUC <-  pROC::auc(label,Avg_Pred$Avg_model_pred) #AUC 

#Best Model metrics based on averaged cross validation predictions
Avg_best_model_metrics <- data.frame(XGBM1[Best_model_ID,], TN = conf_mat[1], FP = conf_mat[2], FN = conf_mat[3], TP = conf_mat[4],
                                     Precision = precision, Recall = recall, F1 = f1, Mean_PR = Mean_PR, AUC = AUC);Avg_best_model_metrics

#TRAINING XGBOOST MODEL AND GETTING VARIABLE (FEATURE) IMPORTANCE ----

##Train XGBoost Model with best parameter combination and average nrounds (iterations) ----

l <- Parms_grid[Best_model_ID,] #Get best hyperparameter combination
l <- lapply(l,function(x)list(x)) #List of hyperparameters

set.seed(42)
xgbData <- xgboost::xgb.DMatrix(data = Data[,-1], label=label, missing = NA)
XGBOOST_Model <- xgboost(data = xgbData,missing = NA,params = l,nrounds = XGBMeanMetrics%>%filter(ModelID==Best_model_ID)%>%pull(iter_avg)%>%round(),verbose = T,scale_pos_weight = sum(label==0)/sum(label==1))
##Getting variable importance ----
importance <- xgb.importance(model = XGBOOST_Model)
#View(importance)#Run to see dataframe
#Plot of variable importance
xgb.ggplot.importance(importance)

#END----




