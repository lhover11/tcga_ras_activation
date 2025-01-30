# Purpose: perform elastic net logistic regression

# load libraries
library(tidyverse)
library(edgeR)
library(data.table)
library(caret)
library(glmnet)
library(pROC)
library(ggrepel)

# read in data
meta <- read.csv("../processed/TCGA_PAAD_meta_for_ML.csv")
data <- read.csv("../processed/TCGA_PAAD_top8000_var_for_ML.csv")
genes <- read.table("../raw-data/PAAD/gencode.v36.annotation.gtf.gene.probemap", skip = 1)

## combine data
data <- meta %>% inner_join(data, by = c("sampleid"))
dim(data)

order <- data$sampleid

## remove sample id from data
row.names(data) <- data$sampleid
data$sampleid <- NULL
dim(data)

# split the data into training and testing sets (75/25 split)
set.seed(25)
trainIndex <- createDataPartition(data$ras_event, p = .75, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

## convert to matrix input and remove intercept term
x_train <- model.matrix(ras_event ~ ., trainData)[, -1]
y_train <- trainData$ras_event


# elastic net logistic regression with 5-fold cross validation
cv_fit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5, # 0.5 for elastic net
  type.measure = "auc", # measure AUROC
  nfolds = 5, # 5-fold cross-validation
  stratified = TRUE # ensure balanced folds
)

plot(cv_fit)
cv_fit

## Summary of AUCs from cross-validation
summary(cv_fit$cvm) # max AUC = 0.89


# Evaluate model on the test data
## get into matrix format, remove intercept term
x_test <- model.matrix(ras_event ~ ., testData)[, -1]
y_test <- testData$ras_event

# Make predictions on the test set
predictions <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")
# predictions <- predict(cv_fit, newx = x_test, s = "lambda.1se", type = "response") # with lambda 1se, a lot fewer features are included with not a huge change in AUC


# Calculate AUC on the test set
roc_obj <- roc(y_test, predictions)
auc_value <- auc(roc_obj)
print(paste("Test AUROC:", auc_value)) # 0.736

# Plot AUC for training and testing data
#### ROC curve for the training data (cross-validation predictions)
train_predictions <- predict(cv_fit, newx = x_train, s = "lambda.min", type = "response")
roc_train <- roc(y_train, train_predictions)
auc(roc_train) # 1 - perfect model on training data, the AUC for the cross-validation data is a more reliable measure of how this model will work on unseen data

df <- data.frame(x = train_predictions, y = y_train, thres = roc_train$thresholds)


### Plot both ROC curves
ggroc(list(train = roc_train, test = roc_obj), size = 1, legacy.axes = TRUE) +  # use legacy.axes to plot the FPR instead of the specificity
  scale_colour_manual("Data", values = c("blue", "darkgreen")) +
  geom_segment(aes(x=0, y=0, xend=1, yend=1), color='black') +
  xlab("FPR") +
  ylab("TPR") +
  theme_classic() +
  theme(text=element_text(size=20))
ggsave("../results/plots/PAAD_test_train_ROC_curve_lamba_min.png", w = 6, h = 6)    


# Find the optimal threshold based on Youden's Index
optimal_threshold <- coords(roc_train, "best", ret = "threshold", best.method = "youden")

class_predictions <- ifelse(predictions > optimal_threshold$threshold[1], 1, 0)
confusion_mat <- confusionMatrix(as.factor(class_predictions), as.factor(y_test), positive = "1")
confusion_mat
## accuracy = 0.641 (0.4718-0.788)
## Specificity = 0.375 - True negative rate
## Sensitivity = 0.8261 - True positive rate, model is better at detecting RAS activated samples (model correctly predicted 19/23 +samples)


# Get coefficients 
feat <- as.data.frame(as.matrix(coef(cv_fit, s = "lambda.min")))
range(feat$s1)

## add in gene names
feat$gene_id <- row.names(feat)
feat <- feat %>% left_join(genes, by = c("gene_id" = "V1"))

## remove intercept
feat <- feat %>% filter(!grepl("Intercept", gene_id))

## add in ranking
feat <- feat %>%
  arrange(s1) %>%
  mutate(rank = 1:nrow(feat))

# Label the top 10 and bottom 10 values
feat$label <- ifelse(feat$rank <= 10, feat$V2,
                     ifelse(feat$rank > (nrow(feat) - 10), feat$V2, NA))

## color top and bottom ranking genes
feat$color <- ifelse(!is.na(feat$label), "top", "not")

ggplot(feat, aes(x = rank, y = s1, label = label, color = color)) +
  geom_point() +
  scale_color_manual(values = c("gray", "red")) +
  geom_text_repel(color = "black") +
  xlab("Ranking") +
  ylab("Lambda min coefficient") +
  ggtitle("PAAD RAS classifier coefficients") +
  theme_minimal() +
  theme(legend.position = "none", 
        text=element_text(size=20))
ggsave("../results/plots/PAAD_genes_impacting_classification_lamba_min.png", w = 6, h = 6)    
