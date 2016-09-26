######################################################################################################################
###
###  PREDICTIVE MODEL FOR SCORING LONG-TERM SURVIVABILITY OF PROSTATE CANCER PATIENTS
###
###  AUTHOR:   Cody Frehr
###  DATE:     7/28/16
###  FILES:    CodyFrehr_ScoringModel.R
###            CodyFrehr_ScoringAnalysis.R
###            training_data.csv
###            CodyFrehr_score.csv
###
###  OUTLINE:  1) Heuristics derivation
###            2) Cross-validation of Bayesian classifier and regression models
###            3) Learned best weights in ensemble method
###            4) Supplementary deterministic (common sense) adjustments
###
######################################################################################################################


### LOAD PACKAGES ####################################################################################################

library(Amelia)
library(arules)
library(caret)
library(e1071)
library(ggplot2)
library(lattice)
library(scatterplot3d)


### READ AND PROCESS DATA ############################################################################################

# read data
setwd("/Users/cfrehr/Job_Search/Interviews/Enova/Interview_Assignment/participant_files/")
patient_data <- read.csv(file="training_data.csv", header=TRUE, sep=",")
score_data <- read.csv(file="CodyFrehr_score.csv", header=TRUE, sep=",")

# preprocess patient_data data types
patient_data$t_score <- ordered(patient_data$t_score, levels(patient_data$t_score)[c(1:10)])
patient_data$m_score <- ordered(patient_data$m_score, levels(patient_data$m_score)[c(1:4)])
patient_data$stage <- ordered(patient_data$stage, levels(patient_data$stage)[c(1:5)])
patient_data$race <- as.factor(patient_data$race)

# preprocess score_data data types
score_data$t_score <- ordered(score_data$t_score, levels(score_data$t_score)[c(1:10)])
score_data$m_score <- ordered(score_data$m_score, levels(score_data$m_score)[c(1:4)])
score_data$stage <- ordered(score_data$stage, levels(score_data$stage)[c(1:5)])
score_data$race <- as.factor(score_data$race)


### DERIVE SURVIVABILITY HEURISTICS (age_group, count_thrpy, change_tumor, change_psa) ###############################
 
prop_data <- patient_data[c("age", "rd_thrpy", "h_thrpy", "chm_thrpy", "cry_thrpy", "brch_thrpy", "rad_rem")]

prop_data["age_group"] <- NA
prop_data["count_thrpy"] <- NA
prop_data["change_tumor"] <- NA
prop_data["change_psa"] <- NA
prop_data["age_group"] <- as.numeric(prop_data$age_group)
prop_data["count_thrpy"] <- as.numeric(prop_data$count_thrpy)
prop_data["change_tumor"] <- as.numeric(prop_data$change_tumor)
prop_data["change_psa"] <- as.numeric(prop_data$change_psa)

prop_data["rd_thrpy"] <- as.numeric(prop_data$rd_thrpy)
prop_data["h_thrpy"] <- as.numeric(prop_data$h_thrpy)
prop_data["chm_thrpy"] <- as.numeric(prop_data$chm_thrpy)
prop_data["cry_thrpy"] <- as.numeric(prop_data$cry_thrpy)
prop_data["brch_thrpy"] <- as.numeric(prop_data$brch_thrpy)
prop_data["rad_rem"] <- as.numeric(prop_data$rad_rem)

# calculate heuristics
for (i in 1:dim(prop_data)[1]) {
    
    # store count_thrpy
    prop_data[i,9] <- prop_data[i,2] + prop_data[i,3] + prop_data[i,4] + prop_data[i,5] + prop_data[i,6] + prop_data[i,7]
    
    # bin ages in age_group
    temp_age <- prop_data[i,1]
    if (!is.na(temp_age)) {
        if (temp_age < 50) {
            prop_data[i,8] <- 1
        } else if (temp_age < 60) {
            prop_data[i,8] <- 2
        } else if (temp_age < 70) {
            prop_data[i,8] <- 3
        } else if (temp_age < 80) {
            prop_data[i,8] <- 4
        } else if (temp_age < 90) {
            prop_data[i,8] <- 5
        } else if (temp_age < 100) {
            prop_data[i,8] <- 6
        } else {
            prop_data[i,8] <- 7
        }
    }
}

# modify data types and unify with existing data 
prop_data[,"age_group"] <- as.factor(prop_data$age_group)
prop_data[,"count_thrpy"] <- as.factor(prop_data$count_thrpy)
patient_data[,12] <- prop_data$age_group
#patient_data[,13] <- prop_data$count_thrpy
names(patient_data)[12]<-paste("age_group")
#names(patient_data)[13]<-paste("count_thrpy")


### FOLD DATA FOR CROSS-VALIDATION ###################################################################################

folds <- createFolds(patient_data$survival_7_years)
folded_data <- lapply(folds, function(ind, dat) dat[ind,], dat = patient_data)

###

solution <- patient_data$survival_7_years
solution[1] <- 10000

######################################################################################################################
### TRAINING PREDICTIVE MODELS
### GOALS:  (1) Naive Bayesian Classifier (gleason_score, t_score, n_score, m_score, stage)
###         (2) Logit Regression Model of Single Input (tumor_x)
###             (i)   Choose tumor_diagnosis, impute if missing
###             (ii)  Choose tumor_1_year, impute if missing
###             (iii) Choose tumor_change, reflective of change over time
###         (3) Logit Regression Model of Single Input (psa_x)
###             (i)   Choose psa_diagnosis, impute if missing
###             (ii)  Choose psa_1_year, impute if missing
###             (iii) Choose best available psa_change
### 
### DETERMINISTIC IMPUTATIONS
###         (1) NA analysis - high accuracy in predicting death
###         (2) 0 tumor/psa - 100% life expectency
###
### 
###
######################################################################################################################


### (1) Naive Bayesian Classifier (gleason_score, t_score, n_score, m_score, stage, age_group, side, count_thrpy) ####

# cross-validation accuracies
a1i <- c(0,0,0,0,0,0,0,0,0,0)



for (i in 1:10) {
    
    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # adjust categorical data types for Bayes
    valid[,"survival_7_years"] <- as.factor(valid$survival_7_years)
    valid[,"gleason_score"] <- ordered(valid$gleason_score)
    valid[,"age_group"] <- ordered(valid$age_group)
    
    #train[,32] <- NULL
    train[,"gleason_score"] <- ordered(train$gleason_score)
    train[,"t_score"] <- ordered(train$t_score)
    #train[,"count_therpy"] <- ordered(train$count_therpy)
    levels(train[,"t_score"]) <- c("T1a", "T1b", "T1c", "T2a", "T2b", "T2c", "T3a", "T3b", "T3c", "T4")
    train[,"n_score"] <- as.factor(train$n_score)
    levels(train[,"n_score"]) <- c("N0", "N1", "NX")
    train[,"side"] <- as.factor(train$side)
    levels(train[,"side"]) <- c("both", "left", "right")
    train[,"m_score"] <- ordered(train$m_score)
    levels(train[,"m_score"]) <- c("M0", "M1a", "M1b", "M1c")
    train[,"stage"] <- ordered(train$stage)
    levels(train[,"stage"]) <- c("I", "IIA", "IIB", "III", "IV")
    train[,"age_group"] <- ordered(train$age_group)
    train[,"survival_7_years"] <- as.factor(train$survival_7_years)
   
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    count <- 0 # count of correctly classified instances
    
    # train and test naive Bayes classifier
    bayesModel <- naiveBayes(survival_7_years ~ gleason_score + t_score + m_score + n_score + stage, data=train)
    surv[,"pred"] <- predict(bayesModel, valid[1:length_valid,])
    
    # evaluate and record accuracy of model
    for (j in 1:length_valid) {
        if (!is.na(surv[j,"pred"])) {
            if (!is.na(surv[j,"true"])) {
                if (surv[j,"pred"] == surv[j,"true"]) {
                    count <- count +1
                }
            }
        }  
    }
    
    a1i[i] <- count / length_valid
}

# final accuracy
Accuracy1i <- sum(a1i)/10


### (2) Logit Regression Model of Single Input (tumor_x) #############################################################
# for prediction intervals: http://stats.stackexchange.com/questions/26650/how-do-i-reference-a-regression-models-coefficients-standard-errors

### (2i) Choose tumor_diagnosis, impute if missing ###################################################################

# psuedo:
# create a vector to store model accuracies
# for each fold
#   copy true survival_7_years values
#   find best fit linear regression (tumor_diagnosis ~ tumor_6_months)
#   find best fit linear regression (tumor_diagnosis ~ tumor_1_year)
#   find best fit linear regression (tumor_diagnosis ~ tumor_6_months + tumor_1_year)
#   find best fit logit regression (survival_7_years ~ tumor_diagnosis)
#   for each case in fold
#     if tumor_diagnosis != NA, plug into logit regression
#       save predicted survival value
#     else, impute tumor_diagnosis
#       plug available tumor values into corresponding regression model
#       plug expected tumor_diagnosis into logit regression
#       save predicted survival value
#     compute and store accuracy of predictions

# cross-validation accuracies
a2i <- c(0,0,0,0,0,0,0,0,0,0)

# train and test 10 models
for (i in 1:10) {

    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"logit"] <- NA # logit values returned from classification of validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    coefs <- data.frame(matrix(ncol=4, nrow=3)) # placeholder for regression model coefficients
    count <- 0 # count of correctly classified instances
    
    # train regression models
    temp_data <- train[c("tumor_diagnosis", "tumor_6_months")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model1 = lm(formula = tumor_diagnosis ~ tumor_6_months, data=temp_data)
    coefs[1:length(coef(model1)),1] <- coef(model1)
    
    temp_data <- train[c("tumor_diagnosis", "tumor_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model2 = lm(formula = tumor_diagnosis ~ tumor_1_year, data=temp_data)
    coefs[1:length(coef(model2)),2] <- coef(model2)
    
    temp_data <- train[c("tumor_diagnosis", "tumor_6_months", "tumor_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model3 = lm(formula = tumor_diagnosis ~ tumor_6_months + tumor_1_year, data=temp_data)
    coefs[1:length(coef(model3)),3] <- coef(model3)
    
    temp_data <- train[c("survival_7_years", "tumor_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model4 = glm(formula = survival_7_years ~ tumor_diagnosis, data=temp_data, family="binomial")
    coefs[1:length(coef(model4)),4] <- coef(model4)
    meanTumor <- mean(temp_data$tumor_diagnosis)
    
    # test classification model
    for (j in 1:length_valid) {
        
        # fields
        td <- valid[j,"tumor_diagnosis"]
        t6 <- valid[j,"tumor_6_months"]
        t1 <- valid[j,"tumor_1_year"]
        
        # conditional block to handle (potential) imputation and classification
        if (!is.na(td)) {
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_diagnosis=td),)
        } else if (!is.na(t6)) {
            if (!is.na(t1)) {
                td <- predict(model3, newdata=data.frame(tumor_6_months=t6, tumor_1_year=t1))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_diagnosis=td))
            } else {
                td <- predict(model1, newdata=data.frame(tumor_6_months=t6))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_diagnosis=td))
            }
        } else if (!is.na(t1)) {
            td <- predict(model2, newdata=data.frame(tumor_1_year=t1))
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_diagnosis=td))
        } else {
            # guess mean value if all tumor_x are NA
            td <- meanTumor
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_diagnosis=td))
        }
    }

    # evaluate and record accuracy of model
    for (j in 1:length_valid) {
        # logit -> prob:
        # logit = L
        # corresponding probability = exp(L)/(1+exp(L))
        L = surv[j,"logit"]
        if (L > 0.0) {
            surv[j,"pred"] <- 1
        } else {
            surv[j,"pred"] <- 0
        }
        if (surv[j,"pred"] == surv[j,"true"]) {
            count <- count + 1
        }
    }
    a2i[i] <- count / length_valid
}

# final accuracy
Accuracy2i <- sum(a2i)/10


### (2ii)  Choose tumor_1_year, impute if missing ####################################################################

# cross-validation accuracies
a2ii <- c(0,0,0,0,0,0,0,0,0,0)

# train and test 10 models
for (i in 1:10) {
    
    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"logit"] <- NA # logit values returned from classification of validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    coefs <- data.frame(matrix(ncol=4, nrow=3)) # placeholder for regression model coefficients
    count <- 0 # count of correctly classified instances
    
    # train regression models
    temp_data <- train[c("tumor_1_year", "tumor_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model1 = lm(formula = tumor_1_year ~ tumor_diagnosis, data=temp_data)
    coefs[1:length(coef(model1)),1] <- coef(model1)
    
    temp_data <- train[c("tumor_1_year", "tumor_6_months")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model2 = lm(formula = tumor_1_year ~ tumor_6_months, data=temp_data)
    coefs[1:length(coef(model2)),2] <- coef(model2)
    
    temp_data <- train[c("tumor_1_year", "tumor_6_months", "tumor_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model3 = lm(formula = tumor_1_year ~ tumor_6_months + tumor_diagnosis, data=temp_data)
    coefs[1:length(coef(model3)),3] <- coef(model3)
    
    temp_data <- train[c("survival_7_years", "tumor_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model4 = glm(formula = survival_7_years ~ tumor_1_year, data=temp_data, family="binomial")
    coefs[1:length(coef(model4)),4] <- coef(model4)
    meanTumor <- mean(temp_data$tumor_1_year)
    
    # test classification model
    for (j in 1:length_valid) {
        
        # fields
        td <- valid[j,"tumor_diagnosis"]
        t6 <- valid[j,"tumor_6_months"]
        t1 <- valid[j,"tumor_1_year"]
        
        # conditional block to handle (potential) imputation and classification
        if (!is.na(t1)) {
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_1_year=t1))
        } else if (!is.na(td)) {
            if (!is.na(t6)) {
                t1 <- predict(model3, newdata=data.frame(tumor_6_months=t6, tumor_diagnosis=td))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_1_year=t1))
            } else {
                t1 <- predict(model1, newdata=data.frame(tumor_diagnosis=td))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_1_year=t1))
            }
        } else if (!is.na(t6)) {
            t1 <- predict(model2, newdata=data.frame(tumor_6_months=t6))
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_1_year=t1))
        } else {
            # guess mean value if all tumor_x are NA
            t1 <- meanTumor
            surv[j,"logit"] <- predict(model4, newdata=data.frame(tumor_1_year=t1))
        }
    }

    # evaluate and record accuracy of model
    for (j in 1:length_valid) {
        # logit -> prob:
        # logit = L
        # corresponding probability = exp(L)/(1+exp(L))
        L = surv[j,"logit"]
        if (L > 0.0) {
            surv[j,"pred"] <- 1
        } else {
            surv[j,"pred"] <- 0
        }
        if (surv[j,"pred"] == surv[j,"true"]) {
            count <- count + 1
        }
    }
    a2ii[i] <- count / length_valid
}

# final accuracy
Accuracy2ii <- sum(a2ii)/10


### (2iii) Choose best available tumor_x #############################################################################

# cross-validation accuracies
a2iii <- c(0,0,0,0,0,0,0,0,0,0)

# train and test 10 models
for (i in 1:10) {
    
    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"logit"] <- NA # logit values returned from classification of validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    coefs <- data.frame(matrix(ncol=4, nrow=3)) # placeholder for regression model coefficients
    count <- 0 # count of correctly classified instances
    
    # train regression models
    temp_data <- train[c("survival_7_years", "tumor_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model1 = glm(formula = survival_7_years ~ tumor_diagnosis, data=temp_data)

    temp_data <- train[c("survival_7_years", "tumor_6_months")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model2 = glm(formula = survival_7_years ~ tumor_6_months, data=temp_data)
    
    temp_data <- train[c("survival_7_years", "tumor_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model3 = glm(formula = survival_7_years ~ tumor_1_year, data=temp_data)
    
    # test classification model
    # tumor_1_year > tumor_diagnosis > tumor_6_months
    for (j in 1:length_valid) {
        
        # fields
        td <- valid[j,"tumor_diagnosis"]
        t6 <- valid[j,"tumor_6_months"]
        t1 <- valid[j,"tumor_1_year"]
        
        # conditional block to handle classification
        if (!is.na(t1)) {
            surv[j,"logit"] <- predict(model3, newdata=data.frame(tumor_1_year=t1),)
        } else if (!is.na(td)) {
            surv[j,"logit"] <- predict(model1, newdata=data.frame(tumor_diagnosis=td))
        } else if (!is.na(t6)) {
            surv[j,"logit"] <- predict(model2, newdata=data.frame(tumor_6_months=t6))
        } else {
            # guess mean value if all tumor_x are NA
            t1 <- mean(temp_data$tumor_1_year)
            surv[j,"logit"] <- predict(model3, newdata=data.frame(tumor_1_year=t1))
        }
    }

    # evaluate and record accuracy of model
    for (j in 1:length_valid) {

        L = surv[j,"logit"]
        if (L > 0.0) {
            surv[j,"pred"] <- 1
        } else {
            surv[j,"pred"] <- 0
        }
        if (surv[j,"pred"] == surv[j,"true"]) {
            count <- count + 1
        }
    }
    a2iii[i] <- count / length_valid
}

# final accuracy
Accuracy2iii <- sum(a2iii)/10

### (3) Logit Regression Model of Single Input (psa_x) ###############################################################

### (3i) Choose psa_diagnosis, impute if missing #####################################################################

# cross-validation accuracies
a3i <- c(0,0,0,0,0,0,0,0,0,0)

# train and test 10 models
for (i in 1:10) {
    
    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"logit"] <- NA # logit values returned from classification of validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    coefs <- data.frame(matrix(ncol=4, nrow=3)) # placeholder for regression model coefficients
    count <- 0 # count of correctly classified instances
    
    # train regression models
    temp_data <- train[c("psa_diagnosis", "psa_6_months")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model1 = lm(formula = psa_diagnosis ~ psa_6_months, data=temp_data)
    coefs[1:length(coef(model1)),1] <- coef(model1)
    
    temp_data <- train[c("psa_diagnosis", "psa_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model2 = lm(formula = psa_diagnosis ~ psa_1_year, data=temp_data)
    coefs[1:length(coef(model2)),2] <- coef(model2)
    
    temp_data <- train[c("psa_diagnosis", "psa_6_months", "psa_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model3 = lm(formula = psa_diagnosis ~ psa_6_months + psa_1_year, data=temp_data)
    coefs[1:length(coef(model3)),3] <- coef(model3)
    
    temp_data <- train[c("survival_7_years", "psa_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model4 = glm(formula = survival_7_years ~ psa_diagnosis, data=temp_data, family="binomial")
    coefs[1:length(coef(model4)),4] <- coef(model4)
    meanPSA <- mean(temp_data$psa_diagnosis)
    
    # test classification model
    for (j in 1:length_valid) {
        
        # fields
        pd <- valid[j,"psa_diagnosis"]
        p6 <- valid[j,"psa_6_months"]
        p1 <- valid[j,"psa_1_year"]
        
        # conditional block to handle (potential) imputation and classification
        if (!is.na(pd)) {
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_diagnosis=pd),)
        } else if (!is.na(p6)) {
            if (!is.na(p1)) {
                pd <- predict(model3, newdata=data.frame(psa_6_months=p6, psa_1_year=p1))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_diagnosis=pd))
            } else {
                pd <- predict(model1, newdata=data.frame(psa_6_months=p6))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_diagnosis=pd))
            }
        } else if (!is.na(p1)) {
            pd <- predict(model2, newdata=data.frame(psa_1_year=p1))
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_diagnosis=pd))
        } else {
            # guess mean value if all psa_x are NA
            pd <- meanPSA
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_diagnosis=pd))
        }
    }

    # evaluate and record accuracy of model
    for (j in 1:length_valid) {
        L = surv[j,"logit"]
        if (L > 0.0) {
            surv[j,"pred"] <- 1
        } else {
            surv[j,"pred"] <- 0
        }
        if (surv[j,"pred"] == surv[j,"true"]) {
            count <- count + 1
        }
    }
    a3i[i] <- count / length_valid
}

# final accuracy
Accuracy3i <- sum(a3i)/10


### (3ii) Choose psa_1_year, impute if missing #######################################################################

# cross-validation accuracies
a3ii <- c(0,0,0,0,0,0,0,0,0,0)

# train and test 10 models
for (i in 1:10) {
    
    # create training and validation sets
    valid <- folded_data[[i]]
    train <- folded_data
    train[[i]] <- NULL
    train <- do.call(Map, c(c, train))
    train <- as.data.frame(train)
    
    # fields
    length_train <- length(train[,31]) # number of training instances
    length_valid <- length(valid[,31]) # number of validation instances
    surv <- data.frame(true=valid[,31]) # true survival values from validation set
    surv[,"logit"] <- NA # logit values returned from classification of validation set
    surv[,"pred"] <- NA # predicted survival values for validation set
    coefs <- data.frame(matrix(ncol=4, nrow=3)) # placeholder for regression model coefficients
    count <- 0 # count of correctly classified instances
    
    # train regression models
    temp_data <- train[c("psa_1_year", "psa_6_months")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model1 = lm(formula = psa_1_year ~ psa_6_months, data=temp_data)
    coefs[1:length(coef(model1)),1] <- coef(model1)
    
    temp_data <- train[c("psa_1_year", "psa_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model2 = lm(formula = psa_1_year ~ psa_diagnosis, data=temp_data)
    coefs[1:length(coef(model2)),2] <- coef(model2)
    
    temp_data <- train[c("psa_1_year", "psa_6_months", "psa_diagnosis")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model3 = lm(formula = psa_1_year ~ psa_6_months + psa_diagnosis, data=temp_data)
    coefs[1:length(coef(model3)),3] <- coef(model3)
    
    temp_data <- train[c("survival_7_years", "psa_1_year")]
    temp_data <- temp_data[complete.cases(temp_data),]  
    model4 = glm(formula = survival_7_years ~ psa_1_year, data=temp_data, family="binomial")
    coefs[1:length(coef(model4)),4] <- coef(model4)
    meanPSA <- mean(temp_data$psa_1_year)
    
    # test classification model
    for (j in 1:length_valid) {
        
        # fields
        pd <- valid[j,"psa_diagnosis"]
        p6 <- valid[j,"psa_6_months"]
        p1 <- valid[j,"psa_1_year"]
        
        # conditional block to handle (potential) imputation and classification
        if (!is.na(p1)) {
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_1_year=p1),)
        } else if (!is.na(p6)) {
            if (!is.na(pd)) {
                p1 <- predict(model3, newdata=data.frame(psa_6_months=p6, psa_diagnosis=pd))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_1_year=p1))
            } else {
                p1 <- predict(model1, newdata=data.frame(psa_6_months=p6))
                surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_1_year=p1))
            }
        } else if (!is.na(pd)) {
            p1 <- predict(model2, newdata=data.frame(psa_diagnosis=pd))
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_1_year=p1))
        } else {
            # guess mean value if all psa_x are NA
            p1 <- meanPSA
            surv[j,"logit"] <- predict(model4, newdata=data.frame(psa_1_year=p1))
        }
    }

    # evaluate and record accuracy of model
    for (j in 1:length_valid) {
        L = surv[j,"logit"]
        if (L > 0.0) {
            surv[j,"pred"] <- 1
        } else {
            surv[j,"pred"] <- 0
        }
        if (surv[j,"pred"] == surv[j,"true"]) {
            count <- count + 1
        }
    }
    a3ii[i] <- count / length_valid
}

# final accuracy
Accuracy3ii <- sum(a3ii)/10



######################################################################################################################
### DETERMINISTIC IMPUTATIONS
### GOALS:  (1) NA analysis - high accuracy in predicting death
###         (2) 0 tumor/psa - 100% life expectency
###
######################################################################################################################


### (1) NA analysis ##################################################################################################

score_data <- read.csv(file="training_data.csv", header=TRUE, sep=",")
tps1_data <- score_data[c("tumor_1_year", "psa_1_year", "survival_1_year")]
tps2_data <- score_data[c("tumor_6_months", "tumor_1_year", "psa_6_months", "psa_1_year", "survival_1_year")]

# 2 NAs present, in both 1 year slots for tumor and psa
for (i in 1:dim(tps1_data)[1]) {
    temp_psa <- tps1_data[i,2]
    temp_tum <- tps1_data[i,1]
    temp_sur <- tps1_data[i,3]
    # if PSA and Tumor are both NA, set survival_1_year=0
    if (is.na(temp_psa) && is.na(temp_tum)) {
        score_data[i,30] <- 0
    }
}

# 4 NA's present in all 6 month/1 year positions for tumor and psa
for (i in 1:dim(tps2_data)[1]) {
    temp_psa1 <- tps2_data[i,3]
    temp_tum1 <- tps2_data[i,1]
    temp_sur <- tps2_data[i,5]
    temp_psa2 <- tps2_data[i,4]
    temp_tum2 <- tps2_data[i,2]
    # if PSA and Tumor are all NA, set survival_1_year=0
    if (is.na(temp_psa1) && is.na(temp_tum1)) {
        if (is.na(temp_psa2) && is.na(temp_tum2)) {
            score_data[i,30] <- 0
        }
    }
}


#-----------------------------------score data
    
    # create training and validation sets
i <-1
valid <- folded_data[[i]]
train <- folded_data
train[[i]] <- NULL
train <- do.call(Map, c(c, train))
train <- as.data.frame(train)

# adjust categorical data types for Bayes
valid[,"survival_7_years"] <- as.factor(valid$survival_7_years)
valid[,"gleason_score"] <- ordered(valid$gleason_score)
valid[,"age_group"] <- ordered(valid$age_group)

#train[,32] <- NULL
train[,"gleason_score"] <- ordered(train$gleason_score)
train[,"t_score"] <- ordered(train$t_score)
#train[,"count_therpy"] <- ordered(train$count_therpy)
levels(train[,"t_score"]) <- c("T1a", "T1b", "T1c", "T2a", "T2b", "T2c", "T3a", "T3b", "T3c", "T4")
train[,"n_score"] <- as.factor(train$n_score)
levels(train[,"n_score"]) <- c("N0", "N1", "NX")
train[,"side"] <- as.factor(train$side)
levels(train[,"side"]) <- c("both", "left", "right")
train[,"m_score"] <- ordered(train$m_score)
levels(train[,"m_score"]) <- c("M0", "M1a", "M1b", "M1c")
train[,"stage"] <- ordered(train$stage)
levels(train[,"stage"]) <- c("I", "IIA", "IIB", "III", "IV")
train[,"age_group"] <- ordered(train$age_group)
train[,"survival_7_years"] <- as.factor(train$survival_7_years)

# fields
length_train <- length(train[,31]) # number of training instances
length_valid <- length(valid[,31]) # number of validation instances
surv <- data.frame(true=valid[,31]) # true survival values from validation set
surv[,"pred"] <- NA # predicted survival values for validation set
count <- 0 # count of correctly classified instances

# train and test naive Bayes classifier
bayesModel <- naiveBayes(survival_7_years ~ gleason_score + t_score + m_score + n_score + stage, data=train)
score_data[,"survival_7_years"] <- predict(bayesModel, score_data[1:11531,])

# evaluate and record accuracy of model
for (j in 1:length_valid) {
    if (!is.na(surv[j,"pred"])) {
        if (!is.na(surv[j,"true"])) {
            if (surv[j,"pred"] == surv[j,"true"]) {
                count <- count +1
            }
        }
    }  
}

