
### LOAD PACKAGES
library(Amelia)
library(arules)
library(rpart)
library(scatterplot3d)
library(e1071)

### READ DATA, PREPROCESS, AND SAMPLE
setwd("/Users/cfrehr/Job_Search/Interviews/Enova/Interview_Assignment/participant_files/")
patient_data <- read.csv(file="training_data.csv", header=TRUE, sep=",")
score_data <- read.csv(file="CodyFrehr_score.csv", header=TRUE, sep=",")

# Preprocess patient_data data types
patient_data$t_score <- ordered(patient_data$t_score, levels(patient_data$t_score)[c(1:10)])
patient_data$m_score <- ordered(patient_data$m_score, levels(patient_data$m_score)[c(1:4)])
patient_data$stage <- ordered(patient_data$stage, levels(patient_data$stage)[c(1:5)])
patient_data$race <- as.factor(patient_data$race)

# Preprocess score_data data types
score_data$t_score <- ordered(score_data$t_score, levels(score_data$t_score)[c(1:10)])
score_data$m_score <- ordered(score_data$m_score, levels(score_data$m_score)[c(1:4)])
score_data$stage <- ordered(score_data$stage, levels(score_data$stage)[c(1:5)])
score_data$race <- as.factor(score_data$race)


######################################################################################################################
### DATA VISUALIZATION: An investigation into the completeness of data sets
### GOALS:  (1) Identify the missingness levels of individual variables
###         (2) Compare completeness of both data sets to sort out irreparably incomplete variables from those
###             worthy of further investigation, imputation, and integration into predictive models
######################################################################################################################

# Missingness maps:
#missmap(patient_data) # training data
#missmap(test_data) # testing data

# Imputation possibilities:
# height <--> weight (regression analysis)
# height & weight --> race (cluster segmentation)
# psa_x <--> tumor_x (regression analysis)
# psa_x1 <--> psa_x2 OR tumor_x1 <--> tumor_x2 (regression analysis)
# psa_1_year OR tumor_1_year --> survival_1_year (logic rules)
# TNM scores & stage --> gleason_score (baeyesian classifier)


######################################################################################################################
### DATA ANALYSIS: height & weight 
### GOALS:  (1) Explore how well height correlates with weight
###         (2) Explore how well height/weight correlate with survival_1_year and survival_7_years
###         (3) Impute missing values?
######################################################################################################################


### (1) Explore how well height correlates with weight

# Create relevant dataframe, excluding entities with NA values
hw_data <- patient_data[c("height", "weight", "survival_1_year", "survival_7_years")]
hw_data <- hw_data[complete.cases(hw_data),]

# Verify normality and create scatterplot
minh <- min(hw_data$height)
maxh <- max(hw_data$height)
minw <- min(hw_data$weight)
maxw <- max(hw_data$weight)
hist(hw_data$height, col="gray", main="Histogram of Patient Heights", xlab = "Height (in)", breaks=seq(minh,maxh,l=10))
hist(hw_data$weight, col="gray", main="Histogram of Patient Weights", xlab = "Weight (lbs)", breaks=seq(minw,maxw,l=20))
plot(hw_data$height, hw_data$weight, type = "p" , main="Scatterplot of Height vs Weight", xlab = "Height (in)", ylab="Weight (lbs)")

# Simple regression model and fitted line
model1 = lm(formula = weight ~ height, data = hw_data)
summary(model1)
coefs1 <- coefficients(model1)
alpha1 <- coefs1[1]
beta1 <- coefs1[2]
abline(alpha1, beta1)

# Takeaways: 
# - height and weight are weakly related in prostate cancer patients


### (2) Explore how well height/weight correlate with survival_1_year & survival_7_years

model2 = glm(formula = survival_1_year ~ height + weight, data = hw_data, family="binomial")
summary(model2)
coefs2 <- coefficients(model2)

model2.5 = glm(formula = survival_7_years ~ height + weight, data = hw_data, family="binomial")
summary(model2.5)
coefs2.5 <- coefficients(model2.5)

# Takeaways:
# - height has no ability to predict survival rates, so it is not worth working to replace missing values
# - weight has a significant effect on survival_1_year, but not survial_7_years; may be worth imputing missing values


### (3) Impute missing values with regression model

# Predict missing height and weight values
mean_height <- round(mean(hw_data$height), digits = 0)
mean_weight <- round(mean(hw_data$weight), digits = 0)
temp_height <- 0
temp_weight <- 0
new_height <- 0
new_weight <- 0

for (i in 1:dim(patient_data)[1]) {
  temp_height <- patient_data[i,9]
  # if height is NA
  if (is.na(temp_height)) {
    temp_weight <- patient_data[i,10]
    # if weight is not NA, fit regression value
    if (!is.na(temp_weight)) {
      new_height <- round((temp_weight-alpha1)/beta1, digits = 0)
      patient_data[i,9] <- new_height
    # else if weight is NA, fill both values with mean vals
    } else {
      patient_data[i,9] <- mean_height
      patient_data[i,10] <- mean_weight
    }
  }
  temp_weight <- patient_data[i,10]
  # else if weight is NA
  if (is.na(temp_weight)) {
    # if height is not NA, fit regression value
    if (!is.na(temp_height)) {
      new_weight <- round(alpha1 + beta1*temp_height, digits = 0)
      patient_data[i,10] <- new_weight
    }
  }
}


######################################################################################################################
### DATA ANALYSIS: tumor_diagnosis & psa_diagnosis
### GOALS:  (1) Explore how well tumor_diagnosis correlates with psa_diagnosis
###         (2) Explore how well change in tumor size and psa levels over time time can predict missing values
###         (3) Explore interactions between tumor size and cancer treatment methods
###         (4) Explore how change in tumor size and psa levels over time affect a patients survivability
###         (5) Impute missing values?
######################################################################################################################


### (1) Explore how well tumor_diagnosis correlates with psa_diagnosis

# Create relevant dataframe, excluding entities with NA values
tp_data <- patient_data[c("tumor_diagnosis", "psa_diagnosis")]
tp_data <- tp_data[complete.cases(tp_data),]

# Verify normality and create scatterplot
mint <- min(tp_data$tumor_diagnosis)
maxt <- max(tp_data$tumor_diagnosis)
minp <- min(tp_data$psa_diagnosis)
maxp <- max(tp_data$psa_diagnosis)
hist(tp_data$tumor_diagnosis, col="gray", main="Histogram of Patient Tumor Sizes at Diagnosis", xlab = "Tumor Size (mm)", breaks=seq(mint,maxt,l=20))
hist(tp_data$psa_diagnosis, col="gray", main="Histogram of Patient PSA Levels at Diagnosis", xlab = "PSA Level (ng/mL)", breaks=seq(minp,maxp,l=20))
plot(tp_data$tumor_diagnosis, tp_data$psa_diagnosis, type = "p" , main="Scatterplot of Tumor Size vs PSA Level", xlab = "Tumor Size (mm)", ylab="PSA Level (ng/mL)")

# Simple regression model and fitted line
model3 = lm(formula = psa_diagnosis ~ tumor_diagnosis, data = tp_data)
summary(model3)
coefs3 <- coefficients(model3)
alpha3 <- coefs3[1]
beta3 <- coefs3[2]
abline(alpha3, beta3)


# Takeaways:
# - Little to no correlation between tumor size and psa level


### (2) Explore how well change in tumor size and psa levels over time time can predict missing values

# Plots of change in tumor size over time
t_data <- patient_data[c("tumor_diagnosis", "tumor_6_months", "tumor_1_year", "survival_1_year", "survival_7_years")]
t_data <- t_data[complete.cases(t_data),]

plot(t_data$tumor_diagnosis, t_data$tumor_1_year, type="p", main="Scatterplot of Tumor Sizes Over 1 Year", xlab="Tumor Size at Diagnosis (mm)", ylab="Tumor Size at 1 Year (mm)")
plot(t_data$tumor_diagnosis, t_data$tumor_6_Months, type="p", main="Scatterplot of Tumor Sizes Over First 6 Months", xlab="Tumor Size at Diagnosis (mm)", ylab="Tumor Size at 6 Months (mm)")
plot(t_data$tumor_6_months, t_data$tumor_1_year, type="p", main="Scatterplot of Tumor Sizes Over Second 6 Months", xlab="Tumor Size at 6 Months (mm)", ylab="Tumor Size at 1 Year (mm)")
scatterplot3d(t_data$tumor_diagnosis, t_data$tumor_1_year, t_data$tumor_6_months, main="3D Scatterplot of Tumor Size Over 1 Year", xlab="Tumor Size at Diagnosis (mm)", ylab="Tumor Size at 1 Year (mm)", zlab="Tumor Size at 6 Months (mm)")
scatterplot3d(t_data$tumor_diagnosis, t_data$tumor_1_year, t_data$tumor_6_months, angle=180, main="3D Scatterplot of Tumor Size Over 1 Year", xlab="Tumor Size at Diagnosis (mm)", ylab="Tumor Size at 1 Year (mm)", zlab="Tumor Size at 6 Months (mm)")

# Pairwise regression models
model4 = lm(formula = tumor_1_year ~ tumor_diagnosis, data=t_data)
summary(model4)
coefs4 <- coefficients(model4)

model5 = lm(formula = tumor_1_year ~ tumor_6_months, data = t_data)
summary(model5)
coefs5 <- coefficients(model5)

# Complete regression model
model6 = lm(formula = tumor_1_year ~ tumor_diagnosis + tumor_6_months, data = t_data)
summary(model6)
coefs6 <- coefficients(model6)

# Takeaways:
# - With only 1 tumor size data point, a relatively strong prediction can be made about a tumors size at 1 year
#   With 2 data points, very strong conclusions can be reached about a tumors size at 1 year
# - Deterministic cutoff-point of 27mm in tumor_diagnosis for 100% accuracy in prediction of 0mm tumor_1_year (6th percentile)
# - Deterministic cutoff-point of 7mm in tumor_6_months for 100% accuracy in prdiction of 0mm tumor_1_year (13th percentile)


# Plots of change in psa level over time
p_data <- patient_data[c("psa_diagnosis", "psa_6_months", "psa_1_year", "survival_1_year", "survival_7_years")]
p_data <- p_data[complete.cases(p_data),]
plot(p_data$psa_diagnosis, p_data$psa_1_year, type="p", main="Scatterplot of PSA Levels Over 1 Year", xlab="PSA Level at Diagnosis (ng/mL)", ylab="PSA Level at 1 Year (ng/mL)")
plot(p_data$psa_diagnosis, p_data$psa_6_Months, type="p", main="Scatterplot of PSA Levels Over First 6 Months", xlab="PSA Level at Diagnosis (ng/mL)", ylab="PSA Level at 6 Months (ng/mL)")
plot(p_data$psa_6_months, p_data$psa_1_year, type="p", main="Scatterplot of PSA Levels Over Second 6 Months", xlab="PSA Level at 6 Months (ng/mL)", ylab="PSA Level at 1 Year (ng/mL)")
scatterplot3d(p_data$psa_diagnosis, p_data$psa_1_year, p_data$psa_6_months, main="3D Scatterplot of PSA Levels Over 1 Year", xlab="PSA Level at Diagnosis (ng/mL)", ylab="PSA Level at 1 Year (ng/mL)", zlab="PSA Level at 6 Months (ng/mL)")
scatterplot3d(p_data$psa_diagnosis, p_data$psa_1_year, p_data$psa_6_months, angle=168, main="3D Scatterplot of PSA Levels Over 1 Year", xlab="PSA Level at Diagnosis (ng/mL)", ylab="PSA Level at 1 Year (ng/mL)", zlab="PSA Level at 6 Months (ng/mL)")

# Pairwise regression models
model8 = lm(formula = psa_1_year ~ psa_diagnosis, data = p_data)
summary(model8)
coefs8 <- coefficients(model8)

model9 = lm(formula = psa_1_year ~ psa_6_months, data = p_data)
summary(model9)
coefs9 <- coefficients(model9)

# Complete regression model
model10 = lm(formula = psa_1_year ~ psa_diagnosis + psa_6_months, data = p_data)
summary(model10)
coefs10 <- coefficients(model10)

# Takeaways:
# - With only 1 psa data point, a relatively strong prediction can be made about a psa level at 1 year
#   With 2 data points, relatively strong conclusions can be reached about psa levels at 1 year
# - Deterministic cutoff-point of 4.3 ng/mL in psa_diagnosis for 100% accuracy in prediction of 0 ng/mL psa_1_year (6th percentile)
# - Deterministic cutoff-point of 1 ng/mL in psa_6_months for 100% accuracy in prdiction of 0 ng/mL psa_1_year (4th percentile)

### (3) Explore interactions between tumor size and cancer treatment methods


### (4) Explore how change in tumor size and psa levels over time affect a patients survivability

t2_data <- patient_data[c("tumor_diagnosis", "survival_1_year", "survival_7_years")]
t2_data <- t2_data[complete.cases(t2_data),]
t3_data <- patient_data[c("tumor_6_months", "survival_1_year", "survival_7_years")]
t3_data <- t3_data[complete.cases(t3_data),]
t4_data <- patient_data[c("tumor_1_year", "survival_1_year", "survival_7_years")]
t4_data <- t4_data[complete.cases(t4_data),]

p2_data <- patient_data[c("psa_diagnosis", "survival_1_year", "survival_7_years")]
p2_data <- p2_data[complete.cases(p2_data),]
p3_data <- patient_data[c("psa_6_months", "survival_1_year", "survival_7_years")]
p3_data <- p3_data[complete.cases(p3_data),]
p4_data <- patient_data[c("psa_1_year", "survival_1_year", "survival_7_years")]
p4_data <- p4_data[complete.cases(p4_data),]

# survival_1_year
model11 = glm(formula = survival_1_year ~ tumor_diagnosis, data = t2_data, family="binomial")
summary(model11)
coefs11 <- coefficients(model11)

model11.2 = glm(formula = survival_1_year ~ tumor_6_months, data = t3_data, family="binomial")
summary(model11.2)
coefs11.2 <- coefficients(model11.2)

model12 = glm(formula = survival_1_year ~ psa_diagnosis, data = p2_data, family="binomial")
summary(model12)
coefs12 <- coefficients(model12)

model12.2 = glm(formula = survival_1_year ~ psa_6_months, data = p3_data, family="binomial")
summary(model12.2)
coefs12.2 <- coefficients(model12.2)

# Takeaways:
# - best individual predictors: tumor_diagnosis and psa_6_months are the best individual predictors of survival_1_year

# survival_7_years
43 29 51 43 21 23 56 35 39 44
0  1  
model13 = glm(formula = survival_7_years ~ tumor_diagnosis, data = t2_data, family="binomial")
summary(model13)
coefs13 <- coefficients(model13)
L <- predict(model13, newdata=data.frame(tumor_diagnosis=111), type=c("link"))
P <- 1/(1+exp(-1*L))
P

model13.2 = glm(formula = survival_7_years ~ tumor_6_months, data = t3_data, family="binomial")
summary(model13.2)
coefs13.2 <- coefficients(model13.2)

model13.3 = glm(formula = survival_7_years ~ tumor_1_year, data = t4_data, family="binomial")
summary(model13.3)
coefs13.3 <- coefficients(model13.3)
L <- predict(model13.3, newdata=data.frame(tumor_diagnosis=111), type=c("link"))
P <- 1/(1+exp(-1*L))
P

model13.4 = glm(formula = survival_7_years ~ tumor_diagnosis + tumor_6_months + tumor_1_year, data = t_data, family="binomial")
summary(model13.4)
coefs13.4 <- coefficients(model13.4)

model13.5 = glm(formula = survival_7_years ~ tumor_diagnosis + tumor_6_months, data = t_data, family="binomial")
summary(model13.5)
coefs13.5 <- coefficients(model13.5)

model14 = glm(formula = survival_7_years ~ psa_diagnosis, data = p2_data, family="binomial")
summary(model14)
coefs14 <- coefficients(model14)

model14.2 = glm(formula = survival_7_years ~ psa_6_months, data = p3_data, family="binomial")
summary(model14.2)
coefs14.2 <- coefficients(model14.2)

model14.3 = glm(formula = survival_7_years ~ psa_1_year, data = p4_data, family="binomial")
summary(model14.3)
coefs14.3 <- coefficients(model14.3)

model14.4 = glm(formula = survival_7_years ~ psa_diagnosis + psa_6_months + psa_1_year, data = p_data, family="binomial")
summary(model14.4)
coefs14.4 <- coefficients(model14.4)


# Takeaways:
# - although there is strong multicollinearity present within these time series, the error is so tightly controlled in
#   this model that I feel comfortable with its accuracy despite the dependencies
# - tumor sizes at all times are very strong predictors of survival. tumor_1_year offers strong predictability with
#   a standard error half the size of tumor_diagnosis and tumor_6_months
# - psa levels at all times are very strong predictors of survival. psa_1_year offers the highest level of significance
#   and lowest standard error between the three
#
# Best individual predictors: tumor_1_year > tumor_diagnosis > tumor_6_months
#                             psa_1_year > psa_6_months > psa_diagnosis

### (6) Impute missing values

# Not worth while because in all cases, there is at least one observation with which to use as a predictor, and all 
# variables are very strong individual predictors. Additionally, there is the concern of multicollinearity present
# between the time-series variables.

# IF using a specific value (like tumor_diagnosis or psa_6_months) is worth having for their predictability of other
# non-class attributes, then imputation may be worth conducting.



######################################################################################################################
### DATA ANALYSIS: gleason_score & TNM staging
### GOALS:  (1) Explore associations between gleason_scores, TNM scores, and stage values
###         (2) Explore how well gleason_scores, TNM scores, and stage predict survivability
###         (3) Explore interactions between cancer stage and cancer treatment methods
###         (4) Impute missing values?
######################################################################################################################


### (1) Explore how strongly related gleason_score, TNM scores, and stage are related

# Create relevant dataframe, excluding entities with NA values
tnm_data <- patient_data[c("gleason_score", "t_score", "n_score", "m_score", "stage", "survival_1_year", "survival_7_years")]
tnm_data <- tnm_data[complete.cases(tnm_data),]

# Get sense of distributions
ming <- min(tnm_data$gleason_score)
maxg <- max(tnm_data$gleason_score)
hist(tnm_data$gleason_score, col="gray", main="Histogram of Gleason Scores", xlab = "Gleason Score")
barplot(table(tnm_data$t_score), main="Histogram of T Scores", ylab = "Frequencies")
barplot(table(tnm_data$m_score), main="Histogram of M Scores", ylab = "Frequencies")
barplot(table(tnm_data$n_score), main="Barplot of N Scores", ylab = "Frequencies")
barplot(table(tnm_data$stage), main="Histogram of Stages", ylab = "Frequencies")

# 


### (2) Explore how well gleason_scores, TNM scores, and stage predict survivability

model15 = glm(formula = survival_1_year ~ gleason_score, data = tnm_data, family="binomial")
summary(model15)
coefs15 <- coefficients(model15)

model16 = glm(formula = survival_7_years ~ gleason_score, data = tnm_data, family="binomial")
summary(model16)
coefs16 <- coefficients(model16)

model17 = glm(formula = survival_1_year ~ t_score, data = tnm_data, family="binomial")
summary(model17)
coefs17 <- coefficients(model17)

model18 = glm(formula = survival_7_years ~ t_score, data = tnm_data, family="binomial")
summary(model18)
coefs18 <- coefficients(model18)

model19 = glm(formula = survival_1_year ~ n_score, data = tnm_data, family="binomial")
summary(model19)
coefs19 <- coefficients(model19)

model20 = glm(formula = survival_7_years ~ n_score, data = tnm_data, family="binomial")
summary(model20)
coefs20 <- coefficients(model20)

model21 = glm(formula = survival_1_year ~ m_score, data = tnm_data, family="binomial")
summary(model21)
coefs21 <- coefficients(model21)

model22 = glm(formula = survival_7_years ~ m_score, data = tnm_data, family="binomial")
summary(model22)
coefs22 <- coefficients(model22)

model23 = glm(formula = survival_1_year ~ stage, data = tnm_data, family="binomial")
summary(model23)
coefs23 <- coefficients(model23)

model24 = glm(formula = survival_7_years ~ stage, data = tnm_data, family="binomial")
summary(model24)
coefs24 <- coefficients(model24)

tnm_data <- patient_data[c("gleason_score", "t_score", "n_score", "m_score", "stage", "survival_1_year", "survival_7_years")]
tnm_data <- tnm_data[complete.cases(tnm_data),]

tnm_data$gleason_score <- as.factor(tnm_data$gleason_score)
tnm_data$gleason_score <- ordered(tnm_data$gleason_score)
tnm_data$t_score <- as.integer(tnm_data$t_score)
tnm_data$m_score <- as.integer(tnm_data$m_score)
tnm_data$stage <- as.integer(tnm_data$stage)

model16 = glm(formula = survival_7_years ~ gleason_score, data = tnm_data, family="binomial")
summary(model16)
coefs16 <- coefficients(model16)

model18 = glm(formula = survival_7_years ~ t_score, data = tnm_data, family="binomial")
summary(model18)
coefs18 <- coefficients(model18)

model22 = glm(formula = survival_7_years ~ m_score, data = tnm_data, family="binomial")
summary(model22)
coefs22 <- coefficients(model22)

model24 = glm(formula = survival_7_years ~ stage, data = tnm_data, family="binomial")
summary(model24)
coefs24 <- coefficients(model24)


# Takeaways:
# - at all levels of values, m_score and stage have the highest significance on predicting survival_7_years.
# - gleason_score, t_score, and n_score have predictive ability, but lack complete significance across levels


######################################################################################################################
### DATA ANALYSIS: survival_1_year
### GOALS:  (1) Verify that non-NA values for tumor_1_year or psa_1_year values imply with full accuracy that 
###             survival_1_year is true.
###         (2) Fill in value of 1 for all patients who meet the logical criteria
######################################################################################################################


### (1) Verify that non-NA values for tumor_1_year or psa_1_year values imply with full accuracy that survival_1_year 
###     is true.

tps_data <- patient_data[c("tumor_1_year", "psa_1_year", "survival_1_year")]
tps_data <- tps_data[complete.cases(tps_data),]

# Look for association between non-NA tumor_1_year, non-NA psa_1_year, and survial_1_year values
temp_tumor <- 0
temp_psa <- 0
temp_survival <- 0
tumor_count <- 0
psa_count <- 0
survival_count <- 0

# check for matching number of cases between tumors and survival
for (i in 1:dim(tps_data)[1]) {
  temp_tumor <- tps_data[i,1]
  # if PSA is non-NA, increase psa_count
  if (!is.na(temp_tumor)) {
    tumor_count <- tumor_count + 1
    temp_survival <- tps_data[i,3]
    # if patient survived, increase survival_count
    if (!is.na(temp_survival)) {
      if (temp_survival == 1) {
        survival_count <- survival_count + 1
      }
    }
  }
}

# if (tumor_count == survival_count) {
#   print("Tumor-Survival Match")
# } else {
#   print("Tumor-Survival Error")
# }

# check for matching number of cases between PSA and survival
survival_count <- 0
for (i in 1:dim(tps_data)[1]) {
  temp_psa <- tps_data[i,2]
  # if PSA is non-NA, increase psa_count
  if (!is.na(temp_psa)) {
    psa_count <- psa_count + 1
    temp_survival <- tps_data[i,3]
    # if patient survived, increase survival_count
    if (!is.na(temp_survival)) {
      if (temp_survival == 1) {
        survival_count <- survival_count + 1
      }
    }
  }
}

# if (psa_count == survival_count) {
#   print("PSA-Survival Match")
# } else {
#   print("PSA-Survival Error")
# }

### (2) Fill in value of 1 for all patients who meet the logical criteria
survival_count <- 0
for (i in 1:dim(tps_data)[1]) {
  temp_tumor <- tps_data[i,1]
  temp_psa <- tps_data[i,2]
  temp_survival <- tps_data[i,3]
  # if tumor or PSA is non-NA, and surival is NA, replace with 1
  if (!is.na(temp_tumor) | !is.na(temp_psa)) {
    if (is.na(temp_survival)) {
      tps_data[i,3] <- 1
      survival_count <- survival_count + 1
    }
  }
}



######################################################################################################################
### DETERMINISTIC IMPUTATIONS
### GOALS:  (1) NA analysis

######################################################################################################################

patient_data <- read.csv(file="training_data.csv", header=TRUE, sep=",")
tps1_data <- patient_data[c("tumor_1_year", "psa_1_year", "survival_1_year")]
tps2_data <- patient_data[c("tumor_6_months", "tumor_1_year", "psa_6_months", "psa_1_year", "survival_1_year")]

na1_count <- 0
match1_count <- 0
for (i in 1:dim(tps1_data)[1]) {
    temp_psa <- tps1_data[i,2]
    temp_tum <- tps1_data[i,1]
    temp_sur <- tps1_data[i,3]
    # if PSA and Tumor are NA, increase na1_count
    if (is.na(temp_psa) && is.na(temp_tum)) {
        na1_count <- na1_count + 1
        # if survival_1_year = 0, also increase match1_count
        if (temp_sur == 0) {
            match1_count <- match1_count + 1
        }
    }
}

accuracyD1 <- na1_count/match1_count

na2_count <- 0
match2_count <- 0
for (i in 1:dim(tps2_data)[1]) {
    temp_psa1 <- tps2_data[i,3]
    temp_tum1 <- tps2_data[i,1]
    temp_sur <- tps2_data[i,5]
    temp_psa2 <- tps2_data[i,4]
    temp_tum2 <- tps2_data[i,2]
    # if PSA and Tumor are all NA, increase na2_count
    if (is.na(temp_psa1) && is.na(temp_tum1)) {
        if (is.na(temp_psa2) && is.na(temp_tum2)) {
            na2_count <- na2_count + 1
            # if survival_1_year = 0, also increase match2_count
            if (temp_sur == 0) {
                match2_count <- match2_count + 1
            }
        }
    }
}

accuracyD2 <- na2_count/match2_count

######################################################################################################################
### PREDICTIVE MODELS
### GOALS:  (1) Naive Bayesian Classifier (gleason_score, t_score, n_score, m_score, stage)
###               i)   Information Gain
###               ii)  Genie Index
###         (3) Logistic Regression Model of Single Input (tumor_x)
###               i)   Choose tumor_diagnosis, impute if missing
###               ii)  Choose tumor_1_year, impute if missing
###               iii) Choose best available tumor_x
###         (4) Logistic Regression Model of Single Input (psa_x)
###               i)   Choose psa_diagnosis, impute if missing
###               ii)  Choose psa_1_year, impute if missing
###               iii) Choose best available psa_x
######################################################################################################################

### Survivability Proportions/Heuristics (stage, age, side, age_group, count_thrpy, change_tumor, change_psa)
patient_data <- read.csv(file="training_data.csv", header=TRUE, sep=",")

prop_data <- patient_data[c("stage", "age", "side", "rd_thrpy", "h_thrpy", "chm_thrpy", "cry_thrpy", "brch_thrpy", "rad_rem", "survival_1_year", "survival_7_years", "t_score", "n_score", "m_score")]
prop_data <- prop_data[complete.cases(prop_data),]

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


for (i in 1:dim(prop_data)[1]) {
    
    # store count_thrpy
    prop_data[i,13] <- prop_data[i,4] + prop_data[i,5] + prop_data[i,6] + prop_data[i,7] + prop_data[i,8] + prop_data[i,9]
    
    # bin ages in age_group
    temp_age <- prop_data[i,2]
    if (temp_age < 50) {
        prop_data[i,12] <- 1
    } else if (temp_age < 60) {
        prop_data[i,12] <- 2
    } else if (temp_age < 70) {
        prop_data[i,12] <- 3
    } else if (temp_age < 80) {
        prop_data[i,12] <- 4
    } else if (temp_age < 90) {
        prop_data[i,12] <- 5
    } else if (temp_age < 100) {
        prop_data[i,12] <- 6
    } else {
        prop_data[i,12] <- 7
    }
}

# calculate prob of surving 1 year and 7 years based on these stage, age, side, and count_thrpy 

prop_data[,"age_group"] <- as.factor(prop_data$age_group)
prop_data[,"count_thrpy"] <- as.factor(prop_data$count_thrpy)
prop_data[,"survival_1_year"] <- as.numeric(prop_data$survival_1_year)

#t_score
temp <- prop_data[c("survival_1_year", "t_score")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("survival_7_years", "t_score")]
(tab <- prop.table(table(temp), 1))

#n_score
temp <- prop_data[c("n_score", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("n_score", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

# m_score
temp <- prop_data[c("m_score", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("m_score", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

# stage
temp <- prop_data[c("stage", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("stage", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

# age_group
temp <- prop_data[c("age_group", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("age_group", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

# side
temp <- prop_data[c("side", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("side", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

# count_thrpy
temp <- prop_data[c("count_thrpy", "survival_1_year")]
(tab <- prop.table(table(temp), 1))
temp <- prop_data[c("count_thrpy", "survival_7_years")]
(tab <- prop.table(table(temp), 1))

### Tumor_Diagnosis and PSA_Diagnosis upper and lower bounds



## SOLUTION POSSIBILITIES:
# -association rule learning (set covering optimization)
# -


## EXPLORING SIGNIFICANT THERAPY ASSOCIATIONS (FREQUENT PATTERN ANALYSIS)

therapy_data <- patient_data[c("rd_thrpy", "h_thrpy", "chm_thrpy", "cry_thrpy", "brch_thrpy", "rad_rem", "multi_thrpy", "survival_1_year", "survival_7_years")]
prop_data <- prop_data[complete.cases(prop_data),]


# convert attributes to factors
for (i in 1:9) {
  therapy_data[,i] <- as.factor(therapy_data[,i])
}
# convert dataframe to transactions
therapy_trans <- as(therapy_data, "transactions")
#inspect(head(therapy_trans))

# Find Frequent Items Sets Using Apriori Algo
capture.output(cat("\nFREQUENT ITEMSET GENERATION:  minSup = 1/100\n"), file = "output.txt")
capture.output(cat("--------------------------------------------\n"), file = "output.txt", append=TRUE)
itemSets <- apriori(therapy_trans, parameter = list(minlen=9, maxlen=9, supp=1/10, target = "frequent itemsets"))
capture.output(cat("\nFREQUENT ITEMSETS\n"), file = "output.txt", append = TRUE)
capture.output(cat("-----------------\n\n"), file = "output.txt", append=TRUE)
capture.output(inspect(sort(itemSets, by="support")), file = "output.txt", append=TRUE)
# Find max frequent itemsets (no frequent superset)
maxSets <- itemSets[is.maximal(itemSets), by="support"]
capture.output(cat("\nMAX FREQUENT ITEMSETS\n"), file = "output.txt", append = TRUE)
capture.output(cat("---------------------\n\n"), file = "output.txt", append=TRUE)
capture.output(inspect(sort(maxSets, by="support")), file = "output.txt", append=TRUE)
# Find closed frequent itemsets (no frequent superset of equal support)
closedSets <- itemSets[is.closed(itemSets), by="support"]
capture.output(cat("\nCLOSED FREQUENT ITEMSETS\n"), file = "output.txt", append=TRUE)
capture.output(cat("------------------------\n\n"), file = "output.txt", append=TRUE)
capture.output(inspect(sort(closedSets, by="support")), file = "output.txt", append=TRUE)

## Find Frequent Association Rules Using Apriori
rules <- apriori(therapy_trans, parameter = list(minlen=7, maxlen =7, supp=1/1000, conf=60/100, ext = F), appearance = list(rhs = c('survival_7_years=0',"survival_7_years=1"), default = 'lhs'))
capture.output(cat("\n\nFREQUENT ASSOCIATION RULE GENERATION:  minSup = 1/1000, minConf = 50/100\n"), file = "output.txt", append = TRUE)
capture.output(cat("-----------------------------------------------------------------------\n\n"), file = "output.txt", append=TRUE)
capture.output(inspect(sort(rules, by="lift")), file = "output.txt", append=TRUE)
