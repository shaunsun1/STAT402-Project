# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2
library(caret)  # used for knn and cross-validation
library(car)  # calculates VIF to check for multicollinearity
library(ResourceSelection)  # Hosmer-Lemeshow goodness-of-fit test

# Bring the chms_2018 dataset into R, and remove the weights as they are not used in this analysis
study_data = read.csv("chms_2018.csv")
study_data = select(study_data, -starts_with("BS"), -starts_with("WGT"))

# Make SMK_12 a factor
study_data$SMK_12 = as.factor(study_data$SMK_12)

# Choose seeds so that results from imputing censored data is fixed.
# These numbers were drawn from a Uniform distribution with boundaries (1, 10^9). This will not be truly random, but will suffice for our purposes
seed1 = 336741157
seed2 = 874515716

#----------------------------------------------------------------------------------------------------------------
# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to 
# be measured (ie. they are below the LOD). First understand how common this is
BCD_LOD = sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  
BCD_LOD  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)
BHG_LOD  # 999.5 appears 599 times in LAB_BHG

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is called jittering
set.seed(seed1)
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
set.seed(seed2)
study_data$LAB_BHG[study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)] = runif(BHG_LOD, 0, 2.1)

# Confirm that all values below the LOD have been replaced.
sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BCD
sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BHG


#----------------------------------------------------------------------------------------------------------------
# Use K-Nearest-Neighbors (KNN) to estimate missing values
# Check how many missing values there are in this data set
sum(is.na(study_data))  # 175 missing values

# Return variables that contain missing values, as well as their respective number of missing values
n_missing = colSums(is.na(study_data)) # number of missing values in each column of data set
matrix(c(names(study_data)[n_missing > 0], n_missing[n_missing > 0]), nrow = 2, byrow = TRUE)

# Are there observations that contain multiple missing values?
multiple_missing = rowSums(is.na(study_data))  # number of missing values in each row of data set
study_data[multiple_missing >= 2, seq(2, 8)] # return observations that have multiple missing values

# To simplify the following estimation, will estimate the LAB_BHG values for these two observations with the 
# sample mean of LAB_BHG. This insures that every observation has no more than 1 missing value
study_data$LAB_BHG[multiple_missing >= 2] = mean(study_data$LAB_BHG, na.rm = TRUE) 

multiple_missing = rowSums(is.na(study_data)) 
sum(multiple_missing > 1) # No observations have more than one missing value now

# Use cross-validation to find the best value of K to use in KNN 
ctrl <- trainControl(method="cv", number = 5) 

# Using HWMDBMI as the response variable, find the best value of K using as the training set only the observations 
# that have no missing values
knn_HWMDBMI <- train(HWMDBMI ~ . -CLINICID, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
# Estimate the missing values in HWMDBMI
study_data$HWMDBMI[is.na(study_data$HWMDBMI)] = predict(knn_HWMDBMI, study_data[is.na(study_data$HWMDBMI),])

# Repeat the above steps for the other variables with missing data
knn_LAB_BCD <- train(LAB_BCD ~ . -CLINICID, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BCD[is.na(study_data$LAB_BCD)] = predict(knn_LAB_BCD, study_data[is.na(study_data$LAB_BCD),])

knn_LAB_BHG <- train(LAB_BHG ~ . -CLINICID, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BHG[is.na(study_data$LAB_BHG)] = predict(knn_LAB_BHG, study_data[is.na(study_data$LAB_BHG),])

knn_SMK_12 <- train(SMK_12 ~ . -CLINICID, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$SMK_12[is.na(study_data$SMK_12)] = predict(knn_SMK_12, study_data[is.na(study_data$SMK_12),])

sum(is.na(study_data))  # 0 missing values


#----------------------------------------------------------------------------------------------------------------
ggplot(study_data, aes(cut_number(CLC_AGE, 20), fill = factor(HIGHBP))) + 
  geom_bar(position = position_fill(reverse = TRUE)) +
  xlab("Age") + 
  ylab("Proportion of HIGHBP") + 
  ggtitle("Age Affects HIGHBP in Three Distinct Sections") +
  scale_fill_discrete(name = "Has Hypertension?", labels = c("Yes", "No")) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))

ggplot(study_data, aes(cut_number(HWMDBMI, 20), fill = factor(HIGHBP))) + 
  geom_bar(position = position_fill(reverse = TRUE)) +
  xlab("HWMDBMI") + 
  ylab("Proportion of HIGHBP") + 
  ggtitle("Relationship between HWMDBMI and HIGHBP is Roughly Linear") +  
  scale_fill_discrete(name = "Has Hypertension?", labels = c("Yes", "No")) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))

ggplot(study_data, aes(cut_number(LAB_BCD, 20), fill = factor(HIGHBP))) + 
  geom_bar(position = position_fill(reverse = TRUE)) +
  xlab("LAB_BCD") + 
  ylab("Proportion of HIGHBP") +
  ggtitle("Relationship between LAB_BCD and HIGHBP is Roughly Linear") +
  scale_fill_discrete(name = "Has Hypertension?", labels = c("Yes", "No")) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))

ggplot(study_data, aes(cut_number(LAB_BHG, 20), fill = factor(HIGHBP))) + 
  geom_bar(position = position_fill(reverse = TRUE)) +
  xlab("LAB_BHG") + 
  ylab("Proportion of HIGHBP") +
  ggtitle("Relationship between LAB_BHG and HIGHBP is Roughly Linear") +
  scale_fill_discrete(name = "Has Hypertension?", labels = c("Yes", "No")) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1))


#----------------------------------------------------------------------------------------------------------------
# Make CLC_AGE into a categorical variable with three levels, and add this to our data frame
study_data = mutate(study_data, CLC_AGE_CAT = as.factor(cut_width(study_data$CLC_AGE, 20, boundary = 20, closed = "left")))

# Uses log odds link function
model_logit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE_CAT+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data)

# Uses inverse CDF of standard Normal distribution link function
model_probit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE_CAT+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'probit'), study_data)

# Uses complimentary log-log link function
model_cloglog = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE_CAT+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'cloglog'), study_data)


#----------------------------------------------------------------------------------------------------------------
# See that the minimum estimated probabilities that HIGHBP=2 are all greater than 10%
min(fitted.values(model_logit))
min(fitted.values(model_probit))
min(fitted.values(model_cloglog))

# The minimum estimated probabilities that HIGHBP=1 are also greater than 10%
max(fitted.values(model_logit))
max(fitted.values(model_probit))
max(fitted.values(model_cloglog))

# Calculate the Hosmer-Lemeshow statistic multiple times for each of the three models
pvalues_logit = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_logit$y, y = fitted.values(model_logit))["p.value",]

pvalues_probit = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_probit$y, y = fitted.values(model_probit))["p.value",]

pvalues_cloglog = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_cloglog$y, y = fitted.values(model_cloglog))["p.value",]

pvalues_logit = unlist(pvalues_logit)
pvalues_probit = unlist(pvalues_probit)
pvalues_cloglog = unlist(pvalues_cloglog)

# Plot distribution of Hosmer-Lemeshow statistics
ggplot(mapping = aes(x = factor(rep(c("Logit", "Probit", "Cloglog"),  each = 51), levels = c("Logit", "Probit", "Cloglog")),  y = c(pvalues_logit, pvalues_probit, pvalues_cloglog))) + 
  geom_boxplot() +
  xlab("Link Functions") + 
  scale_y_continuous("P-values", c(0.2, 0.4, 0.6, 0.8), limits = c(0, 1)) +
  ggtitle("Different Link Functions have Similiar P-values for Hosmer-Lemeshow Statistic") +
  theme(plot.title = element_text(hjust = 0.5))

#----------------------------------------------------------------------------------------------------------------
# Look at two models, one with CLC_AGE and one with CLC_AGE_CAT, and compare their AIC.
summary(model_logit)
summary(update(model_logit, ~ . -CLC_AGE_CAT+CLC_AGE))

# Calculate the Generalized Variance Inflation Factors (GVIF) for each predictor
vif(model_logit) 

# Build residual plot
ggplot(mapping = aes(model_logit$linear.predictors, residuals(model_logit, "pearson"))) + 
  geom_point() +
  xlab("Predicted Log Odds Ratio") + ylab("Pearson Residuals") +
  ggtitle("Residual Plot with No Replicates is Uninformative") +
  theme(plot.title = element_text(hjust = 0.5))

#----------------------------------------------------------------------------------------------------------------
drop1(model_logit, test = "LRT")

# Test whether the interaction terms between CLC_SEX and other variables are non-zero
model_gender = update(model_logit, ~ . + CLC_SEX:SMK_12 + CLC_SEX:CLC_AGE_CAT + CLC_SEX:HWMDBMI + CLC_SEX:LAB_BCD + CLC_SEX:LAB_BHG)
drop1(model_gender, test = "LRT")

# Test whether the interaction terms between CLC_AGE_CAT and other variables are non-zero
model_age = update(model_logit, ~ . + CLC_AGE_CAT:SMK_12 + CLC_AGE_CAT:CLC_SEX + CLC_AGE_CAT:HWMDBMI + CLC_AGE_CAT:LAB_BCD + CLC_AGE_CAT:LAB_BHG)
drop1(model_age, test = "LRT")