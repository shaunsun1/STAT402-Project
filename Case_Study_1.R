# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2
library(caret)  # used for knn and cross-validation
library(ResourceSelection)  # hosmer-lemeshow test is used as a goodness-of-fit test

# Bring the chms_2018 dataset into R
study_data = read_csv("chms_2018.csv")


#----------------------------------------------------------------------------------------------------------------


# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to 
# be measured (ie. they are below the LOD)
# First understand how common this is
BCD_LOD = sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  
BCD_LOD  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)
BHG_LOD  # 999.5 appears 599 times in LAB_BHG

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is 
# called jittering.
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
study_data$LAB_BHG[study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)] = runif(BHG_LOD, 0, 2.1)

# Confirm that all values below the LOD have been replaced.
sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BCD
sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BHG


#----------------------------------------------------------------------------------------------------------------


# Use K-Nearest-Neighbors (KNN) to estimate missing values
# We've assumed here that the missing values have been determined randomly

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
ctrl <- trainControl(method="cv") 

# Using HWMDBMI as the response variable, find the best value of K using as the training set only the observations 
# that have no missing values
knn_HWMDBMI <- train(HWMDBMI ~ ., data = study_data[, seq(2, 8)], subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
# Estimate the missing values in HWMDBMI
study_data$HWMDBMI[is.na(study_data$HWMDBMI)] = predict(knn_HWMDBMI, study_data[is.na(study_data$HWMDBMI), seq(2, 8)])

# Repeat the above steps for the other variables with missing data
knn_LAB_BCD <- train(LAB_BCD ~ ., data = study_data[, seq(2, 8)], subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BCD[is.na(study_data$LAB_BCD)] = predict(knn_LAB_BCD, study_data[is.na(study_data$LAB_BCD), seq(2, 8)])

knn_LAB_BHG <- train(LAB_BHG ~ ., data = study_data[, seq(2, 8)], subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BHG[is.na(study_data$LAB_BHG)] = predict(knn_LAB_BHG, study_data[is.na(study_data$LAB_BHG), seq(2, 8)])

knn_SMK_12 <- train(as.factor(SMK_12) ~ ., data = study_data[, seq(2, 8)], subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$SMK_12[is.na(study_data$SMK_12)] = predict(knn_SMK_12, study_data[is.na(study_data$SMK_12), seq(2, 8)])

sum(is.na(study_data))  # 0 missing values


#----------------------------------------------------------------------------------------------------------------


# Now that all missing values are estimated, we can create a model that predicts HIGHBP
# Since we are interested only in determining which factors have a significant effect on HIGHBP 
# (given the other variables in the model), we can ignore interaction terms

# Model the response (HIGHBP) as having a Bernoulli distribution, and assume link function has the logit form
yfit = glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+as.factor(CLC_SEX)+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data)

# The null hypothesis for this test is that our proposed model fits the sample data as well as the saturated model
# The alternative hypothesis is that the saturated model fits the sample data better than the proposed model

hoslem.test(yfit$y, fitted.values(yfit), 10)

# Using a significance level of 0.05, we cannot reject the null hypothesis. This suggests that our model does 
# fit the sample data well

# Another test should be done before being satisfied with the fit of the model. Perhaps the residuals can be examined?
# Also, models with different link functions could be compared (ex. extreme value distribution, log-log function) 
# to see which one fits best.
# Additionally, even if our model does fit the sample data well, that does not mean it has predictive power 
# (ex. overfitting) This is a lower priority

# Check results of model
summary(yfit)


#----------------------------------------------------------------------------------------------------------------


# Test whether male and females share the same risk factors that affect HIGHBP
# Split observations by gender, and test whether the same risk factors hold for both males and females
male_fit<- glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_SEX == 1))
summary(male_fit) 

# For men, the only significant factor that influences HIGHBP (using a SL of 0.05) is CLC_AGE

female_fit<- glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_SEX == 2))
summary(female_fit) 

# For women, the significant factors are now CLC_AGE as well as HWMDBMI
# Thus, we see that women and men do indeed experience different risk factors

# Do different ages experience different risk factors for HIGHBP?
# Split the observations into three categories: Younger, Middle_Aged, and Older. Then test by category
younger_fit<- glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+as.factor(CLC_SEX)+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE < 40))
summary(younger_fit) 

# For the younger group, the significant factors are CLC_SEX and HWMDBMI (at a signficant level of 0.05)

middle_aged_fit<- glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+as.factor(CLC_SEX)+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE < 60&study_data$CLC_AGE >= 40))
summary(middle_aged_fit) 

# For the middle aged group, the only significant factor is now HWMDBMI

older_fit<- glm(as.factor(HIGHBP) ~ as.factor(SMK_12)+as.factor(CLC_SEX)+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE >= 60))
summary(older_fit) 

# For the older group, none of these factors are significant!
# Thus, different age groups do experience different risk factors for HIGHBP.

# We still need to do the analysis with the survey weights being taken into account. I think this can actually 
# be done using the regular glm() function