# We want to create a model that predicts whether a patient is hypertensive (HIGHBP = 1) or not (HIGHBP = 0)

# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2
library(class)  # a package of recommended priority used for knn

# First bring the chms_2018 dataset into R
study_data = read_csv("chms_2018.csv")

# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to be measured (ie. they are below the LOD)
# First understand how common this is
BCD_LOD = sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  
BCD_LOD  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)
BHG_LOD  # 999.5 appears 599 times in LAB_BHG

# Since data set contains about 3000 observations, the proportion of times when LAB_BHG = 999.5 is quite large, and estimating these
# values incorrectly could have a big (negative) impact on our analysis

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is called jittering.
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
study_data$LAB_BHG[study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)] = runif(BHG_LOD, 0, 2.1)

# Confirm that all LOD values have been replaced.
sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BCD
sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BHG



#---------------------------------------------------------------------------------------------------------------
# Use K-Nearest-Neighbors (KNN) to estimate missing values

# Check how many missing values there are in this data set
sum(is.na(study_data))  # 175 missing values

# Return variables that contain missing values, as well as their respective number of missing values
n_missing = colSums(is.na(study_data)) # number of missing values in each column of data set
matrix(c(names(study_data)[n_missing > 0], n_missing[n_missing > 0]), nrow = 2, byrow = TRUE)

# Are there observations that contain multiple missing values?
multiple_missing = rowSums(is.na(study_data))  # number of missing values in each row of data set
sum(multiple_missing >= 2)  # number of observations that contain at least 2 missing values equals two
sum(multiple_missing == 2)  # number of observations that contain exactly 2 missing values equals two

# Before using KNN, the data should be standardized so that each variable in model is given equal weight
# This can be done with base::scale()

# Next, the value for K (# of neighbors) needs to be decided. This selection should be as rigorous as possible
# Cross-validation is a techique that estimates the "test error" (in contrast to "training error") of estimates using a particular model 
# Test error measures how well the model will perform on new data (which is what we care about!), while training error measures how well model fits training data (ie. original sample)
# The K value that produces the model that has the smallest estimate of test error will be used in our model

# Building the knn model can be done with class:knn() and I think cross-validation can be done with class::knn.cv()

# Now that we have our value of K, we predict our missing values 
# There are multiple variables that contain missing values however. One possible approach is to use an iterative procedure 
# For example, can start by imputing all missing values using some simple technique, such as replacement by sample mean of respective variable
# These can be considered initial guesses for the true values of the missing data
# At this point, the variables that had missing values are given some order, and one of the variables is chosen
# For this chosen variable, the previous guesses for the values of the missing data are replaced with the KNN estimates
# The next variable in line is chosen, and the process repeats until (hopeful!) convergence of the estimates of missing values



#------------------------------------------------------------------------------------------------------------------
# Now that all missing values are estimated, can finally create model that predicts HIGHBP
# Don't know if the rest of this works yet!

# A model is created using the simple replacement approach described above
count(study_data, SMK_12)  # SMK_12 has one NA value, and has a mode of 3
study_data = study_data %>%
  replace_na(list(SMK_12 = 3, HWMDBMI = mean(.$HWMDBMI, na.rm = TRUE), LAB_BCD = mean(.$LAB_BCD, na.rm = TRUE), LAB_BHG = mean(.$LAB_BHG, na.rm = TRUE)))
  
# Confirm that there are no more missing values in data set
sum(is.na(study_data))  # 0 missing values

# Will model the response (HIGHBP) as having a Bernoulli distribution, and assume link function has the logit form
model_full = glm(as.factor(HIGHBP) ~ as.factor(SMK_12)*as.factor(CLC_SEX)*CLC_AGE*HWMDBMI*LAB_BCD*LAB_BHG, binomial, study_data)
summary(model_full)

# The p-values given by the summary function tell us the significance of the variable given the other variables included in the model
# We could now reduce the variables included in the final model using backwards selection.
# There are however 2 flaws with this method:
#   1) Given there are so many variables and interaction terms in the full model, it would take a long time to perform this backwards selection
#   2) Performing so many tests increases both type 1 and type 2 error

# Thus, it may be a good idea to use ANOVA to test the significance of multiple variables/interaction terms at once
# Notice that the 'anova' function in R is capable of computing anaylsis of deviance tables for glms, which is what we learned in class

# Start by comparing the full model (every term and interaction term) with the model without interaction terms of "length" 5 or 6
model_1 = glm(as.factor(HIGHBP) ~ (as.factor(SMK_12)+as.factor(CLC_SEX)+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG)^4, binomial, study_data)
anova(model_1, model_full, test = "Chisq")

# Repeat until final model has been selected

# Notice the errors we get when model_full and model_1 are created. I'm not sure why we are getting these errors, but it may go away
# after properly estimating the values below the LOD and the missing values