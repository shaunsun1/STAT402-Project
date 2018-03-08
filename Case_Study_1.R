# We want to create a model that predicts whether a patient is hypertensive (HIGHBP = 1) 

# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2

# First bring the chms_2018 dataset into R
study_data = read_csv("chms_2018.csv")

# The values for LAB_BCD and LAB_BHG are set to 999.5 when the patients recorded values are too low to be measured (ie. they are below the LOD)
# First understand how common this is
BCD_LOD = sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  
BCD_LOD  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)
BHG_LOD  # 999.5 appears 599 times in LAB_BHG

# Since data set contains about 3000 observations, the proportion of times when LAB_BHG = 999.5 is quite large, and estimating these
# values incorrectly could have a big (negative) impact on our analysis

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD) as per the professor's suggestion; this is called jittering.
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
study_data$LAB_BHG[study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)] = runif(BHG_LOD, 0, 2.1)

# Confirm that all LOD values have been replaced.
sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BCD
sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BHG



#---------------------------------------------------------------------------------------------------------------
# Dealing with missing values (work in progress)

# Check how many missing values there are in this data set
sum(is.na(study_data))  # 175 missing values

# Return variables that contain missing values, as well as their respective number of missing values
n_missing = colSums(is.na(study_data)) # number of missing values in each column of data set
matrix(c(names(study_data)[n_missing > 0], n_missing[n_missing > 0]), nrow = 2, byrow = TRUE)

# One of these is a categorical variable with 3 levels. Others are continuous variables
# How should we model these variables?

# Start by ploting HWMDBMI. This may help us figure out how to model it (error message is okay, it just refers to the na values)
ggplot(study_data, aes(HWMDBMI)) + geom_histogram(binwidth = 0.1, boundary = 0)

# This sampling distribution seems very odd!?! Why does it look like this?
# Perhaps the missing values are not random at all, but explain the 'dips' of this graph?
# Hmm, on the other hand, the 79 missing observations cannot affect this graph too much since there are ~3000 observations
# This is how it would basically look even with those missing values

# How does this affect analysis?

# Plot LAB_BCD. Notice how LAB_BCD exclusively takes on integer values when > 10
ggplot(study_data, aes(LAB_BCD)) + geom_histogram(binwidth = 0.5, boundary = 0)

#Plot LAB_BHG
ggplot(study_data, aes(LAB_BHG)) + geom_histogram(binwidth = 0.5, boundary = 0)

# Will try using non-parametric techniques to model data. This allows us to avoid guessing the distributions of the explanatory variables
# 


#------------------------------------------------------------------------------------------------------------------
# Once missing values are estimated...
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