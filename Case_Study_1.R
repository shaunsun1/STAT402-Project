# We want to create a model that predicts whether a patient is hypertensive (HIGHBP = 1) 

# Reasons not to use simple subsitution when estimating values below the limit of detection (LOD) 
# (Talks about LOD in case study). More about this later on
# https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S0033845111067287 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3122101/
# http://analytics.ncsu.edu/sesug/2003/SD08-Croghan.pdf
# https://community.jmp.com/t5/JMPer-Cable/When-Responses-Are-Below-the-Limit-of-Detection/ba-p/28973

# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2 and dplyr

setwd("../Jonathan/R/Stat 402/Project") 

# First bring the chms_2018 dataset into R
study_data = read_csv("chms_2018.csv")

# The values for LAB_BCD and LAB_BHG are set to 999.5 when the patients recorded values are too low to be measured (ie. they are below the LOD)
# First understand how common this is
sum(study_data$LAB_BCD == 999.5)  # Occurs 54 times
sum(study_data$LAB_BHG == 999.5)  # Occurs 599 times 

# Since data set contains about 3000 observations, the proportion of times when LAB_BHG = 999.5 is quite large, and estimating these
# values incorrectly could have a big (negative) impact on our analysis

# Check how many missing values there are in this data set
sum(is.na(study_data))  # 175 missing values

# How to deal with these values below the LOD? Could estimate them all with a single value (ex. let them all equal 0 or LOD/2 or LOD/sqrt(2))
# This is very simple, but based on some reading (shown at the top), this approach does not seem recommended (makes estimators unbiased?)
# The alternative would be to somehow use the data points with LAB_BCD/LAB_BHG values above the LOD to estimate those values below the LOD
# For example, it sounds like the method of maximum likelihood can be used to estimate these values below the LOD
# Time series analysis can also apparently be used, but I don't know anything about this!


# After these values below the LOD are estimated, we can estimate the missing values in the data set (imputation)
# A simple approach is to replace missing values in continuous variables with the sample mean of that variable, and missing values in categorical variables with the mode (since mean would not make sense here)
# However, if the sample variance of LAB_BCD and LAB_BHG is not small (which it seems to be), then this is not an accurate approach 
# A better solution may be to replace missing values with randomly selected non-missing values of the same variable
# This would be easy to do, and I think would better take into account the variance of the data

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