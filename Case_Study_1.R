# Load some useful packages
library(tidyverse)  # a set of packages including ggplot2
library(caret)  # used for knn and cross-validation
library(car)  # calculates VIF to check for multicollinearity
library(survey)  # weighted glm
library(ResourceSelection)  # Hosmer-Lemeshow goodness-of-fit test

# Bring the chms_2018 dataset into R, and remove the bootstrap weights as they are not used in this analysis
study_data = read.csv("chms_2018.csv")
study_data = select(study_data, -starts_with("BS"))

# Make SMK_12 a factor
study_data$SMK_12 = as.factor(study_data$SMK_12)

# Choose seeds so that results from imputing censored data is fixed.
# These numbers were drawn from a Uniform distribution with boundaries (1, 10^9). This will not be truly random, but will suffice for our purposes
seed1 = 336741157
seed2 = 874515716


#----------------------------------------------------------------------------------------------------------------
# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to 
# be measured (ie. they are below the LOD)
# First understand how common this is
BCD_LOD = sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  
BCD_LOD  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)
BHG_LOD  # 999.5 appears 599 times in LAB_BHG

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is called jittering.
set.seed(seed1)
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
set.seed(seed2)
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
ctrl <- trainControl(method="cv", number = 5) 

# Using HWMDBMI as the response variable, find the best value of K using as the training set only the observations 
# that have no missing values
knn_HWMDBMI <- train(HWMDBMI ~ . -CLINICID-WGT_FULL, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
# Estimate the missing values in HWMDBMI
study_data$HWMDBMI[is.na(study_data$HWMDBMI)] = predict(knn_HWMDBMI, study_data[is.na(study_data$HWMDBMI),])

# Repeat the above steps for the other variables with missing data
knn_LAB_BCD <- train(LAB_BCD ~ . -CLINICID-WGT_FULL, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BCD[is.na(study_data$LAB_BCD)] = predict(knn_LAB_BCD, study_data[is.na(study_data$LAB_BCD),])

knn_LAB_BHG <- train(LAB_BHG ~ . -CLINICID-WGT_FULL, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BHG[is.na(study_data$LAB_BHG)] = predict(knn_LAB_BHG, study_data[is.na(study_data$LAB_BHG),])

knn_SMK_12 <- train(SMK_12 ~ . -CLINICID-WGT_FULL, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$SMK_12[is.na(study_data$SMK_12)] = predict(knn_SMK_12, study_data[is.na(study_data$SMK_12),])

sum(is.na(study_data))  # 0 missing values

study_data$HIGHBP = study_data$HIGHBP - 1  # change levels for simplicity

#----------------------------------------------------------------------------------------------------------------
# Now that all missing values are estimated, we can create a model that predicts HIGHBP
# Since we are interested only in determining which factors have a significant effect on HIGHBP 
# (given the other variables in the model), we can ignore interaction terms

# Three models are built each using a different link function. All of them model the response (HIGHBP) as a bernoulli random variable

# Assumes log odds link function
model_logit = glm(HIGHBP ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data)

# Assumes inverse CDF of standard Normal distribution link function
model_probit = glm(HIGHBP ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'probit'), study_data)

# Assumes complimentary log-log link function
model_cloglog = glm(HIGHBP ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'cloglog'), study_data)

# To check how well these models fit the sample data, the Deviance or the Chi-Squared statistic would normally be used
# However, our model includes continuous explanatory variables, which means these statistics will not follow their 
# theoretical asymptotic distributions. To get around this, the Hosmer-Lemeshow statistic is used instead, 
# as suggested by our class textbook (page 136). 

# Unfortunately, the Hosmer-Lemeshow statistic is sensitive to the number of groups the data is split into
# Thus, we calculate this statistic repeatedly, letting the number of groups vary between 10 and 60
# If this statistic is to be used at all, it seems like a process simliliar to this would be necessary

# Calculate the Hosmer-Lemeshow statistic multiple times for each of the three models
pvalues_logit = sapply(X = seq(10, 50, 1), FUN = hoslem.test, x = model_logit$y, y = fitted.values(model_logit))["p.value",]
pvalues_probit = sapply(X = seq(10, 50, 1), FUN = hoslem.test, x = model_probit$y, y = fitted.values(model_probit))["p.value",]
pvalues_cloglog = sapply(X = seq(10, 50, 1), FUN = hoslem.test, x = model_cloglog$y, y = fitted.values(model_cloglog))["p.value",]

pvalues_logit = unlist(pvalues_logit)
pvalues_probit = unlist(pvalues_probit)
pvalues_cloglog = unlist(pvalues_cloglog)

# Plot distribution of Hosmer-Lemeshow statistics
ggplot(mapping = aes(x = "Logit",  y = pvalues_logit)) + geom_boxplot()
ggplot(mapping = aes(x = "Probit", y = pvalues_probit)) + geom_boxplot()
ggplot(mapping = aes(x = "Cloglog", y = pvalues_cloglog)) + geom_boxplot()

# As can be seen from the graphs, all three models produce similiar Hosmer-Lemeshow statistics, suggesting that they all 
# fit the sample about as well. Additionally, there AIC values are all quite similiar, again suggesting that they all fit 
# roughly the same. However, most of the p-values are quite small, providing evidence that none of these models 
# fit the sample data as well as the saturated model.

# This is concerning, although maybe not surprising because no interaction terms or higher order terms have been fitted
# 

# Residuals would usually be checked as well. However...

# Since the Logit model gives the highest pvalues and also the smallest AIC value, we will use exclusively this model
# for the rest of the analysis


# Calculate the variance inflation factors (VIF) for each predictor to check for multicollinearity
vif(model_logit)  # These values are all close to 1, so we do not need to worry about multicollinearity

# To check whether a variable has a "significant" effect on the response HIGHBP, we remove each variable in turn 
# and compare the difference in deviance statistics. Even with sparse data such as ours, this difference in 
# deviance should still approximate a chi-square distribution with 1 d.f

drop1(model_logit, test = "LRT")


# These p-values give another perspective on the strength of the relationship between a regressor and HIGHBP.
# We sort the variables from those with the smallest to largest p-values: CLC_AGE, HWBDBMI, CLC_SEX, LAB_BHG, and then SMK_12 or LAB_BCD.


# Using a significance level of 0.05, we cannot reject the null hypothesis' that SMK_12 and LAB_BCD affect HIGHBP,
# and thus leave them out of our model on the grounds of parsimony


#----------------------------------------------------------------------------------------------------------------
# Test whether male and females share the same risk factors that affect HIGHBP
# Split observations by gender, and test whether the same risk factors hold for both males and females
model_male<- glm(HIGHBP ~ SMK_12+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_SEX == 1))
 drop1(model_male, test = "LRT")


model_female<- glm(HIGHBP ~ SMK_12+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_SEX == 2))
drop1(model_female, test = "LRT")

# For women, the significant factors are now CLC_AGE as well as HWMDBMI
# Thus, we see that women and men do indeed experience different risk factors




# Do different ages experience different risk factors for HIGHBP?
# Split the observations into three categories: Younger, Middle_Aged, and Older. Then test by category
model_younger<- glm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE < 40))
drop1(model_younger, test = "LRT")
# For the younger group, the significant factors are CLC_SEX and HWMDBMI (at a signficant level of 0.05)

model_middleAged<- glm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE >= 40&study_data$CLC_AGE < 60))
drop1(model_middleAged, test = "LRT")

# For the middle aged group, the only significant factor is now HWMDBMI

model_older<- glm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data, subset = (study_data$CLC_AGE >= 60))
drop1(model_older, test = "LRT")

# For the older group, none of these factors are significant!
# Thus, different age groups do experience different risk factors for HIGHBP.




#----------------------------------------------------------------------------------------------------------------
# Redo analysis with survey weights
# A quasibinomial model is fit to avoid an error about non-integer responses. This solution was recommended by the help files
weighted_design = svydesign(ids = ~1, data = study_data, weights = study_data$WGT_FULL)  # used in design variable below
model_weighted_logit = svyglm(HIGHBP ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data)
drop1(model_weighted_logit, test = "LRT")

# Test whether male and females share the same risk factors that affect HIGHBP
# Split observations by gender, and test whether the same risk factors hold for both males and females
model_weighted_male = svyglm(HIGHBP ~ SMK_12+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data, subset = (study_data$CLC_SEX == 1))
drop1(model_weighted_male, test = "LRT")

model_weighted_female = svyglm(HIGHBP ~ SMK_12+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data, subset = (study_data$CLC_SEX == 2))
drop1(model_weighted_female, test = "LRT")

# For women, the significant factors are now CLC_AGE as well as HWMDBMI
# Thus, we see that women and men do indeed experience different risk factors




# Do different ages experience different risk factors for HIGHBP?
# Split the observations into three categories: Younger, Middle_Aged, and Older. Then test by category
model_weighted_younger = svyglm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data, subset = (study_data$CLC_AGE < 40))
drop1(model_weighted_younger, test = "LRT")
# For the younger group, the significant factors are CLC_SEX and HWMDBMI (at a signficant level of 0.05)

model_weighted_middleAged = svyglm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data, subset = (study_data$CLC_AGE >= 40&study_data$CLC_AGE < 60))
drop1(model_weighted_middleAged, test = "LRT")

# For the middle aged group, the only significant factor is now HWMDBMI

model_weighted_older = svyglm(HIGHBP ~ SMK_12+CLC_SEX+HWMDBMI+LAB_BCD+LAB_BHG, design = weighted_design, family = quasibinomial(link = 'logit'), data = study_data, subset = (study_data$CLC_AGE >= 60))
drop1(model_weighted_older, test = "LRT")

# For the older group, none of these factors are significant!
# Thus, different age groups do experience different risk factors for HIGHBP.


