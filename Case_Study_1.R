# Load some useful packages
library(tidyverse)  # Set of packages including ggplot2
library(survival)  # Survival analysis
library(caret)  # Used for knn and cross-validation
library(car)  # Calculates VIF to check for multicollinearity
library(ResourceSelection)  # Hosmer-Lemeshow goodness-of-fit test
library(survey)  # Build glm using weights

# Bring the chms_2018 dataset into R
study_data = read.csv("chms_2018.csv")

# Make SMK_12 a factor
study_data$SMK_12 = as.factor(study_data$SMK_12)

# Choose seeds so that results from imputing censored data is fixed.
# These numbers were drawn from a Uniform distribution with boundaries (1, 10^9). This will not be truly random, but will suffice for our purposes
seed1 = 336741157
seed2 = 874515716


#----------------------------------------------------------------------------------------------------------------
# Estimate censored data

# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to be measured (ie. they are below the LOD). 
# First understand how common this is
BCD_LOD = study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)
sum(BCD_LOD)  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)
sum(BHG_LOD)  # 999.5 appears 599 times in LAB_BHG

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is called jittering
set.seed(seed1)
study_data$LAB_BCD[study_data$LAB_BCD == 999.5 & !is.na(study_data$LAB_BCD)] = runif(BCD_LOD, 0, 0.71)
set.seed(seed2)
study_data$LAB_BHG[study_data$LAB_BHG == 999.5 & !is.na(study_data$LAB_BHG)] = runif(BHG_LOD, 0, 2.1)




# Alternative method: Use survival analysis to estimate missing values!
# Need to estimate both the distribution of the response, and also come up with an appropriate link function

# Plot non-missing and non-censored values of LAB_BHG
#ggplot(study_data[study_data$LAB_BHG < 900, ], aes(LAB_BHG)) + geom_histogram(binwidth = 2)

# Do the same for LAB_BCD
#ggplot(study_data[study_data$LAB_BCD < 900, ], aes(LAB_BCD)) + geom_histogram(binwidth = 1)

# Both of these look like they may roughly follow a Weibull distribution

# Plot weibull distributions with different parameters
#ggplot(mapping = aes(rweibull(3000, shape = 1, scale = 15))) + geom_histogram(binwidth = 2)


# Create a survival object for LAB_BHG
# Let start_time_BHG equal study_data$LAB_BHG unless LAB_BHG equals 999.5 (ie. unless BHG_LOD = TRUE). Then, should equal NA
#start_time_BHG = study_data$LAB_BHG
#is.na(start_time_BHG) = BHG_LOD

# Let end_time_BHG equal study_data$LAB_BHG unless LAB_BHG equals 999.5. Then, should equal 2.1
#end_time_BHG = study_data$LAB_BHG
#end_time_BHG[BHG_LOD] = 2.1

#surv_object_BHG = Surv(time = start_time_BHG, time2 = end_time_BHG, type = 'interval2')

#study_data$surv_object_BHG = surv_object_BHG  # Add surv_object_BHG to dataframe
#study_data = study_data[, c(1:8, 510, 9:509)]


# Use Kaplan-Meier algorithm to estimate Survival function for LAB_BHG
#survival_BHG = survfit(surv_object_BHG ~ 1, type = "kaplan-meier")

# Plot -log(Kaplan-Meier) estimate against response. A straight line is added as reference 
#ggplot(mapping = aes(survival_BHG$time, survival_BHG$surv)) + geom_point() + geom_smooth(method = "lm", se = FALSE)

# This is not a straight line, which suggests that the response cannot be approximated well with an exponential distribution

# Model LAB_BHG as a Weibull distribution and use log(?) link function
# Use as training set all observations with no missing data or LAB_BCD censored data
#model_surv_BHG = survreg(surv_object_BHG ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+HIGHBP+LAB_BCD, study_data, subset = (rowSums(is.na(study_data)) == 0) & !BCD_LOD)

# Use model to predict censored LAB_BHG data
#predict(model_surv_BHG, study_data[BHG_LOD, ])
#ggplot(mapping = aes(predict(model_surv_BHG, study_data[BHG_LOD, ]))) + geom_histogram()


#surv_design = svrepdesign(study_data[, c(2:7, 9)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
#new_design = svydesign(~1, data = study_data)
#svysurvreg(HIGHBP ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+surv_object_BHG, design = new_design)




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
knn_HWMDBMI <- train(HWMDBMI ~ as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE+LAB_BCD+LAB_BHG, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
# Estimate the missing values in HWMDBMI
study_data$HWMDBMI[is.na(study_data$HWMDBMI)] = predict(knn_HWMDBMI, study_data[is.na(study_data$HWMDBMI),])

# Repeat the above steps for the other variables with missing data
knn_LAB_BCD <- train(LAB_BCD ~ HWMDBMI+as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE+LAB_BHG, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BCD[is.na(study_data$LAB_BCD)] = predict(knn_LAB_BCD, study_data[is.na(study_data$LAB_BCD),])

knn_LAB_BHG <- train(LAB_BHG ~ LAB_BCD+HWMDBMI+as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BHG[is.na(study_data$LAB_BHG)] = predict(knn_LAB_BHG, study_data[is.na(study_data$LAB_BHG),])

knn_SMK_12 <- train(SMK_12 ~ LAB_BHG+LAB_BCD+HWMDBMI+as.factor(HIGHBP)+CLC_SEX+CLC_AGE, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$SMK_12[is.na(study_data$SMK_12)] = predict(knn_SMK_12, study_data[is.na(study_data$SMK_12),])

sum(is.na(study_data))  # 0 missing values


#----------------------------------------------------------------------------------------------------------------
# Do some EDA

# Check correlations between continuous explanatory variables (CLC_AGE, HWMDBMI, LAB_BCD, and LAB_BHG)
cor(study_data[, 4:8])
pairs(study_data[2:8])

# Variables are not linearly correlated. Thus, we can maybe ignore other explanatory variables when plotting log odds against
# a single variable

# Group by categorical CLC_AGE and calculate log odds
eda = study_data[, 1:8] %>%
  mutate(CLC_AGE_CAT = cut_number(CLC_AGE, 30)) %>%
  group_by(CLC_AGE_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))))

# Plot logit(HIGHBP) against our categorical CLC_AGE
ggplot(eda, aes(CLC_AGE_CAT, logit)) +
  geom_point() 

# No higher order terms for CLC_AGE are needed

# Repeat for HWMDBMI
eda = study_data[, 1:8] %>%
  mutate(HWMDBMI_CAT = cut_number(HWMDBMI, 30)) %>%
  group_by(HWMDBMI_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))))

# Plot logit(HIGHBP) against our categorical HWMDBMI
ggplot(eda, aes(HWMDBMI_CAT, logit)) +
  geom_point() 

# A quadratic term for HWMDBMI should be added

# Repeat for LAB_BCD
eda = study_data[, 1:8] %>%
  mutate(LAB_BCD_CAT = cut_number(LAB_BCD, 20)) %>%
  group_by(LAB_BCD_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))))

# Plot logit(HIGHBP) against our categorical HWMDBMI
ggplot(eda, aes(LAB_BCD_CAT, logit)) +
  geom_point()  

# A quadratic term for LAB_BCD should be added

# Repeat for LAB_BHG
eda = study_data[, 1:8] %>%
  mutate(LAB_BHG_CAT = cut_number(LAB_BHG, 30)) %>%
  group_by(LAB_BHG_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))))

# Plot logit(HIGHBP) against our categorical HWMDBMI
ggplot(eda, aes(LAB_BHG_CAT, logit)) +
  geom_point() 

# A quadratic term for LAB_BHG could be added


#----------------------------------------------------------------------------------------------------------------
# Build three models each using a different link function

# Uses log odds link function
model_logit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data)

# Uses inverse CDF of standard Normal distribution link function
model_probit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'probit'), study_data)

# Uses complimentary log-log link function
model_cloglog = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'cloglog'), study_data)


#----------------------------------------------------------------------------------------------------------------
# Hosmer-Lemeshow statistic

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
# Continue model adequacy checking

# Calculate the Generalized Variance Inflation Factors (GVIF) for each predictor
vif(model_logit) 

# Build residual plot
ggplot(mapping = aes(model_logit$linear.predictors, residuals(model_logit, "pearson"))) + 
  geom_point() +
  xlab("Predicted Log Odds Ratio") + ylab("Pearson Residuals") +
  ggtitle("Residual Plot with No Replicates is Uninformative") +
  theme(plot.title = element_text(hjust = 0.5))

# Add some quadratic terms to logit model
model_logit_quad = update(model_logit, ~ .+I(HWMDBMI^2)+I(LAB_BCD^2)+I(LAB_BHG^2))


#----------------------------------------------------------------------------------------------------------------
# Use the package 'survey' to integrate weights into model
# Here, 'weights' represent the sample weights, and 'repweights' are the replicate weights 
glm_design = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
model_weighted = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = glm_design, family = quasibinomial(), data = study_data)

# Now add a few quadratic terms
model_weighted_quad = update(model_weighted, ~ .+I(HWMDBMI^2)+I(LAB_BCD^2)+I(LAB_BHG^2))


#----------------------------------------------------------------------------------------------------------------
# Draw conclusions about individual parameters for unweighted and weighted model

# For the given model, returns a dataframe with three columns (parameter estimates, std. errors, and p-values) 
param_estimates = function(model) {
  return(as.data.frame(summary(model)$coefficients) %>%
           select('Estimate', 'Std. Error') %>%
           rownames_to_column("Variable") %>%
           mutate('p_value' = rep(NA, nrow(.))))
}

# For the given model, adds to the given dataframe three columns 
# (weighted parameter estimates, weighted std. errors, and weighted p-values) 
weighted_param_estimates = function(model, data_frame) {
  weighted_data_frame = as.data.frame(summary(model)$coefficients) %>%
                    select('Estimate', 'Std. Error') %>%
                    rownames_to_column("Variable") %>%
                    rename('Weighted Estimate' = 'Estimate', 'Weighted Std. Error' = 'Std. Error') %>%
                    mutate('Weighted p_value' = rep(NA, nrow(.)))
  
  complete_data_frame = full_join(data_frame, weighted_data_frame, 'Variable')
  complete_data_frame = complete_data_frame[, c(1, 2, 5, 3, 6, 4, 7)]
}

# For logit model, create dataframe containing MLE estimates and s.e's for parameters
logit_likelihood_tests = param_estimates(model_logit)

# Add likelihood ratio statistic p-values to the above dataframe
logit_likelihood_tests[3:8, 'p_value'] = as.data.frame(drop1(model_logit, test = "LRT"))[2:7, 'Pr(>Chi)']


# Since it seems likely that a weighted maximum likelihood method is being used to estimate parameters (see p. 268 of Fitting 
# Regression Models to Survey Data), we cannot just use the regular likelihood ratio test. Use regTermTest() from survey 
# package to perform hypothesis testing using the Rao-Scott working likelihood ratio test instead (see p.273 of same paper)

# For weighted logit model, create dataframe containing MLE estimates and s.e's for parameters
logit_likelihood_tests = weighted_param_estimates(model_weighted, logit_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
logit_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_weighted, 'SMK_12', method = 'LRT')$p
for(i in c(4:8)) {
  logit_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_weighted, logit_likelihood_tests$Variable[i], method = 'LRT')$p
}
logit_likelihood_tests

# For both unweighted and weighted models, even when taking into account the other variables in the model, there is strong 
# evidence that sex, age, and body mass index influence average hypertension, and weaker evidence that blood mercury levels 
# (LAB_BHG) has an effect. It is inconclusive whether smoking status or blood cadmium levels (LAB_BCD) affect average hypertension.

# The p-values for weighted model do seem to be a bit more 'extreme' than those from unweighted model
# ie. large p-values from unweighted model tend to be even larger in weighted model, and small p-values from unweighted 
# model tend to be even smaller

# Notice that both unweighted and weighted models have similiar parameter estimates and std.error estimates,
# though std. errors for weighted model seem to be a bit larger

# Repeat with quadratic terms
# First unweighted logit model
logit_quad_likelihood_tests = param_estimates(model_logit_quad)
logit_quad_likelihood_tests[3:11, 'p_value'] = as.data.frame(drop1(model_logit_quad, test = "LRT"))[2:10, 'Pr(>Chi)']


# Now weighted model
logit_quad_likelihood_tests = weighted_param_estimates(model_weighted_quad, logit_quad_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
logit_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_weighted_quad, 'SMK_12', method = 'LRT')$p
for(i in c(4:11)) {
  logit_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_weighted_quad, logit_quad_likelihood_tests$Variable[i], method = 'LRT')$p
}
logit_quad_likelihood_tests


#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_SEX and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_SEX interaction terms
model_gender = update(model_logit, ~ . +CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
gender_likelihood_tests = param_estimates(model_gender)

# Add likelihood ratio statistic p-values to the above dataframe
gender_likelihood_tests[c(3:8, 10:14), 'p_value'] = as.data.frame(drop1(model_gender, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']


# Repeat for weighted model
# Build weighted model with CLC_SEX interaction terms
model_gender_weighted = update(model_weighted, ~ .+CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
gender_likelihood_tests = weighted_param_estimates(model_gender_weighted, gender_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
gender_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_gender_weighted, 'SMK_12', method = 'LRT')$p
gender_likelihood_tests[10, 'Weighted p_value'] = regTermTest(model_gender_weighted, 'SMK_12:CLC_SEX', method = 'LRT')$p
for(i in c(5:8, 11:14)) {
  gender_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_gender_weighted, gender_likelihood_tests$Variable[i], method = 'LRT')$p
}
gender_likelihood_tests


# Repeat with quadratic terms
# First unweighted model with sex interaction terms
model_gender_quad = update(model_logit_quad, ~ . +CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)
gender_quad_likelihood_tests = param_estimates(model_gender_quad)
gender_quad_likelihood_tests[c(3:11, 13:17), 'p_value'] = as.data.frame(drop1(model_gender_quad, scope = ~ ., test = "LRT"))[2:15, 'Pr(>Chi)']

# Now weighted model
model_gender_weighted_quad = update(model_weighted_quad, ~ .+CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)
gender_quad_likelihood_tests = weighted_param_estimates(model_gender_weighted_quad, gender_quad_likelihood_tests)

gender_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_gender_weighted_quad, 'SMK_12', method = 'LRT')$p
gender_quad_likelihood_tests[13, 'Weighted p_value'] = regTermTest(model_gender_weighted_quad, 'SMK_12:CLC_SEX', method = 'LRT')$p
for(i in c(5:11, 14:17)) {
  gender_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_gender_weighted_quad, gender_quad_likelihood_tests$Variable[i], method = 'LRT')$p
}
gender_quad_likelihood_tests


#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_AGE and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_AGE interaction terms
model_age = update(model_logit, ~ . +CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
age_likelihood_tests = param_estimates(model_age)

# Add likelihood ratio statistic p-values to the above dataframe
age_likelihood_tests[c(3:8, 10:14), 'p_value'] = as.data.frame(drop1(model_age, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']


# Repeat for weighted model
# Build weighted model with CLC_AGE interaction terms
model_age_weighted = update(model_weighted, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
age_likelihood_tests = weighted_param_estimates(model_age_weighted, age_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
age_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_age_weighted, 'SMK_12', method = 'LRT')$p
age_likelihood_tests[10, 'Weighted p_value'] = regTermTest(model_age_weighted, 'SMK_12:CLC_AGE', method = 'LRT')$p
for(i in c(4, 6:8, 11:14)) {
  age_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_age_weighted, age_likelihood_tests$Variable[i], method = 'LRT')$p
}
age_likelihood_tests


# Repeat with quadratic terms
# First unweighted model with age interaction terms
model_age_quad = update(model_logit_quad, ~ . +CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
age_quad_likelihood_tests = param_estimates(model_age_quad)
age_quad_likelihood_tests[c(3:11, 13:17), 'p_value'] = as.data.frame(drop1(model_age_quad, scope = ~ ., test = "LRT"))[2:15, 'Pr(>Chi)']

# Now weighted model
model_age_weighted_quad = update(model_weighted_quad, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
age_quad_likelihood_tests = weighted_param_estimates(model_age_weighted_quad, age_quad_likelihood_tests)

age_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, 'SMK_12', method = 'LRT')$p
age_quad_likelihood_tests[13, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, 'SMK_12:CLC_AGE', method = 'LRT')$p
for(i in c(4, 6:11, 14:17)) {
  age_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, age_quad_likelihood_tests$Variable[i], method = 'LRT')$p
}
age_quad_likelihood_tests

#----------------------------------------------------------------------------------------------------------------
# Some Questions to Think About

# 1) In svrepdesign function, are the replication weights (our bootstrap weights) really of type 'bootstrap'? (see help page)
#    Need to better understand what Canadian Health Measures Survey (CHMS) means by 'bootstrap weights', and what author of survey package means by 'boostrap' replication weights

#    First talk about how variance is calculated (p. 240 and in section 2.3 of text)
#    From 'The Rao-Wu Rescaling Bootstrap: From theory to practice' (on github), seems likely that StatsCan uses the Rao-Wu rescaled bootstrap for the CHMS
#    as.svyrepdesign() seems to confirm that type = 'bootstrap' in svyrepdesign() can refer to Rao and Wu's rescaled boostrap
#    Section 2.3 mentions what scale in svrepdesign() should be for bootstrap. From experiments below, scale seems to equal 500-1 for bootstrap method
#    This agrees with 'Rescaled Bootstrap for stratified multistage sampling'

# 2) In svrepdesign function, not specifying repweights produces the warning: 'You must provide replication weights' 
#    Why is this?

#    I'm guessing that since replication weights are needed to estimate the variance of estimators, instead of just returning
#    point estimates without any kind of measure of error (which would be useless), the function just creates a warning

# 3) In svyglm function, letting family = binomial() produces the warning: 'In eval(family$initialize) : non-integer #successes in a binomial glm!'
#    Both the svyglm help file and the survey author's textbook (p.110) recommend using family = quasibinomial() to avoid these warnings, but no clear rationale is given

# 4) According to survey author's textbook (p.98), svyglm() does not seem to use maximum likelihood to estimate parameters,
#    and author concludes that we can't use likelihood ratio tests to compare nested models, but should use Wald statistic instead
#    Not sure if this is just for linear regression, or for logistic regression as well.

#    Actually, I think svyglm() uses pseudo-likelihood estimation (see p. 268 of Fitting Regression Models to Survey Data)
