# Load some useful packages
library(tidyverse)  # Set of packages including ggplot2
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


# Alternative method: Use survival analysis to estimate missing values!
# Need to estimate both the distribution of the response, and also come up with an appropriate link function

# Plot non-missing and non-censored values of LAB_BHG
#ggplot(study_data[study_data$LAB_BHG < 900, ], aes(LAB_BHG)) + geom_histogram(binwidth = 2)

# Do the same for LAB_BCD
#ggplot(study_data[study_data$LAB_BCD < 900, ], aes(LAB_BCD)) + geom_histogram(binwidth = 1)

# Both of these look like they may roughly follow a Weibull distribution

# Plot weibull distributions with different parameters
#ggplot(mapping = aes(rweibull(3000, shape = 1, scale = 15))) + geom_histogram(binwidth = 2)

# Not done yet...

#simple_random_sample = svydesign(~1, variables = study_data[, 2:8], data = study_data)
#svykm(LAB_BHG ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+HIGHBP+LAB_BCD, simple_random_sample)


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


#----------------------------------------------------------------------------------------------------------------
# Use the package 'survey' to integrate weights into model
# Here, 'weights' represent the sample weights, and 'repweights' are the replicate weights 
glm_design = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
model_weighted = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = glm_design, family = quasibinomial(), data = study_data)


#----------------------------------------------------------------------------------------------------------------
# Draw conclusions about individual parameters for unweighted and weighted model

# Returns a dataframe with two columns: parameter estimates and std. errors for the given model
param_estimates = function(model) {
  return(as.data.frame(summary(model)$coefficients) %>%
           select('Estimate', 'Std. Error') %>%
           rownames_to_column("Variable") %>%
           mutate('p_value' = rep(NA, nrow(.))))
}

# For logit model, create dataframe containing MLE estimates and s.e's for parameters
logit_likelihood_tests = param_estimates(model_logit)

# Add likelihood ratio statistic p-values to the above dataframe
logit_likelihood_tests[3:8, 'p_value'] = as.data.frame(drop1(model_logit, test = "LRT"))[2:7, 'Pr(>Chi)']
logit_likelihood_tests



# Since it seems likely that a weighted maximum likelihood method is being used to estimate parameters (see p. 268 of Fitting 
# Regression Models to Survey Data), we cannot just use the regular likelihood ratio test. Use regTermTest() from survey 
# package to perform hypothesis testing using the Rao-Scott working likelihood ratio test instead (see p.273 of same paper)

# For weighted logit model, create dataframe containing MLE estimates and s.e's for parameters
weighted_likelihood_tests = param_estimates(model_weighted)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
weighted_likelihood_tests[3, 'p_value'] = regTermTest(model_weighted, 'SMK_12', method = 'LRT')$p
for(i in c(4:8)) {
  weighted_likelihood_tests[i, 'p_value'] = regTermTest(model_weighted, weighted_likelihood_tests$Variable[i], method = 'LRT')$p
}
weighted_likelihood_tests

# For both unweighted and weighted models, even when taking into account the other variables in the model, there is strong 
# evidence that sex, age, and body mass index influence average hypertension, and weaker evidence that blood mercury levels 
# (LAB_BHG) has an effect. It is inconclusive whether smoking status or blood cadmium levels (LAB_BCD) affect average hypertension.

# The p-values for weighted model do seem to be a bit more 'extreme' than those from unweighted model
# ie. large p-values from unweighted model tend to be even larger in weighted model, and small p-values from unweighted 
# model tend to be even smaller

# Notice that both unweighted and weighted models have similiar parameter estimates and std.error estimates,
# though std. errors for weighted model seem to be a bit larger


#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_SEX and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_SEX interaction terms
model_gender = update(model_logit, ~ . +CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
gender_likelihood_tests = param_estimates(model_gender)

# Add likelihood ratio statistic p-values to the above dataframe
gender_likelihood_tests[c(3:8, 10:14), 'p_value'] = as.data.frame(drop1(model_gender, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']
gender_likelihood_tests



# Repeat for weighted model
# Build weighted model with CLC_SEX interaction terms
model_gender_weighted = update(model_weighted, ~ .+CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
gender_weighted_likelihood_tests = param_estimates(model_gender_weighted)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
gender_weighted_likelihood_tests[3, 'p_value'] = regTermTest(model_gender_weighted, 'SMK_12', method = 'LRT')$p
gender_weighted_likelihood_tests[10, 'p_value'] = regTermTest(model_gender_weighted, 'SMK_12:CLC_SEX', method = 'LRT')$p
for(i in c(5:8, 11:14)) {
  gender_weighted_likelihood_tests[i, 'p_value'] = regTermTest(model_gender_weighted, gender_weighted_likelihood_tests$Variable[i], method = 'LRT')$p
}
gender_weighted_likelihood_tests


#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_AGE and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_AGE interaction terms
model_age = update(model_logit, ~ . +CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
age_likelihood_tests = param_estimates(model_age)

# Add likelihood ratio statistic p-values to the above dataframe
age_likelihood_tests[c(3:8, 10:14), 'p_value'] = as.data.frame(drop1(model_age, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']
age_likelihood_tests



# Repeat for weighted model
# Build weighted model with CLC_AGE interaction terms
model_age_weighted = update(model_weighted, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
age_weighted_likelihood_tests = param_estimates(model_age_weighted)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
age_weighted_likelihood_tests[3, 'p_value'] = regTermTest(model_age_weighted, 'SMK_12', method = 'LRT')$p
age_weighted_likelihood_tests[10, 'p_value'] = regTermTest(model_age_weighted, 'SMK_12:CLC_AGE', method = 'LRT')$p
for(i in c(4, 6:8, 11:14)) {
  age_weighted_likelihood_tests[i, 'p_value'] = regTermTest(model_age_weighted, age_weighted_likelihood_tests$Variable[i], method = 'LRT')$p
}
age_weighted_likelihood_tests


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

#----------------------------------------------------------------------------------------------------------------
# Experiments

# 1) Build different svyglm models using both sampling and replication weights, but using different kinds of replication weights
# Should get same point estimates, but different variances

test_design_1a = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
test_model_1a = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = test_design_1a, family = quasibinomial(), data = study_data)
summary(test_model_1a)

test_design_1b = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'BRR')
test_model_1b = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = test_design_1b, family = quasibinomial(), data = study_data)
summary(test_model_1b)

test_design_1c = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'other')
test_model_1c = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = test_design_1c, family = quasibinomial(), data = study_data)
summary(test_model_1c)

# And that is exactly what we see!

# 2) It was suggested that the sampling weights WGT_FULL might be measured in thousands or something, thus making them integers.
#    To test this, remember that these weights represent the number of Canadians that each person in the sample represent.
#    Thus, if a weight of 1 represents one person, the sum of these weights should equal about 30 million

sum(study_data$WGT_FULL)

# Which is what we get. So it seems like the weights are not scaled and are not integers. So using glm() with weights is
# not appropriate

# 3) As described on p. 25 of survey textbook, using type = 'bootstrap' in svyrepdesign() should automatically set scale ~ 1/500
# Notice that using type = 'other' and scale 1/499 gives same summary() results as type = 'bootstrap'

test_design_3a = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
test_model_3a = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = test_design_3a, family = quasibinomial(), data = study_data)
summary(test_model_3a)

test_design_3b = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'other', scale = 1/499, rscales = rep.int(1, 500))
test_model_3b = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = test_design_3b, family = quasibinomial(), data = study_data)
summary(test_model_3b)