# Workshop on mixtures, Columbia August 2020

# install the packages needed for the workshop
#install.packages("gWQS")
# import the packages just installed
library(gWQS)


# define the path
#directory_path = "/Users/stefanorenzetti/Documents/Columbia Workshop Mixtures/"
#directory_path = "/Users/gennic01/desktop/temp teaching/mailman mixtures workshop 2018/"
#directory_path_out = "/Users/gennic01/desktop/temp teaching/mailman mixtures workshop 2020/"
# import the dataset
#dataset = read.csv(paste0(directory_path, "studypop.csv"))

dataset = read.csv("Data/studypop.csv")

# define the chemicals to include in the mixture
#mixture = c("LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA", "LBX180LA", "LBX187LA",
             "LBX194LA", "LBXD03LA", "LBXD05LA", "LBXD07LA", "LBXF03LA", "LBXF04LA", "LBXF05LA", "LBXF08LA",
             "LBXHXCLA", "LBXPCBLA")
## for workshop lecture example
mixture = c("LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA", "LBX180LA", "LBX187LA",
            "LBX194LA")

# log-transform the outcome
dataset$log_TELOMEAN = log(dataset$TELOMEAN)
summary(dataset)
# covariate only model
cov_only = glm(log_TELOMEAN ~ LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male,  data = dataset)
summary(cov_only)

# fit a first unadjusted model to look at the association between the mixture and the outcome
# TELOMEAN = Mean Telomere Length
results1 = gwqs(log_TELOMEAN ~ wqs, mix_name = mixture, data = dataset, q = 10, validation = 0.6, 
                b = 100, b1_pos = FALSE, b1_constr = TRUE, family = "gaussian", seed = 123)
# bar plot
gwqs_barplot(results1)
# scatter plot y vs wqs
gwqs_scatterplot(results1)
# scatter plot residuals vs fitted values
gwqs_fitted_vs_resid(results1)
#summary table
gwqs_summary_tab(results1)

summary(results1$fit)
results1$final_weights

# adjusting for covariates:
# blood data: LBXWBCSI LBXLYPCT LBXMOPCT LBXEOPCT LBXBAPCT LBXNEPCT
# demographics: age_cent age_sq race_cat bmi_cat3 ln_lbxcot edu_cat male
# positive direction
result2 = gwqs(log_TELOMEAN ~ wqs+LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)
gwqs_barplot(result2)
gwqs_summary_tab(result2)
gwqs_scatterplot(result2)
# number of positive betas
sum(result2$bres$b1 > 0)
# number of negative betas
sum(result2$bres$b1 < 0)

# negative direction
result3 = gwqs(log_TELOMEAN ~ wqs+LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, 
               data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = FALSE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)
gwqs_summary_tab(result3)
gwqs_barplot(result3)
gwqs_scatterplot(result3)
sum(result3$bres$b1 < 0)

#compare to random subset WQS
result4 = gwqs(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123, rs=TRUE, n_var=5) #3 for set of 9 and 5 for set of 18

gwqs_summary_tab(result4) 
gwqs_barplot(result4)
gwqs_scatterplot(result4)
gwqs_weights_tab(result4)
sum(result4$bres$b1 > 0)

# WQS with interaction
result2int = gwqs(log_TELOMEAN ~ wqs*male +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result2int) 
gwqs_barplot(result2int)
gwqs_scatterplot(result2int)
gwqs_weights_tab(result2int)


# run the wqs model using the stratified variables in the mixtures  
dataset$sex = factor(dataset$male) 
is.factor(dataset$male)
is.factor(dataset$sex) 

result5 = gwqs(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result5)  
gwqs_barplot(result5) 
gwqs_scatterplot(result5) 
gwqs_fitted_vs_resid(result5) 
gwqs_weights_tab(result5) 

# run the wqs model using the stratified variables in the mixtures with interaction
result6 = gwqs(log_TELOMEAN ~ wqs*male +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)

gwqs_summary_tab(result6) 
gwqs_scatterplot(result6) 
gwqs_fitted_vs_resid(result6) 
gwqs_barplot(result6) 

################################ takes FOREVER ##########################
#repeated holdout WQS
# positive direction
result5rh = gwqsrh(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
                 stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123, rh=20)
gwqs_summary_tab(result5rh)
gwqsrh_boxplot(result5rh)
gwqs_weights_tab(result5rh)
summary(result5rh$fit)
result5rh$final_weights


#write.csv(result5rh$wmat, paste0(directory_path_out, "WQSwtsrh.csv"))
