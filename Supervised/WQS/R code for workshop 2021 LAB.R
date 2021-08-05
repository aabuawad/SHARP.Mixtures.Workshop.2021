# Workshop on mixtures, Columbia August 2021

#rm(list = ls())
# install the packages needed for the workshop
install.packages("gWQS")
# import the packages just installed
library(gWQS)

# define the path 
directory_path = "/Users/gennic01/desktop/temp teaching/mailman mixtures workshop 2018/"
directory_path_out = "/Users/gennic01/desktop/temp teaching/mailman mixtures workshop 2021/"
# import the dataset
dataset = read.csv(paste0(directory_path_out, "studypop.csv"))
dim(dataset)

# define the chemicals to include in the mixture
mixture=names(dataset)[grep("LA",names(dataset))]
#mixture = c("LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA", "LBX180LA", "LBX187LA",
#             "LBX194LA", "LBXD03LA", "LBXD05LA", "LBXD07LA", "LBXF03LA", "LBXF04LA", "LBXF05LA", "LBXF08LA",
#             "LBXHXCLA", "LBXPCBLA")
mixture
## for workshop lecture example
#mixture = c("LBX074LA", "LBX099LA", "LBX118LA", "LBX138LA", "LBX153LA", "LBX170LA", "LBX180LA", "LBX187LA","LBX194LA")

# log-transform the outcome
dataset$log_TELOMEAN = log(dataset$TELOMEAN)
summary(dataset)
# covariate only model
cov_only = glm(log_TELOMEAN ~ LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                as.factor(race_cat) + as.factor(bmi_cat3) + ln_lbxcot + as.factor(edu_cat) + male,  data = dataset)
summary(cov_only)

# fit a first unadjusted model to look at the association between the mixture and the outcome
# TELOMEAN = Mean Telomere Length
## there's a signal only in the negative direction
#summary(dataset[,mixture])
results1 = gwqs(log_TELOMEAN ~ wqs, mix_name = mixture, data = dataset, q = 10, validation = 0.6, 
                b = 100, b1_pos = FALSE, b1_constr = FALSE, family = "gaussian", seed = 123)
gwqs_barplot(results1)
gwqs_scatterplot(results1)
gwqs_fitted_vs_resid(results1)
gwqs_summary_tab(results1)
sum(results1$bres$b1 < 0) 
## all 100 bootstraps had negative beta estmates in the GLM

summary(results1$fit)
results1$final_weights

# adjusting for covariates:
# blood data: LBXWBCSI LBXLYPCT LBXMOPCT LBXEOPCT LBXBAPCT LBXNEPCT
# demographics: age_cent age_sq race_cat bmi_cat3 ln_lbxcot edu_cat male
dataset$race_cat=as.factor(dataset$race_cat)
dataset$bmi_cat3=as.factor(dataset$bmi_cat3)
dataset$edu_cat=as.factor(dataset$edu_cat)
is.factor(dataset$race_cat)
is.factor(dataset$bmi_cat3)
is.factor(dataset$edu_cat)
is.factor(dataset$male)

# positive direction -  when adjusted the association is in the positive direction: F08 F03 PCB
result2 = gwqs(log_TELOMEAN ~ wqs + 
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = FALSE, family = "gaussian", 
               seed = 123)  ## signal="abst","t2"(default),"expt"
## can require stricter convergence:  control = list(trace = FALSE, maxit = 2000, reltol = 1e-10) 
## different optimization algorithms, optim.method="BFGS"
## plan_strategy = "multisession"
gwqs_barplot(result2)
gwqs_summary_tab(result2)
gwqs_scatterplot(result2) 
sum(result2$bres$b1 > 0)   ## 95 of 100 are positive beta estimates
sum(result2$bres$b1 < 0)

# negative direction - although constrained to be negative, only 23 of 100 bootstraps resulted in a negative coeff; 
#                      LBX074 099 194 180 170 and final estimate is positive and sig
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
sum(result3$bres$b1 > 0)

#compare to random subset WQS - F08 F03 PCB D05 D07 pos and sig
result4 = gwqs(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123,rs=TRUE, n_var=5) #3 for set of 9 and 5 for set of 18

gwqs_summary_tab(result4) 
gwqs_barplot(result4)
gwqs_scatterplot(result4)
gwqs_weights_tab(result4)
sum(result4$bres$b1 > 0)

# WQS with interaction - b12 NS (p=0.378)
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
## females(0): F08 F03 PCB 194 180 138 D05
## males(1): F03 F08 D07 PCB (~D05)
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
## the interaction is NS (p=0.851)
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
result5rh = gwqsrh(log_TELOMEAN ~ wqs*male +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + sex, 
                 stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123, rh=100,
               plan_strategy = "multisession")
gwqs_summary_tab(result5rh)
gwqsrh_boxplot(result5rh)
gwqs_weights_tab(result5rh)
summary(result5rh$fit)
result5rh$final_weights


write.csv(result5rh$wmat, paste0(directory_path_out, "labstWQSwtsrh.csv"))
