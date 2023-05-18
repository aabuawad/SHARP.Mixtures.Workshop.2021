# Workshop on mixtures, Columbia August 2021

#rm(list = ls())
# install the packages needed for the workshop
install.packages("gWQS")
# import the packages just installed
library(gWQS)

# define the path  
#directory_path_out = "/Users/gennic01/Desktop/RESEARCH/Columbia Mixtures Workshop/"
# import the dataset
dataset = read_csv(here::here('Data/studypop.csv'))
dim(dataset)
dataset = dataset[complete.cases(dataset),]
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
               seed = 123)  
## signal="abst","t2"(default),"expt"
## can require stricter convergence:  control = list(trace = FALSE, maxit = 2000, reltol = 1e-10) 
## different optimization algorithms, optim.method="BFGS"
## plan_strategy = "multisession"
gwqs_barplot(result2)
gwqs_summary_tab(result2)
gwqs_scatterplot(result2)
sum(result2$bres$b1 > 0)   ## 95 of 100 are positive beta estimates
sum(result2$bres$b1 < 0)

# negative direction - although constrained to be negative, only 23 of 100 bootstraps resulted in a negative coeff; 
#                      final estimate is positive and sig
result3 = gwqs(log_TELOMEAN ~ wqs+LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + male, 
               mix_name = mixture, 
               data = dataset, q = 10, 
               validation = 0.6,  b = 100, b1_pos = FALSE, b1_constr = TRUE, family = "gaussian", 
               seed = 123)
gwqs_summary_tab(result3)
gwqs_barplot(result3)
gwqs_scatterplot(result3)
sum(result3$bres$b1 < 0)  #23 of 100 are negative
sum(result3$bres$b1 > 0)

#compare to random subset WQS - F08 F03 PCB D05 PCB D07 pos and sig
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


#####################################################################
## divergent plots
## load libraries
library(gapminder)
library(tidyverse)
library(stringr)

forplot = result5rh$final_weights  ## specify the data to be used in the plot; i.e., WQS str or WQS str int model

forplot$strata = str_sub(forplot$mix_name,-1)
forplot$chem = str_sub(forplot$mix_name,1,nchar(as.character(forplot$mix_name))-2)
forplot$strataname = ifelse(forplot$strata == "1", "Male", "Female")

resultstr <- forplot %>%
  mutate(mean_weight = ifelse(strata == "1",
                              mean_weight,
                              -1*mean_weight)) 
resultstr

## calculate breaks values
breaks_values <- pretty(resultstr$mean_weight)

ggplot(data=resultstr,
       aes(reorder(x = chem, mean_weight), y = mean_weight, fill = strataname))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values))+
  theme_minimal()+
  geom_hline(yintercept = 1/nrow(resultstr), linetype='dotted', col='black')+  
  geom_hline(yintercept = -1/nrow(resultstr), linetype='dotted', col='black')+
  geom_hline(yintercept = 0, col='black')+  
  labs(x="Component", y="Mean Weight", fill="Strata")

wt_perc_males = round(100*sum(resultstr[resultstr$strata == 1,]$mean_weight),1)
wt_perc_males

################################ takes FOREVER ##########################
#repeated holdout WQS
# positive direction
result5rh = gwqsrh(log_TELOMEAN ~ wqs +
                 LBXWBCSI + LBXLYPCT + LBXMOPCT + LBXEOPCT + LBXBAPCT + LBXNEPCT + age_cent + age_sq + 
                 race_cat + bmi_cat3 + ln_lbxcot + edu_cat + sex, 
                 stratified="sex", mix_name = mixture, data = dataset, q = 10, 
               validation = 0.6, b = 100, b1_pos = TRUE, b1_constr = TRUE, family = "gaussian", 
               seed = 123, rh=30,
               plan_strategy = "multisession")
gwqs_summary_tab(result5rh)
gwqsrh_boxplot(result5rh)
gwqs_weights_tab(result5rh)
summary(result5rh$fit)
result5rh$final_weights
result5rh$final_weights$mean_weight = result5rh$final_weights$Estimate


write.csv(result5rh$wmat, paste0(directory_path_out, "labstWQSwtsrh30.csv"))


#### boxplots ################################################################
## for WQS st rh model 
boxplot(result5rh$coefmat[,"wqs"],
       xlab="WQS",  ylab = "Beta Estimates, rh=30", 
       col = c("orange"),  border = "brown", horizontal = FALSE,  notch = TRUE)


## for WQS st int rh model
WQS = result6rh$coefmat[,"wqs"]
WQS_SEX = result6rh$coefmat[,"wqs:sex1"]

boxplot(WQS, WQS_SEX,
        at = c(1,2),
        names = c("WQS", "WQS_SEX"),
        ylab = "Beta Estimates, rh=30",
        col = c("orange","red"),
        border = "brown",
        horizontal = FALSE,
        notch = TRUE
)
