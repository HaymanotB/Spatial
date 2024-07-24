install.packages("emdi")
library(emdi)
data("eusilcA_smp")
View(eusilcA_smp)
data("eusilcA_pop")
table(eusilcA_smp$state)
table(eusilcA_pop$state)
summary(as.numeric(table(eusilcA_smp$district)))
mdi_direct <- direct(y = "eqIncome", smp_data = eusilcA_smp,
                      smp_domains = "district", weights = "weight", threshold = 10885.33,
                      var = TRUE)
mdi_direct
summary(mdi_direct)
#######
The R object emdi_direct is of classes emdi and direct.
An example of using the ebp method for computing point and MSE estimates for the predened
indicators and two custom indicators, namely the minimum and maximum equivalized
income is provided below
mdi_model <- ebp(fixed = eqIncome ~ gender + eqsize + cash + self_empl +
                   + unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
                   + fam_allow + house_allow + cap_inv + tax_adj, pop_data = eusilcA_pop,
                 + pop_domains = "district", smp_data = eusilcA_smp,
                 + smp_domains = "district", threshold = 10885.33, MSE = TRUE,
                 + custom_indicator = list(my_max = function(y, threshold){max(y)},
                                           +
                                             my_min = function(y, threshold){min(y)}))