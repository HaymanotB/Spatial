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
#The R object emdi_direct is of classes emdi and direct.
#An example of using the ebp method for computing point and MSE estimates for the predened
#indicators and two custom indicators, namely the minimum and maximum equivalized
#income is provided below
 
mdi_model <- ebp(fixed = eqIncome ~ gender + eqsize + cash + self_empl +
                    unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent +
                    fam_allow + house_allow + cap_inv + tax_adj, pop_data = eusilcA_pop,
                 pop_domains = "district", smp_data = eusilcA_smp,
                  smp_domains = "district", threshold = 10885.33, MSE = TRUE,
                  custom_indicator = list(my_max = function(y, threshold){max(y)},
                                           
                                            my_min = function(y, threshold){min(y)}))

mdi_model
summary(mdi_model)
#The plots functionused for residual analyses of the object emdi_model

plot(mdi_model, label = "no_title", color = c("red3", "red4"))

#estimate the median of equalized income and the Gini coefficient and
# the corresponding CV estimates for the first 6 districts in Austria.

head(estimators(mdi_model, indicator = c("Gini", "Median"),
                    MSE = FALSE, CV = TRUE))
compare_plot(mdi_direct, mdi_model, indicator = c("Gini", "Median"),
             label = "no_title", color = c("red3", "blue"))

#maps of point estimates and CVs of the median income for 94 districts in Austria.
map_plot(mdi_model, MSE = FALSE, CV = TRUE, map_obj = shape_austria_dis,
          indicator = "Median", map_dom_id = "PB")
# Export of the summary output and estimates to Exce
write.excel(mdi_model, file = "excel_output.xlsx", indicator = "Median",
             MSE = FALSE, CV = TRUE)
write.csv(eusilcA_pop,"eusilcA_pop.CSV")
