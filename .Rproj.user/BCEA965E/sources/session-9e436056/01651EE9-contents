#################################################################################
##         PROGRAM: ZA_sae.R
##
##         PROJECT: Small area estimation of district-level indicators
##                  for South Africa
##
## INVESTIGATOR(S): Ehimario Igumbor, Steve Gutreuter, Njiri Wabiri,
##                  Khangelani Zuma, Mitesh Desai and Lizette Durand
##
##    REFERENCE(S): Fay RE II and Herriot RA.  Estimates of income
##                     for small places: An application of James-Stein
##                     procedures to census data.  Journal of the
##                     American Statistical Association. 1979; 74(366):269-277.
##                  Marhuenda Y., Molina I., Morales D. Small area estimation
##                     with spatio-temporal Fay-Herriot models Computational
##                     Statistics and Data Analysis. 2013; 58:308-325.
##                  Molina I, Marhuenda Y.  sae: An R package for small area
##                     estimation. The R Journal. 2015; 7(1):81-98.
##
##           INPUT: South_Africa_District_HIV_Data.xlsx
##
## REQUIRED R PKGS: sae, readxl, ggplot2, reshape, rgdal, scales, car
##
##       IMPORTANT: The working directory and data directory MUST be
##                  set as appropriate for your file structure. That
##                  is done in the "IMPORTANT" block below.
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/HIDMSB Statistics, Estimation
##                  and Modeling Team
##                  E-mail: sgutreuter@cdc.gov
##
##      DISCLAIMER: Although this program has been used by the Centers
##                  for Disease Control & Prevention (CDC), no warranty,
##                  expressed or implied, is made by the CDC or the U.S.
##                  Government as to the accuracy and functioning of the
##                  program and related program material nor shall the
##                  fact of distribution constitute any such warranty,
##                  and no responsibility is assumed by the CDC in
##                  connection therewith.
##
##            DATE: 2015-06-05
##
#################################################################################

#################################################################################
## IMPORTANT: Define path to the base location of project files here
##
## Example:
basepath <- "F:/Spatial"
#################################################################################
## No additional changes should be needed below this line
#################################################################################

workpath <- file.path(basepath, "code")
datapath <- file.path(basepath, "data")
setwd(workpath)

library(sae)
library(readxl)
library(ggplot2)
library(reshape)
library(rgdal)
library(sf) # instead of rgdal
library(scales)
library(car)

#################################################################################
## Create file structure as needed and test for presence of data
#################################################################################
if(!file.exists(file.path(workpath, "output"))){
    dir.create(file.path(workpath, "output"))}
outdir <- file.path(workpath, "output")
datafile <- file.path(datapath, "South_Africa_District_HIV_Data.xlsx")
if(!file.exists(datafile)) stop("Data file not found")

#########################################################################
## Define inverse-logit (expit) functions
#########################################################################
expit <- function (x){1/(1 + exp(-x))}

#########################################################################
## Import the data from the spreadsheet as a simple data frame
## NOTE: The tibbles produced by read_xlsx break mseFH.
#########################################################################
ZAdata <- data.frame(read_xlsx(path = datafile,
                               sheet = "data",
                               col_types = c(rep("text", 5),
                               rep("numeric", 18))))

#########################################################################
## Create a variable containing the survey variance estimates and
## the mean design effect from the survey.
#########################################################################
ZAdata$varHHS <- ZAdata$se_HHS^2
DEFF <- ZAdata$deft_HHS^2                       # DEFF from Kish Deft
DEFF[DEFF < 1] <- 1                             # Truncate implausibles
mean.deff <- mean(DEFF, na.rm = TRUE)
DEFF[is.na(DEFF)] <- mean.deff

#########################################################################
## Compute HHS missing values for Overberg (DC3) and Central Karoo (DC5).
## Those values were omitted from the HHS report based upon the concern
## that precision was too low, and therefore the estimates were not
## trusted.  The missing prevalences are computed as
## 100*HIVpos_HHS/Ntested_HHS.  The missing variances are approximated
## as the sample variance under SRS times the mean Deff over all
## districts.
#########################################################################
for(jj in c("DC3", "DC5")){
    ZAdata$p_HHS[ZAdata$code == jj] <- {100 *
        ZAdata$HIVpos_HHS[ZAdata$code == jj] /
            ZAdata$Ntested_HHS[ZAdata$code == jj] }
    ZAdata$varHHS[ZAdata$code == jj] <- 10000 *
        (ZAdata$HIVpos_HHS[ZAdata$code == jj]/100 * (1 -
                           ZAdata$HIVpos_HHS[ZAdata$code == jj]/100) /
                           ZAdata$Ntested_HHS[ZAdata$code == jj] *
                           mean.deff)
    ZAdata$se_HHS[ZAdata$code == jj] <-
        sqrt(ZAdata$varHHS[ZAdata$code == jj])
}
rm(jj, mean.deff)
#########################################################################
## Transform prevalence percentages and their variances to the logit
## scale for proportions.  The variances are delta-method approximations.
#########################################################################
ZAdata$eta_HHS <- logit(ZAdata$p_HHS/100)
p2 <- (ZAdata$p_HHS/100)^2
ZAdata$var_eta_HHS <- 1e-04*ZAdata$varHHS/(p2 * (1 - p2))
ZAdata$eta_ANC <- logit(ZAdata$p_anc2012/100)
#########################################################################
## Read the spatial adjacency matrix.  The sheet "proxmat" was created
## manually from a district map of South Africa.  Rows identify districts
## (labeled by district codes) and the adjacent neighbors of those
## districts.  Column names are district codes.  The non-zero elements in
## each row equal 1/n, where n is the number of adjacent neighboring
## districts of the row district.  Rows sum to 1.
#########################################################################
ZAprox <- read_xlsx(path = datafile,
                    sheet = "proxmat",
                    col_types = c(rep("text", 3),
                                  rep("numeric", 52)))
ZAprox <- as.matrix(ZAprox[, 4:55], nrow = 52)

#########################################################################
#########################################################################
## Frequentist small-area estimation for 13 model variations
#########################################################################
#########################################################################
## FH1: Auxilliary predictor is logit(p_anc2012).
FH1 <- mseFH(eta_HHS ~ eta_ANC,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH2: Auxilliary predictors are logit(p_anc2012), p_anc2012 and fdwell
FH2 <- mseFH(eta_HHS ~ eta_ANC + fdwell,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH3: Auxilliary predictors are logit(p_anc2012), and dep_ratio
FH3 <- mseFH(eta_HHS ~ eta_ANC + dep_ratio,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH4: Auxilliary predictors are logit(p_anc2012), and se_qtile
FH4 <- mseFH(eta_HHS ~ eta_ANC + as.factor(se_qtile),
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH5: Auxilliary predictors are logit(p_anc2012), and mmr
FH5 <- mseFH(eta_HHS ~ eta_ANC + mmr,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH6: Auxilliary predictors are logit(p_anc2012), fdwell and dep_ratio
FH6 <- mseFH(eta_HHS ~ eta_ANC + fdwell + dep_ratio,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH7: Auxilliary predictors are logit(p_anc2012), fdwell dep_ratio and
##      mmr
FH7 <- mseFH(eta_HHS ~ eta_ANC + fdwell + dep_ratio + mmr,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH8: Auxilliary predictors are logit(p_anc2012), fdwell dep_ratio,
##      mmr and fac_del
FH8 <- mseFH(eta_HHS ~ eta_ANC + fdwell + dep_ratio + mmr + fac_del,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## FH9: Auxilliary predictors are  fdwell dep_ratio,  mmr and fac_del
FH9 <- mseFH(eta_HHS ~ fdwell + dep_ratio + mmr + fac_del,
             vardir = var_eta_HHS, method = "ML", data = ZAdata)
## SFH1: Spatial Fay-Herriot model.  Auxilliary predictor is
##       logit(p_anc2012)
SFH1 <- mseSFH(eta_HHS ~ eta_ANC,
               vardir = var_eta_HHS, method = "ML", proxmat = ZAprox,
               data = ZAdata)
## SFH2: Spatial Fay-Herriot model.  Auxilliary predictors are
##       logit(p_anc2012) and dep_ratio.
SFH2 <- mseSFH(eta_HHS ~ eta_ANC + dep_ratio,
               vardir = var_eta_HHS, method = "ML",
               proxmat = ZAprox, data = ZAdata)
## SFH3: Spatial Fay-Herriot model.  Auxilliary predictors are
##       fdwell dep_ratio,  mmr and fac_del.
SFH3 <- mseSFH(eta_HHS ~ fdwell + dep_ratio + mmr + fac_del,
               vardir = var_eta_HHS, method = "ML",
               proxmat = ZAprox, data = ZAdata)
## SFH4: Spatial Fay-Herriot model without non-spatial auxilliary
##       predictors.
SFH4 <- mseSFH(eta_HHS ~ 1,
               vardir = var_eta_HHS, method = "ML",
               proxmat = ZAprox, data = ZAdata)
## Assemble a table of information criteria
models <- c("FH1", "FH2", "FH3", "FH4", "FH5", "FH6", "FH7", "FH8",
            "FH9", "SFH1", "SFH2", "SFH3", "SFH4")
ICtable <- data.frame(NULL)
for(jj in models){
    fit_ <- eval(parse(text = paste(jj, "$est$fit$goodness", sep = "")))
    names(fit_) <- NULL
    tmpIC <- data.frame(model = jj, fit_[1], fit_[2], fit_[3] )
    names(tmpIC) <- c("Model", "Log-likelihood", "AIC", "BIC")
    ICtable <- rbind(ICtable, tmpIC)
}
ICtable$Predictors <- c("eta_ANC",
                        "eta_ANC + fdwell",
                        "eta_ANC + dep_ratio",
                        "eta_ANC + se_qtile",
                        "eta_ANC + mmr",
                        "eta_ANC + fdwell + dep_ratio",
                        "eta_ANC + fdwell + dep_ratio + mmr",
                        "eta_ANC + fdwell + dep_ratio + mmr + fac_del",
                        "fdwell + dep_ratio + mmr + fac_del",
                        "eta_ANC + SAR",
                        "eta_ANC + dep_ratio + SAR",
                        "fdwell + dep_ratio + mmr + fac_del + SAR",
                        "SAR")
print(ICtable)
## Compute Wald approximate confidence intervals for naive (design-based)
## and Fay-Herriot estimates of district level prevalence and numbers of
## people living with HIV (PLHIV) from the AIC-best model (FH3).
p_hat_FH <- 100*expit(FH3$est$eblup)
lclFH <- 100*expit(FH3$est$eblup - 1.96*sqrt(FH3$mse))
uclFH <- 100*expit(FH3$est$eblup + 1.96*sqrt(FH3$mse))
lclRAW <- ZAdata$p_HHS - 1.96*ZAdata$se_HHS
lclRAW[which(lclRAW < 0)] <- 0
uclRAW <- ZAdata$p_HHS + 1.96*ZAdata$se_HHS
plhiv_HHS <- round(ZAdata$pop2014 * ZAdata$p_HHS/100)
lcl.plhiv_HHS <-
    round(plhiv_HHS - 1.96*ZAdata$pop2014*ZAdata$se_HHS/100)
lcl.plhiv_HHS[which(lcl.plhiv_HHS < 0)] <- 0
ucl.plhiv_HHS <-
    round(plhiv_HHS + 1.96*ZAdata$pop2014*ZAdata$se_HHS/100)
plhiv_sae <- round(ZAdata$pop2014 * p_hat_FH/100)
lcl.plhiv_sae <- round(ZAdata$pop2014*lclFH/100)
ucl.plhiv_sae <- round(ZAdata$pop2014*uclFH/100)
vp_sae <- (exp(FH3$est$eblup)/((1 + exp(FH3$est$eblup))^2))^2 * FH3$mse
rse_sae <- round(100*sqrt(vp_sae) / expit(FH3$est$eblup), digits = 1)
rse_HHS <- round(100*ZAdata$se_HHS / ZAdata$p_HHS, digits = 1)
rse <- data.frame(ZAdata$district, ZAdata$code, rse_sae, rse_HHS)
names(rse) <- c("District", "Code", "RSE_sae", "RSE_HHS")

#########################################################################
## Tabulate and plot RSEs and HIV prevalence by district
#########################################################################
theme_set(theme_bw())
## Plot RSEs
rseplot <- ggplot(rse, aes(x = RSE_HHS, y = RSE_sae)) + geom_point()
rseplot <- rseplot + scale_x_continuous(limits = c(0, 70)) +
          scale_y_continuous(limits = c(0, 70))
rseplot <- rseplot + xlab("Relative MSE direct (%)") + ylab("Relative MSE Fay-Herriot (%)")
rseplot <- rseplot + geom_abline(intercept = 0, slope = 1)
plot(rseplot)
## Construct prevalence data frame
est.tmp <- data.frame(ZAdata$district, ZAdata$code, ZAdata$p_HHS, lclRAW, uclRAW,
                      p_hat_FH, lclFH, uclFH)
names(est.tmp) <- c("District", "Code", "p_HHS", "lcl_HHS", "ucl_HHS", "p_sae",
                    "lcl_sae", "ucl_sae")
## Sort in order of increasing p_sae:
est.tmp <- est.tmp[order(est.tmp$p_sae), ]
est.tmp$District <- factor(as.character(est.tmp$District),
                           levels = est.tmp$District,
                           ordered = TRUE)
## Reshape est.tmp to enable helpful graphics:
mlt.tmp <- melt(est.tmp, id = c("District", "Code"))
mlt.tmp <- cbind(mlt.tmp, colsplit(mlt.tmp$variable, "_",
                                   names = c("stat", "Source")))
est_prev <- cast(mlt.tmp, District + Code + Source ~ stat)
est_prev$Estimate <- factor(est_prev$Source,
                            labels = c("Survey", "Fay-Heriott"))
est_prev$rownum <- as.numeric(rownames(est_prev))
## Create plots of  point estimates and confidence intervals by district
theme_set(theme_bw())
## Plot prevalence
prevplot <- ggplot(est_prev, aes(x = District, y = p,
                                 group = Estimate,
                                 colour = Estimate,)) +
                   coord_flip() +  ylab("HIV Prevalence (%)")
prevplot <- prevplot + scale_color_manual(values =
                                          c("#D55E00", "#009E73"))
prevplot <- prevplot + geom_point(size = 2)
prevplot <- prevplot + geom_errorbar(aes(ymin = lcl,
                                         ymax = ucl,
                                         group = Estimate))
prevplot <- prevplot + theme(legend.position="top")
plot(prevplot)
#########################################################################
## Tablulate and plot numbers (100,000s) of PLHIV by district
#########################################################################
est.tmp <- data.frame(ZAdata$district, ZAdata$code, plhiv_HHS, lcl.plhiv_HHS,
                      ucl.plhiv_HHS, plhiv_sae, lcl.plhiv_sae, ucl.plhiv_sae)
names(est.tmp) <- c("District", "Code",
                    "PLHIV_HHS", "lcl.PLHIV_HHS", "ucl.PLHIV_HHS",
                    "PLHIV_sae", "lcl.PLHIV_sae", "ucl.PLHIV_sae")
## Sort in order of increasing PLHIV:
est.tmp <- est.tmp[order(est.tmp$PLHIV_sae), ]
est.tmp$District <- factor(as.character(est.tmp$District),
                           levels = est.tmp$District,
                           ordered = TRUE)
## Reshape est.tmp dataframe to enable helpful graphics:
mlt.tmp <- melt(est.tmp, id = c("District", "Code"))
mlt.tmp <- cbind(mlt.tmp, colsplit(mlt.tmp$variable, "_",
                                   names = c("stat", "Source")))
est_plhiv <- cast(mlt.tmp, District + Code +
                  Source ~ stat)
est_plhiv$Estimate <- factor(est_plhiv$Source, labels =
                             c("Survey", "Fay-Heriott"))
est_plhiv$rownum <- as.numeric(rownames(est_plhiv))
## Create plots of numbers of PLHIV by district
PLHIVplot <- ggplot(est_plhiv, aes(x = District, y = PLHIV/100000,
                                   group = Estimate,
                                   colour = Estimate)) +
                    coord_flip() +  ylab("PLHIV (hundred-thousands)")
PLHIVplot <- PLHIVplot + scale_color_manual(values =
                                            c("#D55E00", "#009E73"))
PLHIVplot <- PLHIVplot + geom_point(size = 2)
PLHIVplot <- PLHIVplot + geom_errorbar(aes(ymin = lcl.PLHIV/100000,
                                           ymax = ucl.PLHIV/100000,
                                           group = Estimate))
PLHIVplot <- PLHIVplot + theme(legend.position="top")
plot(PLHIVplot)

## Write all estimates to a csv file
est_prev <- est_prev[order(est_prev$Code, est_prev$Estimate) ,]
est_plhiv <- est_plhiv[order(est_plhiv$Code, est_plhiv$Estimate) ,]
all.est <- data.frame(subset(est_prev,
                             select = c(District, Code, Estimate,
                                 p, lcl, ucl)),
                      subset(est_plhiv,
                             select = c(PLHIV, lcl.PLHIV, ucl.PLHIV)))
names(all.est) <- c("District", "Code", "Estimate", "Prevalence_%",
                    "Prevalence_LCL", "Prevalence_UCL",
                    "NoPLHIV", "NoPLHIV_LCL", "NoPLHIV_UCL")
write.csv(all.est, file = file.path(workpath, "output", "estimates.csv"),
          row.names = FALSE)


#########################################################################
## Extra plotting
#########################################################################
beta <- matrix(FH3$est$fit$estcoef$beta, ncol = 1)
mm <- model.matrix(eta_HHS ~ eta_ANC + dep_ratio, data = ZAdata)
pred <- mm %*% beta
resid1 <- ZAdata$eta_HHS - pred
resid2 <- FH3$est$eblup - pred
png(filename = file.path(workpath, "Residual_plot_domain_vs_synthetic.png"))
plot(pred, resid1, pch = 20, ylab = "Residual", xlab = "Predicted",
     ylim = c(-1, 1))
abline(a = 0, b = 0)
dev.off()
png(filename = file.path(workpath, "Residual_plot_EBLUP_vs_synthetic.png"))
plot(pred, resid2, pch = 20, ylab = "Residual", xlab = "Predicted",
     ylim = c(-1, 1))
abline(a = 0, b = 0)
dev.off()


#########################################################################
## Clean up after R.
#########################################################################
rm(all.est, DEFF, est_plhiv, est_prev, p_hat_FH, p2, lcl.plhiv_HHS,
   lcl.plhiv_sae, maploc, mlt.tmp, p_hat_FH, plhiv_sae, plot_data,
   ucl.plhiv_HHS, ucl.plhiv_sae, varSRS, ZAdists, ZAprox, lclFH, lclRAW,
   mapfiles, plot_data, rse, uclFH, uclRAW)
