# ----------------------------------------------------------------------------- #
# Relicate Vision Zero example from "A (Flexible) Synthetic Control Method ..." #
# ----------------------------------------------------------------------------- #

# Get helper functions (change ... to relevant folder)

source(".../CSCM_helper_functions.R")

# Load/install required packages

packages <- c("ggplot2", "Synth", "glmnet",
              "dplyr", "osqp", "optimx")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
sapply(packages, require, character.only=TRUE)

# Version check

ifelse(.Machine$sizeof.pointer==4,warning("It appears you are running a 32-bit version of R. Please switch to 64-bit to ensure that the optimization does not produce different results due to rounding errors."),
       print("64-bit version check: OK."))

# Load the data (change ... to relevant folder)

viszero <- read.csv(".../viszero.csv", header=TRUE)

# --------------------- #
# Run the main analysis #
# --------------------- #

road.main.res <- countSynth(data=viszero,
                            predictors=NULL, #Auxiliary covariates can be included here as well (passed to dataprep synth, averaged across pre-period with na.rm=TRUE). Caution: needs to have non-missing values in all holdout samples.
                            dependent=c("deathrate_mln"),
                            unit.variable="ID",
                            time.variable = "TIME",
                            treatment.identifier = 25, # 25 = Sweden
                            controls.identifier = c(2,3,5,7,9,10,13,14,16), #We exclude UK, Norway, Netherlands due to similar policies. Est, Rou, Lva due to data issues. Other countries due to missing values.
                            t_int=27, #Treatment time point (1997 in real years)
                            min_1se=F, #Use 1se rule for lambda CV? (try changing to T for sensitivity analysis)
                            K=2) #Number of holdout periods for cross-fitting (try changing to 3 for sensitivity analysis)

# Print the cross-fitted ATT estimate (with 95% CIs)

road.main.res$ATT

# -------------------------- #
# Produce unit weights table #
# -------------------------- #

synth.w <- round(road.main.res$unit.weight.SCM, 2)
csynth.w <- round(road.main.res$unit.weight.full.sample, 2)
cunit.names <- unique(viszero %>% dplyr::filter(ID %in% c(2,3,5,7,9,10,13,14,16)) %>% dplyr::select(ID, COUNTRYNAME))
cunit.ordered <- as.data.frame(cunit.names[order(cunit.names$ID),])
cunit.synth = cunit.csynth = cunit.ordered
cunit.synth$Country <- c("Belgium", "Bulgaria", "Germany", "Spain", "Finland","France","Croatia","Hungary","Italy")
cunit.synth$SCM.W <- synth.w
cunit.synth$CSCM.W <- csynth.w

cunit.synth[,3:5]


# -------------------- #
# Counterfactuals plot #
# -------------------- #

synth.orig <- road.main.res$dataprep.main$Y0plot %*% as.numeric(road.main.res$unit.weight.SCM)
cscm.cf <- road.main.res$dataprep.main$Y0plot %*% as.numeric(road.main.res$unit.weight.full.sample)
year <- 1970:2015

adf <- as.data.frame(cbind(year,synth.orig,rep("SCM",length(year))))
bdf <- as.data.frame(cbind(year,cscm.cf,rep("CSCM",length(year))))
cdf <- as.data.frame(cbind(year,road.main.res$dataprep.main$Y1plot,rep("Observed",length(year))))
colnames(adf)=colnames(bdf)=colnames(cdf)=c("Year","Y","Method")
cfdf <- rbind(adf,bdf)
cfdf <- cfdf %>% mutate(Y = as.numeric(as.character(Y)),
                        Year = as.numeric(as.character(Year)))
cdf <- cdf %>% mutate(Y = as.numeric(as.character(Y)),
                      Year = as.numeric(as.character(Year)))

g_m <- ggplot(cfdf, aes(y=Y, x=Year, group=Method, linetype=Method, color=Method)) + geom_vline(xintercept=1997) + geom_point(y=rep(cdf$Y,2),x=rep(cdf$Year,2),color="black",size=2) + scale_color_manual(values=c("black","gray60")) + theme_bw() + geom_line(size=1.2)
g_p <- g_m + theme_bw() +
  theme(text = element_text(size=15)) +
  ylab("Fatality rate") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom")

g_p

# --------------------------- #
# Effect estimate using synth #
# --------------------------- #

RR.Synth.raw = mean(road.main.res$dataprep.main$Y1plot[27:length(synth.orig)])/mean(synth.orig[27:length(synth.orig)])
RR.CSCM.raw = mean(road.main.res$dataprep.main$Y1plot[27:length(synth.orig)])/mean(cscm.cf[27:length(synth.orig)])

RR.Synth.raw #RR using synth counterfactual
RR.CSCM.raw #RR using (not debiased) CSCM counterfactual
