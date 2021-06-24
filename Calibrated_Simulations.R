# -------------------------------------------------------------------------------- #
# Relicate calibrated simulations from "A (Flexible) Synthetic Control Method ..." #
# -------------------------------------------------------------------------------- #

# Get helper functions (change ... to relevant folder)

source(".../CSCM_helper_functions.R")

# Load/install required packages

packages <- c("ggplot2", "Synth", "glmnet",
              "dplyr", "osqp", "optimx","MASS","cowplot", "gsynth")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
sapply(packages, require, character.only=TRUE)

# Set seed and random number generation settings for reproducibility of the simulation studies.
# Changes were made to random number generator in R in version 3.6.0. 
# Make sure to update R to reproduce results if you have and older version than that.

ifelse(.Machine$sizeof.pointer==4,warning("It appears you are running a 32-bit version of R. Please switch to 64-bit to ensure that the optimization does not produce different results due to rounding errors."),
       print("64-bit version check: OK."))
RNGversion("3.6.0")
set.seed(131348)

# Load the data (change ... to relevant folder)

viszero <- read.csv(".../viszero.csv", header=TRUE)
sygdat <- read.csv(".../syg_data.csv") %>% dplyr::select(-X) #Change ... to the folder that contains the data file.
sygdat <- sygdat %>% mutate(State=as.character(State))

## Prepare data for simulation

# Vision zero data
simcaldf <- viszero %>% filter(ID %in% c(25,2,3,5,7,9,10,13,14,16))
simcaldf <- simcaldf %>% mutate(D = ifelse(TIME>1996 & ID==25, 1, 0), logY = log(deathrate_mln))
testsim <- factor.init(simcaldf,t_int=27,id="ID",time="TIME",outcome="logY",trvar="D",r=3)

# SYG data
simcaldf_long <- sygdat %>% filter(State.Code %in% c(12,5,9,10,15,19,23,24,25,31,34,36,38,39,44,56))
simcaldf_long <- simcaldf_long  %>% mutate(D = ifelse(time>80 & State.Code==12, 1, 0), logY = log(Homicide_count+1))
testsim_long <- factor.init(simcaldf_long,t_int=81,id="State.Code",time="time",outcome="logY",trvar="D",r=3)


# Note: If you want to run simulations using your own dataset, replace
#       the data in the factor.init function(s) above with a balanced panel
#       dataset (make sure you log-transform the outcomes before doing so
#       since gsynth assumes a linear model).

#       t_int is the intervention time point (relative to start of period)
#       trvar is a treatment*post interaction variable
#       id is the id variable
#       time is time variable

# Number of simulation replicates
nreps <- 1000

# ---------- Run a calibrated simulation study with empirical residual covariances ------------- #

sim_att_list_2 = list()
sim_est_list_2 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim)
  simres  <- countSynth(data=simdats[["simdata"]],
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=27,
                        min_1se=F,
                        K=2,
                        full.model=FALSE)
  sim_att_list_2[[i]] <- simres$ATT
  sim_est_list_2[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
simres_att_2 <- do.call("rbind",sim_att_list_2)
simres_est_2 <- do.call("rbind",sim_est_list_2)

coverage_2 <- simres_att_2 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                        ifelse(RR.upper>1 & RR.lower>1,0,1)))
# K3

sim_att_list_3 = list()
sim_est_list_3 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim)
  simres  <- countSynth(data=simdats[["simdata"]],
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=27,
                        min_1se=F,
                        K=3,
                        full.model=FALSE)
  sim_att_list_3[[i]] <- simres$ATT
  sim_est_list_3[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
simres_att_3 <- do.call("rbind",sim_att_list_3)
simres_est_3 <- do.call("rbind",sim_est_list_3)

coverage_3 <- simres_att_3 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                        ifelse(RR.upper>1 & RR.lower>1,0,1)))


# Plot histograms of effects

library(ggplot2)

coverplotdf <- as.data.frame(cbind(simres_att_2, simres_att_3))
colnames(coverplotdf) <- c("RR_2", "RR_l2", "RR_u2", "RR_3", "RR_l3", "RR_u3")
bw1 <-  with(coverplotdf, 2 * IQR(log(RR_2)) / length(log(RR_2))^(1/3))
hist1 <- ggplot(coverplotdf, aes(x=log(RR_2))) + geom_histogram(aes(y=..density..),binwidth = bw1,
                                                                fill="gray70",
                                                                color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(coverplotdf$RR_2))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(coverplotdf$RR_2)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("Empirical residuals, K=2")

bw2 <-  with(coverplotdf, 2 * IQR(log(RR_3)) / length(log(RR_3))^(1/3))
hist2 <- ggplot(coverplotdf, aes(x=log(RR_3))) + geom_histogram(aes(y=..density..),binwidth = bw2,
                                                                fill="gray70",
                                                                color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(coverplotdf$RR_3))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(coverplotdf$RR_3)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("Empirical residuals, K=3")


# ----------- Repeat the simulation simulating from an AR(3) error distribution ----------- #

# Run 1000 simulations for K=2 and then K=3

ar3_att_list_2 = list()
ar3_est_list_2 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim, ar=c(0.4,0.2,0.1))
  simres  <- countSynth(data=simdats[["simdata"]],
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=27,
                        min_1se=F,
                        K=2,
                        full.model=FALSE)
  ar3_att_list_2[[i]] <- simres$ATT
  ar3_est_list_2[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
ar3res_att_2 <- do.call("rbind",ar3_att_list_2)
ar3res_est_2 <- do.call("rbind",ar3_est_list_2)

ar3_coverage_2 <- ar3res_att_2 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                        ifelse(RR.upper>1 & RR.lower>1,0,1)))
ar3_att_list_3 = list()
ar3_est_list_3 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim)
  simres  <- countSynth(data=simdats[["simdata"]],
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=27,
                        min_1se=F,
                        K=3,
                        full.model=FALSE)
  ar3_att_list_3[[i]] <- simres$ATT
  ar3_est_list_3[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
ar3res_att_3 <- do.call("rbind",ar3_att_list_3)
ar3res_est_3 <- do.call("rbind",ar3_est_list_3)

ar3_coverage_3 <- ar3res_att_3 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                        ifelse(RR.upper>1 & RR.lower>1,0,1)))


# Plot histogram of effects

arcoverplotdf <- as.data.frame(cbind(ar3res_att_2, ar3res_att_3))
colnames(arcoverplotdf) <- c("RR_2", "RR_l2", "RR_u2", "RR_3", "RR_l3", "RR_u3")
bw1 <-  with(arcoverplotdf, 2 * IQR(log(RR_2)) / length(log(RR_2))^(1/3))
ar_hist1 <- ggplot(arcoverplotdf, aes(x=log(RR_2))) + geom_histogram(aes(y=..density..),binwidth = bw1,
                                                                fill="gray70",
                                                                color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(coverplotdf$RR_2))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(coverplotdf$RR_2)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("AR(3), K=2")

bw2 <-  with(arcoverplotdf, 2 * IQR(log(RR_3)) / length(log(RR_3))^(1/3))
ar_hist2 <- ggplot(arcoverplotdf, aes(x=log(RR_3))) + geom_histogram(aes(y=..density..),binwidth = bw2,
                                                                fill="gray70",
                                                                color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(arcoverplotdf$RR_3))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(arcoverplotdf$RR_3)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("AR(3), K=3")



# ------------ Simulate using another, longer dataset -------------- #


long_att_list_2 = list()
long_est_list_2 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim_long)
  simres  <- countSynth(data=simdats[["simdata"]],
                        predictors=c("Y"),
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=81,
                        min_1se=F,
                        K=2,
                        full.model=FALSE)
  long_att_list_2[[i]] <- simres$ATT
  long_est_list_2[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
longres_att_2 <- do.call("rbind",long_att_list_2)
longres_est_2 <- do.call("rbind",long_est_list_2)

long_coverage_2 <- longres_att_2 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                              ifelse(RR.upper>1 & RR.lower>1,0,1)))

long_att_list_3 = list()
long_est_list_3 = list()
for (i in 1:nreps) {
  print(i)
  simdats <- factor.sim(testsim_long)
  simres  <- countSynth(data=simdats[["simdata"]],
                        predictors=c("Y"),
                        dependent=c("Y"),
                        unit.variable="ID",
                        time.variable = "TIME",
                        treatment.identifier = simdats[["treated.id"]],
                        controls.identifier = simdats[["controls.id"]],
                        t_int=81,
                        min_1se=F,
                        K=3,
                        full.model=FALSE)
  long_att_list_3[[i]] <- simres$ATT
  long_est_list_3[[i]] <- simres$Estimates$est.counterfactual.cross.fitted
  
}
longres_att_3 <- do.call("rbind",long_att_list_3)
longres_est_3 <- do.call("rbind",long_est_list_3)

long_coverage_3 <- longres_att_3 %>% mutate(coverage = ifelse(RR.lower<1 & RR.upper<1,0,
                                                              ifelse(RR.upper>1 & RR.lower>1,0,1)))

# Plots

lcoverplotdf <- as.data.frame(cbind(longres_att_2, longres_att_3))
colnames(lcoverplotdf ) <- c("RR_2", "RR_l2", "RR_u2", "RR_3", "RR_l3", "RR_u3")
bw1 <-  with(lcoverplotdf , 2 * IQR(log(RR_2)) / length(log(RR_2))^(1/3))
l_hist1 <- ggplot(lcoverplotdf , aes(x=log(RR_2))) + geom_histogram(aes(y=..density..),binwidth = bw1,
                                                                     fill="gray70",
                                                                     color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(lcoverplotdf$RR_2))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(coverplotdf$RR_2)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("Longer pre-period, K=2")

bw2 <-  with(lcoverplotdf, 2 * IQR(log(RR_3)) / length(log(RR_3))^(1/3))
l_hist2 <- ggplot(lcoverplotdf , aes(x=log(RR_3))) + geom_histogram(aes(y=..density..),binwidth = bw2,
                                                                     fill="gray70",
                                                                     color="black") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sd(log(lcoverplotdf$RR_3))), size=1) + theme_bw() +
  geom_vline(xintercept=mean(log(lcoverplotdf$RR_3)), linetype="dashed", size=1) +
  ylab("Density") + xlab("ln(RR)") + scale_x_continuous(limits=(c(-3,3))) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + ggtitle("Longer pre-period, K=3")


################################
# PRINT COVERAGE PROBABILITIES #
################################

#-- Coverage rates with CIs

## Vision zero data (t0 = 26, empirical residuals version)

# K = 2
binom.test(sum(coverage_2$coverage), nreps)

# K = 3
binom.test(sum(coverage_3$coverage), nreps)

## Vision zero data (t0 = 26, AR(3) version)

# K = 2
binom.test(sum(ar3_coverage_2$coverage), nreps)

# K = 3
binom.test(sum(ar3_coverage_3$coverage), nreps)


## SYG data (t0 = 81)

# K = 2
binom.test(sum(long_coverage_2$coverage), nreps)

# K = 3
binom.test(sum(long_coverage_3$coverage), nreps)


##################
# Produce Fig C1 #
##################

plot_grid(hist1, hist2, ar_hist1, ar_hist2, l_hist1, l_hist2, ncol=2, nrow=3, align="hv", axis="b")

