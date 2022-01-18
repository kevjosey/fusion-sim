################################
### PURPOSE: Simulation Code ###
### BY: Kevin Josey          ###
################################

library(parallel)

iter <- 1000
n <- c(500, 1000, 2000)
scenario <- c("base", "ps-mis", "out-mis", "ps-exchange",
              "out-exchange", "sample-overlap", "treat-overlap")

simConditions <- expand.grid(n, scenario, stringsAsFactors = FALSE)
names(simConditions) <- c("n", "scenario")
index <- 1:nrow(simConditions)

## multicore simulation

set.seed(42)
  
library(sandwich)
library(survey)

source("~/Github/fusion-sim/transport.R")
source("~/Github/fusion-sim/tmle.R")
source("~/Github/fusion-sim/aug.R")
source("~/Github/fusion-sim/simfun.R")

start <- Sys.time()

mclapply(index, function(i,...) {
  
  print(i)
  
  dat <- simConditions[i,]
  
  sig2 <- 1
  n <- dat$n
  scen <- dat$scenario

  simDat <- replicate(iter, gen_data(n = n, sig2 = sig2, scenario = scen))
  
  datFilename <- paste("~/Dropbox/JoseyDissertation/Data/fusion/simData/", 
                       n, scen, ".RData", sep = "_")
  save(simDat, file = datFilename)
  
  idx <- 1:iter # simulation iteration index
  estList <- sapply(idx, simfit, simDat = simDat)
  
  PATE <- mean(do.call(c, estList[3,]))
  misc_out <- data.frame(n = rep(n, times = iter),
                         scenario = rep(scen, times = iter),
                         PATE = rep(PATE, times = iter),
                         stringsAsFactors = FALSE)
  
  tau_tmp <- do.call(rbind, estList[1,])
  se_tmp <- do.call(rbind, estList[2,])
  tau_tmp[is.infinite(tau_tmp) | is.nan(tau_tmp)] <- NA
  se_tmp[is.infinite(se_tmp) | is.nan(se_tmp)] <- NA
  colnames(tau_tmp) <- colnames(se_tmp) <- c("TMLE-T", "AUG-T", "CAL-T", "TMLE-F", "AUG-F","CAL-F")
  
  cp_tmp <- matrix(NA, nrow = iter, ncol = 6)
  cp_tmp[,1] <- as.numeric(tau_tmp[,1] - se_tmp[,1]*1.96 <= PATE & tau_tmp[,1] + se_tmp[,1]*1.96 >= PATE)
  cp_tmp[,2] <- as.numeric(tau_tmp[,2] - se_tmp[,2]*1.96 <= PATE & tau_tmp[,2] + se_tmp[,2]*1.96 >= PATE)
  cp_tmp[,3] <- as.numeric(tau_tmp[,3] - se_tmp[,3]*1.96 <= PATE & tau_tmp[,3] + se_tmp[,3]*1.96 >= PATE)
  cp_tmp[,4] <- as.numeric(tau_tmp[,4] - se_tmp[,4]*1.96 <= PATE & tau_tmp[,4] + se_tmp[,4]*1.96 >= PATE)
  cp_tmp[,5] <- as.numeric(tau_tmp[,5] - se_tmp[,5]*1.96 <= PATE & tau_tmp[,5] + se_tmp[,5]*1.96 >= PATE)
  cp_tmp[,6] <- as.numeric(tau_tmp[,6] - se_tmp[,6]*1.96 <= PATE & tau_tmp[,6] + se_tmp[,6]*1.96 >= PATE)
  names(cp_tmp) <- colnames(tau_tmp)
  
  tau <- data.frame(misc_out, tau_tmp, stringsAsFactors = FALSE)
  cp <- data.frame(misc_out, cp_tmp, stringsAsFactors = FALSE)
  
  tauFilename <- paste("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/", n, scen, ".RData", sep = "_")
  coverageFilename <- paste("~/Dropbox/JoseyDissertation/Data/fusion/coverageProb/", n, scen, ".RData", sep = "_")
  
  save(tau, file = tauFilename)
  save(cp, file = coverageFilename)
  
}, mc.cores = 7)

stop <- Sys.time()
stop - start

### Output

library(ggplot2)
library(gridExtra)

dir_1 <- "~/Dropbox/JoseyDissertation/Data/fusion/tauHat/"
dir_2 <- "~/Dropbox/JoseyDissertation/Data/fusion/coverageProb/"

files <- list.files(dir_1)
out_1 <- matrix("", nrow = length(files), ncol = 9)
out_2 <- matrix("", nrow = length(files), ncol = 9)
out_3 <- matrix("", nrow = length(files), ncol = 9)

colnames(out_1) <- colnames(out_2) <- colnames(out_3) <- 
  c("n", "scenario", "PATE",  "TMLE-T", "AUG-T", "CAL-T", "TMLE-F", "AUG-F", "CAL-F")
j <- 1

for (fdx in files) {
  
  file_est <- paste0(dir_1, fdx)
  file_coverage <- paste0(dir_2, fdx)
  load(file_est)
  load(file_coverage)
  
  lbl <- do.call(c, lapply(tau[1,1:2], as.character))
  
  out_1[j,1:2] <- out_2[j,1:2] <- out_3[j,1:2] <- lbl
  out_1[j,3] <- out_2[j,3] <- out_3[j,3] <- round(mean(tau[,3]), 2)
  PATE <- mean(tau[,3])
  
  est <- apply(tau[,4:ncol(tau)], 2, mean, na.rm = TRUE)
  
  #Replace outliers according to temp_range
  
  se <- apply(tau[,4:ncol(tau)], 2, sd, na.rm = TRUE)
  est_se <- sapply(1:length(est), function(i,...) 
    paste(round(est[i], 2), " (", round(se[i], 2), ")", sep = ""))
  
  cover <- round(colMeans(cp[,4:ncol(cp)], na.rm = TRUE), 3)

  mse <- sqrt(colMeans((tau[,4:ncol(tau)] - PATE)^2, na.rm = TRUE))
  bias <- colMeans(tau[,4:ncol(tau)] - PATE, na.rm = TRUE)
  mse_bias <- sapply(1:length(mse), function(i,...) 
    paste(round(bias[i], 2), " (", round(mse[i], 2), ")", sep = ""))
  
  out_1[j,4:ncol(out_1)] <- est_se
  out_2[j,4:ncol(out_2)] <- mse_bias
  out_3[j,4:ncol(out_3)] <- cover
  
  j <- j + 1
  
}

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_out-mis_.RData")
dat1 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat1$ind <- factor(dat1$ind, 
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p1 <- ggplot(dat1) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-5, 5)  +
  xlab("") +
  ggtitle("Outcome Model Misspecification") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_ps-mis_.RData")
dat2 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat2$ind <- factor(dat2$ind, 
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p2 <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-1, 1) +
  xlab("") +
  ggtitle("Treatment & Sample Misspecification") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_ps-exchange_.RData")
dat3 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat3$ind <- factor(dat3$ind, 
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p3 <- ggplot(dat3) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-5, 5)  +
  xlab("") +
  ggtitle("Propensity Score Exchangeability Violation") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_out-exchange_.RData")
dat4 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat4$ind <- factor(dat4$ind,
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p4 <- ggplot(dat4) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-3, 3)  +
  xlab("") +
  ggtitle("Potential Outcome Exchangeability Violation") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_sample-overlap_.RData")
dat5 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat5$ind <- factor(dat5$ind,
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p5 <- ggplot(dat5) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-15, 15)  +
  xlab("") +
  ggtitle("Sample Overlap Violation") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("~/Dropbox/JoseyDissertation/Data/fusion/tauHat/_1000_treat-overlap_.RData")
dat6 <- stack(as.data.frame(tau[,4:ncol(tau)] - mean(tau[,3])))
dat6$ind <- factor(dat6$ind,
                   levels = c("TMLE.T", "AUG.T", "CAL.T",
                              "TMLE.F", "AUG.F", "CAL.F"),
                   labels = c("TMLE-T", "AUG-T", "CAL-T",
                              "TMLE-F", "AUG-F", "CAL-F"))
p6 <- ggplot(dat6) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("Bias") +
  ylim(-12, 12)  +
  xlab("") +
  ggtitle("Treatment Overlap Violation") +
  geom_hline(yintercept = 0, colour = "red", linetype = 3, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

png("~/Dropbox/JoseyDissertation/Output/fusion/Figures/ATE_plot.png", 
    width = 3000, 
    height = 3000,
    res = 300, 
    units = "px")

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)

dev.off()

out_1 <- as.data.frame(out_1)
out_1[,1] <- as.numeric(as.character(out_1[,1]))
out_1 <- out_1[order(out_1$n, out_1$scenario),]

out_2 <- as.data.frame(out_2)
out_2[,1] <- as.numeric(as.character(out_2[,1]))
out_2 <- out_2[order(out_2$n, out_2$scenario),]

out_3 <- as.data.frame(out_3)
out_3[,1] <- as.numeric(as.character(out_3[,1]))
out_3 <- out_3[order(out_3$n, out_3$scenario),]

filename1 <- "~/Dropbox/JoseyDissertation/Output/fusion/Tables/estimates.csv"
filename2 <- "~/Dropbox/JoseyDissertation/Output/fusion/Tables/mse_bias.csv"
filename3 <- "~/Dropbox/JoseyDissertation/Output/fusion/Tables/coverageProbs.csv"

write.csv(out_1, file = filename1)
write.csv(out_2, file = filename2)
write.csv(out_3, file = filename3)
