################################
### PURPOSE: Simulation Code ###
### BY: Kevin Josey          ###
################################

library(snow)

iter <- 1000
n <- c(1000, 2000)
y_scen <- c("a", "b")
z_scen <- c("a", "b")
s_scen <- c("a", "b")

simConditions <- expand.grid(n, y_scen, z_scen, s_scen, stringsAsFactors = FALSE)
names(simConditions) <- c("n", "y_scen", "z_scen", "s_scen")
index <- 1:nrow(simConditions)

## multicore simulation
cl <- makeCluster(3, type = "SOCK")

clusterEvalQ(cl, {
  
  set.seed(09301987)
  
  library(sandwich)
  library(survey)
  library(cbal)
  
  source("D:/Github/combine-sim/transport.R")
  source("D:/Github/combine-sim/tmle.R")
  source("D:/Github/combine-sim/simfun.R")
  
})

clusterExport(cl = cl, list = list("simConditions", "iter"), envir = environment())

start <- Sys.time()

clusterApply(cl, index, function(i) {
  
  dat <- simConditions[i,]
  
  sig2 <- 5
  n <- dat$n
  z_scen <- dat$z_scen
  y_scen <- dat$y_scen
  s_scen <- dat$s_scen

  simDat <- replicate(iter, hte_data(n = n, sig2 = sig2, y_scen = y_scen, z_scen = z_scen, s_scen = s_scen))
  
  datFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/simData/", 
                       n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  save(simDat, file = datFilename)
  
  idx <- 1:iter # simulation iteration index
  estList <- sapply(idx, simfit, simDat = simDat)
  
  PATE <- mean(do.call(c, estList[3,]))
  misc_out <- data.frame(n = rep(n, times = iter),
                         y_scen = rep(y_scen, times = iter),
                         z_scen = rep(z_scen, times = iter),
                         s_scen = rep(s_scen, times = iter),
                         PATE = rep(PATE, times = iter),
                         stringsAsFactors = FALSE)
  
  tau_tmp <- do.call(rbind, estList[1,])
  se_tmp <- do.call(rbind, estList[2,])
  cp_tmp <- matrix(NA, nrow = iter, ncol = 2)
  cp_tmp[,1] <- as.numeric(tau_tmp[,4] - se_tmp[,1]*1.96 <= PATE & tau_tmp[,4] + se_tmp[,1]*1.96 >= PATE)
  cp_tmp[,2] <- as.numeric(tau_tmp[,6] - se_tmp[,2]*1.96 <= PATE & tau_tmp[,6] + se_tmp[,2]*1.96 >= PATE)
  colnames(tau_tmp) <- c("GLM", "OUT", "TMLE", "CAL-T", "ACAL", "CAL-F")
  names(cp_tmp) <- c("CAL-T", "CAL-F")
  
  tau <- data.frame(misc_out, tau_tmp, stringsAsFactors = FALSE)
  cp <- data.frame(misc_out, cp_tmp, stringsAsFactors = FALSE)
  
  tauFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/", n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  coverageFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/coverageProb/", n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  
  save(tau, file = tauFilename)
  save(cp, file = coverageFilename)
  
} )

stopCluster(cl)

stop <- Sys.time()
stop - start

### Output

library(ggplot2)
library(gridExtra)

dir_1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/"
dir_2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/coverageProb/"

files <- list.files(dir_1)
out_1 <- matrix("", nrow = length(files), ncol = 11)
out_2 <- matrix("", nrow = length(files), ncol = 11)
out_3 <- matrix("", nrow = length(files), ncol = 7)

colnames(out_1) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "GLM", "OUT", "TMLE", "CAL-T", "ACAL", "CAL-F")
colnames(out_2) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "GLM", "OUT", "TMLE", "CAL-T", "ACAL", "CAL-F")
colnames(out_3) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "CAL-T", "CAL-F")
j <- 1

for (fdx in files) {
  
  file_est <- paste0(dir_1, fdx)
  file_coverage <- paste0(dir_2, fdx)
  load(file_est)
  load(file_coverage)
  
  lbl <- do.call(c, lapply(tau[1,1:4], as.character))
  
  out_1[j,1:4] <- out_2[j,1:4] <- out_3[j,1:4] <- lbl
  out_1[j,5] <- out_2[j,5] <- out_3[j,5] <- round(mean(tau[,5]), 2)
  
  est <- apply(tau[,6:ncol(tau)], 2, mean, na.rm = TRUE)
  se <- apply(tau[,6:ncol(tau)], 2, sd, na.rm = TRUE)
  est_se <- sapply(1:length(est), function(i,...) 
    paste(round(est[i], 2), " (", round(se[i], 2), ")", sep = ""))
  
  cover <- round(colMeans(cp[,6:7], na.rm = TRUE), 3)
  PATE <- round(mean(tau[,5]), 2)
  
  mse <- apply(tau[,6:ncol(tau)], 2, function(x, PATE) mean((x - PATE)^2, na.rm = TRUE), PATE = PATE)
  bias <- apply(tau[,6:ncol(tau)], 2, function(x, PATE) mean((x - PATE), na.rm = TRUE), PATE = PATE)
  mse_bias <- sapply(1:length(mse), function(i,...) 
    paste(round(mse[i], 2), " (", round(bias[i], 2), ")", sep = ""))
  
  out_1[j,6:ncol(out_1)] <- est_se
  out_2[j,6:ncol(out_2)] <- mse_bias
  out_3[j,6:ncol(out_3)] <- cover
  
  j <- j + 1
  
}

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_a_a_a_.RData")
dat1 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat1$ind <- factor(dat1$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p1 <- ggplot(dat1) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50)  +
  xlab("") +
  ggtitle("outcome: a, treatment: a, sampling: a") +
  geom_hline(yintercept = 24.25, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_a_b_a_.RData")
dat2 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat2$ind <- factor(dat2$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p2 <- ggplot(dat2) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50) +
  xlab("") +
  ggtitle("outcome: a, treatment: b, sampling: a") +
  geom_hline(yintercept = 24.9, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_b_a_a_.RData")
dat3 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat3$ind <- factor(dat3$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p3 <- ggplot(dat3) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50)  +
  xlab("") +
  ggtitle("outcome: b, treatment: a, sampling: a") +
  geom_hline(yintercept = 23.5, colour = "red", linetype = 3, size = 1, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_b_b_a_.RData")
dat4 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat4$ind <- factor(dat4$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p4 <- ggplot(dat4) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50)  +
  xlab("") +
  ggtitle("outcome: b, treatment: b, sampling: a") +
  geom_hline(yintercept = 23.5, colour = "red", linetype = 3, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_b_a_b_.RData")
dat5 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat5$ind <- factor(dat5$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p5 <- ggplot(dat5) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50)  +
  xlab("") +
  ggtitle("outcome: b, treatment: a, sampling: b") +
  guides(fill =  FALSE) +
  geom_hline(yintercept = 23, colour = "red", linetype = 3, show.legend = FALSE) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

load("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/fusion/tauHat/_1000_a_b_b_.RData")
dat6 <- stack(as.data.frame(tau[,6:ncol(tau)])[,-4])
dat6$ind <- factor(dat6$ind, labels = c("IOS", "OM", "AIOS", "CAL"))
p6 <- ggplot(dat6) + 
  geom_boxplot(aes(x = ind, y = values, fill = ind)) + 
  ylab("PATE") +
  ylim(0, 50)  +
  xlab("") +
  ggtitle("outcome: a, treatment: b, sampling: b") +
  geom_hline(yintercept = 24.3, colour = "red", linetype = 3, show.legend = FALSE) + 
  guides(fill =  FALSE) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

png("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/fusion/Figures/ATE_plot.png", 
    width = 3000, 
    height = 3000,
    res = 300, 
    units = "px")

grid.arrange(p1, p3, p6, p4, ncol = 2, nrow = 2)

dev.off()

out_1 <- as.data.frame(out_1)
out_1[,1] <- as.numeric(as.character(out_1[,1]))
out_1 <- out_1[order(out_1$n, out_1$y_scen, out_1$z_scen, out_1$s_scen),]

out_2 <- as.data.frame(out_2)
out_2[,1] <- as.numeric(as.character(out_2[,1]))
out_2 <- out_2[order(out_2$n, out_2$y_scen, out_2$z_scen, out_2$s_scen),]

out_3 <- as.data.frame(out_3)
out_3[,1] <- as.numeric(as.character(out_3[,1]))
out_3 <- out_3[order(out_3$n, out_3$y_scen, out_3$z_scen, out_3$s_scen),]

filename1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/fusion/Tables/estimates.csv"
filename2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/fusion/Tables/mse_bias.csv"
filename3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/fusion/Tables/coverageProbs.csv"

write.csv(out_1, file = filename1)
write.csv(out_2, file = filename2)
write.csv(out_3, file = filename3)
