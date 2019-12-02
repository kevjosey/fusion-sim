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
  
  set.seed(07271989)
  
  library(sandwich)
  library(survey)
  library(cbal)
  
  source("D:/Github/target-sim/simfun.R")
  
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
  
  datFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cgen/simData/", 
                       n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  save(simDat, file = datFilename)
  
  idx <- 1:iter # simulation iteration index
  estList <- sapply(idx, simfit, simDat = simDat)
  
  misc_out <- data.frame(n = rep(n, times = iter),
                         y_scen = rep(y_scen, times = iter),
                         z_scen = rep(z_scen, times = iter),
                         s_scen = rep(s_scen, times = iter),
                         PATE = do.call(c, estList[3,]),
                         stringsAsFactors = FALSE)
  
  tau_tmp <- do.call(rbind, estList[1,])
  cp_tmp <- do.call(c, estList[2,])
  colnames(tau_tmp) <- c("GLM", "OUT", "AIPW", "ENT")
  names(cp_tmp) <- c("CAL")
  
  tau <- data.frame(misc_out, tau_tmp, stringsAsFactors = FALSE)
  cp <- data.frame(misc_out, cp_tmp, stringsAsFactors = FALSE)
  
  tauFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cgen/tauHat/", n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  coverageFilename <- paste("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cgen/coverageProb/", n, y_scen, z_scen, s_scen, ".RData", sep = "_")
  
  save(tau, file = tauFilename)
  save(cp, file = coverageFilename)
  
} )

stopCluster(cl)

stop <- Sys.time()
stop - start

### Output

library(ggplot2)
library(gridExtra)

dir_1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cgen/tauHat/"
dir_2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Data/cgen/coverageProb/"

files <- list.files(dir_1)
out_1 <- matrix("", nrow = length(files), ncol = 9)
out_2 <- matrix("", nrow = length(files), ncol = 9)
out_3 <- matrix("", nrow = length(files), ncol = 6)

colnames(out_1) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "GLM", "OUT", "AIPW", "CAL")
colnames(out_2) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "GLM", "OUT", "AIPW", "CAL")
colnames(out_3) <- c("n", "y_scen", "z_scen", "s_scen", "PATE", "CAL")
j <- 1

for (fdx in files) {
  
  file_est <- paste0(dir_1, fdx)
  file_coverage <- paste0(dir_2, fdx)
  load(file_est)
  load(file_coverage)
  
  lbl <- do.call(c, lapply(tau[1,1:4], as.character))
  
  out_1[j,1:4] <- out_2[j,1:4] <- out_3[j,1:4] <- lbl
  out_1[j,5] <- out_2[j,5] <- out_3[j,5] <- round(mean(tau[,5]), 2)
  
  est <- apply(tau[,6:ncol(tau)], 2, function(x) mean(x[!is.infinite(x)], na.rm = TRUE))
  se <- apply(tau[,6:ncol(tau)], 2, function(x) sd(x[!is.infinite(x)], na.rm = TRUE))
  est_se <- sapply(1:length(est), function(i,...) 
    paste(round(est[i], 2), " (", round(se[i], 2), ")", sep = ""))
  
  cover <- round(mean(cp[,6], na.rm = TRUE), 3)
  PATE <- round(mean(tau[,5]), 2)
  
  mse <- apply(tau[,6:ncol(tau)], 2, function(x, PATE) mean((x[!is.infinite(x)] - PATE)^2, na.rm = TRUE), PATE = PATE)
  bias <- apply(tau[,6:ncol(tau)], 2, function(x, PATE) mean((x[!is.infinite(x)] - PATE), na.rm = TRUE), PATE = PATE)
  mse_bias <- sapply(1:length(mse), function(i,...) 
    paste(round(mse[i], 2), " (", round(bias[i], 2), ")", sep = ""))
  
  out_1[j,6:ncol(out_1)] <- est_se
  out_2[j,6:ncol(out_2)] <- mse_bias
  out_3[j,6] <- cover
  
  j <- j + 1
  
}

# plot outcomes

# png("D:/Dropbox (ColoradoTeam)/Projects/Transportability/Output/Figures/ATE_plot.png", 
#     width = 1000, 
#     height = 1000,
#     res = 100, 
#     units = "px")
# 
# grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
# 
# dev.off()

out_1 <- as.data.frame(out_1)
out_1[,1] <- as.numeric(as.character(out_1[,1]))
out_1 <- out_1[order(out_1$n, out_1$y_scen, out_1$z_scen, out_1$s_scen),]

out_2 <- as.data.frame(out_2)
out_2[,1] <- as.numeric(as.character(out_2[,1]))
out_2 <- out_2[order(out_2$n, out_2$y_scen, out_2$z_scen, out_2$s_scen),]

out_3 <- as.data.frame(out_3)
out_3[,1] <- as.numeric(as.character(out_3[,1]))
out_3 <- out_3[order(out_3$n, out_3$y_scen, out_3$z_scen, out_3$s_scen),]

filename1 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cgen/Tables/estimates.csv"
filename2 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cgen/Tables/mse_bias.csv"
filename3 <- "D:/Dropbox (ColoradoTeam)/JoseyDissertation/Output/cgen/Tables/coverageProbs.csv"

write.csv(out_1, file = filename1)
write.csv(out_2, file = filename2)
write.csv(out_3, file = filename3)
