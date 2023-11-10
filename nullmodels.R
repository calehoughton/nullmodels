library(vegan)
install.packages("gam")
library(gam)
install.packages("scam")
library(scam)
install.packages("gdm")
library(gdm)
library(betapart)
library(tidyverse)
install.packages("zetadiv")
library(zetadiv)

load("islands_null_regression_models_practical.RData")
beta.pair(dat.soc.pa)

help("permatfull")
dat.soc.pa.perm.none <- permatfull(dat.soc.pa, mtype = "prab", fixedmar = "none", times = 99)
dat.soc.pa.perm.row <- permatfull(dat.soc.pa, mtype = "prab", fixedmar = "rows", times = 99)
dat.soc.pa.perm.col<- permatfull(dat.soc.pa, mtype = "prab", fixedmar = "columns", times = 99)
dat.soc.pa.perm.both <- permatfull(dat.soc.pa, mtype = "prab", fixedmar = "both", times = 99)

load("permutations_hawaii.RData") 

beta.mean <- function(dat){ ##this function computes the average for Sorensen and Simpson beta diversity for a givent site-by-species matrix dat
  beta.dat <- beta.pair(dat)
  return(c(mean(beta.dat$beta.sor),mean(beta.dat$beta.sim)))
}

beta.rand.soc.none <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.none$perm,beta.mean)),99,2,byrow = TRUE)) ## this applies the beta.mean() function above to each permutated site-by-species matrix
names(beta.rand.soc.none) <- c("Sorensen","Simpson") ## rename the columns of the data frame

beta.rand.soc.row <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.row$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.row) <- c("Sorensen","Simpson")

beta.rand.soc.col <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.col$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.col) <- c("Sorensen","Simpson")

beta.rand.soc.both <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.both$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.both) <- c("Sorensen","Simpson")

help(lapply)
help(unlist)
help(matrix)
help(data.frame)

beta.soc.obs <- beta.pair(dat.soc.pa)
beta.soc.sim <- mean(beta.soc.obs$beta.sim)
beta.soc.sor <- mean(beta.soc.obs$beta.sor)

par(mfrow=c(2,4))
hist(beta.rand.soc.none$Sorensen,breaks=seq(0,1,0.05),main="None fixed",xlab="Sorensen")
abline(v=beta.soc.sor,col="red")
hist(beta.rand.soc.row$Sorensen,breaks=seq(0,1,0.05),main="Rows fixed",xlab="Sorensen")
abline(v=beta.soc.sor,col="red")
hist(beta.rand.soc.col$Sorensen,breaks=seq(0,1,0.05),main="Columns fixed",xlab="Sorensen")
abline(v=beta.soc.sor,col="red")
hist(beta.rand.soc.both$Sorensen,breaks=seq(0,1,0.05),main="Both fixed",xlab="Sorensen")
abline(v=beta.soc.sor,col="red")

hist(beta.rand.soc.none$Simpson,breaks=seq(0,1,0.05),main="",xlab="Simpson")
abline(v=beta.soc.sim,col="red")
hist(beta.rand.soc.row$Simpson,breaks=seq(0,1,0.05),main="",xlab="Simpson")
abline(v=beta.soc.sim,col="red")
hist(beta.rand.soc.col$Simpson,breaks=seq(0,1,0.05),main="",xlab="Simpson")
abline(v=beta.soc.sim,col="red")
hist(beta.rand.soc.both$Simpson,breaks=seq(0,1,0.05),main="",xlab="Simpson")
abline(v=beta.soc.sim,col="red")

#alpha
islands.pred.soc$alpha <- rowSums(dat.soc.pa)
islands.pred.haw$alpha <- rowSums(dat.haw.pa)

mod.glm.soc <- glm(alpha~IslandArea,data=islands.pred.soc,family = poisson())
mod.glm.haw <- glm(alpha~IslandArea,data=islands.pred.haw,family = poisson())
summary(mod.glm.soc)
summary(mod.glm.haw)

cor(islands.pred.soc$alpha,predict(mod.glm.soc,type="response"))
cor(islands.pred.haw$alpha,predict(mod.glm.haw,type="response"))

mod.gam.soc <- gam(alpha~s(IslandArea),data=islands.pred.soc,family = poisson(),method = "REML")
mod.gam.haw <- gam(alpha~s(IslandArea),data=islands.pred.haw,family = poisson(),method = "REML")

cor(islands.pred.soc$alpha,predict(mod.gam.soc,type="response"))
cor(islands.pred.haw$alpha,predict(mod.gam.haw,type="response"))

par(mfrow=c(1,2))
plot(mod.gam.soc)
points(islands.pred.soc$IslandArea,log(islands.pred.soc$alpha)-mod.gam.soc$coefficients[1])
plot(mod.gam.haw,ylim=c(-6,3))
points(islands.pred.haw$IslandArea,log(islands.pred.haw$alpha)-mod.gam.haw$coefficients[1])

mod.gam.haw <- gam(alpha~s(IslandArea),data=islands.pred.haw,family = poisson(),method = "REML")
mod.gam.haw2 <- gam(alpha~s(IslandArea,k=4),data=islands.pred.haw,family = poisson(),method = "REML")
mod.gam.soc2 <- gam(alpha~s(IslandArea,k=4),data=islands.pred.soc,family = poisson(),method = "REML")


cor(islands.pred.soc$alpha,predict(mod.gam.soc2,type="response"))
cor(islands.pred.haw$alpha,predict(mod.gam.haw2,type="response"))

par(mfrow=c(1,2))
plot(mod.gam.soc2)
points(islands.pred.soc$IslandArea,log(islands.pred.soc$alpha)-mod.gam.soc2$coefficients[1])
plot(mod.gam.haw2,ylim=c(-6,3))
points(islands.pred.haw$IslandArea,log(islands.pred.haw$alpha)-mod.gam.haw$coefficients[1])

mod.glm.soc <- glm(alpha~elev_max+log(IslandArea),data=islands.pred.soc,family = poisson())
mod.glm.haw <- glm(alpha~elev_max+log(IslandArea),data=islands.pred.haw,family = poisson())

mod.gam.soc <- gam(alpha~s(log(IslandArea)),data=islands.pred.soc,family = poisson(),method = "REML")
mod.gam.soc2 <- gam(alpha~s(log(IslandArea),k=5),data=islands.pred.soc,family = poisson(),method = "REML")
mod.gam.haw <- gam(alpha~s(log(IslandArea)),data=islands.pred.haw,family = poisson(),method = "REML")
mod.gam.haw2 <- gam(alpha~s(log(IslandArea),k=5),data=islands.pred.haw,family = poisson(),method = "REML")

summary(mod.glm.soc)
summary(mod.glm.haw)

cor(islands.pred.soc$alpha,predict(mod.glm.soc,type="response"))
cor(islands.pred.haw$alpha,predict(mod.glm.haw,type="response"))

summary(mod.gam.soc)
summary(mod.gam.haw)

cor(islands.pred.soc$alpha,predict(mod.gam.soc,type="response"))
cor(islands.pred.haw$alpha,predict(mod.gam.haw,type="response"))

cor(islands.pred.soc$alpha,predict(mod.gam.soc2,type="response"))
cor(islands.pred.haw$alpha,predict(mod.gam.haw2,type="response"))

par(mfrow=c(2,2))
plot(mod.gam.soc)
points(log(islands.pred.soc$IslandArea),log(islands.pred.soc$alpha)-mod.gam.soc$coefficients[1])
plot(mod.gam.haw,ylim=c(-6,3))
points(log(islands.pred.haw$IslandArea),log(islands.pred.haw$alpha)-mod.gam.haw$coefficients[1])
plot(mod.gam.soc2)
points(log(islands.pred.soc$IslandArea),log(islands.pred.soc$alpha)-mod.gam.soc2$coefficients[1])
plot(mod.gam.haw2,ylim=c(-6,3))
points(log(islands.pred.haw$IslandArea),log(islands.pred.haw$alpha)-mod.gam.haw$coefficients[1])

