#Data Summarization
library(readr)
streamflow <- read_csv("streamflow.csv")
View(streamflow)
colnames(streamflow) <- c("X1","Y","X2","X3","X4","X5","X6","X7","X8")
summary(streamflow)

## HISTOGRAMS
par(mfrow=c(3,3))
hist(streamflow$X1,main="Histogram of X1",xlab="Staid")
hist(streamflow$Y,main="Histogram of Y",xlab="90th  percentile  of  the  time  series")
hist(streamflow$X2,main="Histogram of X2",xlab="the drainage area ")
hist(streamflow$X3,main="Histogram of X3",xlab="average basin precipitation ")
hist(streamflow$X4,main="Histogram of X4",xlab="average basin temperature  ")
hist(streamflow$X5,main="Histogram of X5",xlab="average temperature at the stream location ")
hist(streamflow$X6,main="Histogram of X6",xlab="average relative humidity across the basin")
hist(streamflow$X7,main="Histogram of X7",xlab="average March precipitation")
hist(streamflow$X8,main="Histogram of X8",xlab="median relief ratio")

## BOX PLOTS
boxplot(streamflow$Y, xlab="Max90", main="Max90")
boxplot(streamflow$X1, xlab="STAID", main="STAID")
boxplot(streamflow$X2, xlab="DRAIN_SQKM", main=" The drainage area")
boxplot(streamflow$X3, xlab="PPTAVG_BASIN", main="The average basinprecipitation")
boxplot(streamflow$X4, xlab="T_AVG_BASIN", main="Boxplot of the average basintemperature")
boxplot(streamflow$X5, xlab="T_AVG_SITE", main="Boxplot of the averagetemperature at the stream location")
boxplot(streamflow$X6, xlab="RH_BASIN", main="Boxplot of the average relativehumidity across the basin ")
boxplot(streamflow$X7, xlab="MAR_PPT7100_CM", main="Boxplot of the averageMarch precipitation ")
boxplot(streamflow$X8, xlab="RRMEDIAN", main="Boxplot of the median reliefratio")

## SCATTER PLOTS
plot(streamflow$Y, xlab="Max90",ylab="Max90", main="Boxplot of Stream Flow")
plot(streamflow$X1, xlab="STAID",ylab="STAID", main="STAID")
plot(streamflow$X2, xlab="DRAIN_SQKM", ylab="DRAIN_SQKM", main=" The drainage area")
plot(streamflow$X3, xlab="PPTAVG_BASIN",ylab="PPTAVG_BASIN", main="The average basinprecipitation")
plot(streamflow$X4, xlab="T_AVG_BASIN",ylab="T_AVG_BASIN", main="Boxplot of the average basintemperature")
plot(streamflow$X5, xlab="T_AVG_SITE", ylab="T_AVG_SITE", main="Boxplot of the averagetemperature at the stream location")
plot(streamflow$X6, xlab="RH_BASIN",ylab="RH_BASIN", main="Boxplot of the average relativehumidity across the basin ")
plot(streamflow$X7, xlab="MAR_PPT7100_CM",ylab="MAR_PPT7100_CM", main="Boxplot of the averageMarch precipitation ")
plot(streamflow$X8, xlab="RRMEDIAN",ylab="RRMEDIAN", main="Boxplot of the median reliefratio")

## ADDED VARIABLE PLOTS
library(car)
avPlots(fitstream)

## CORRELATION MATRIX
cor(streamflow)
pairs(streamflow)
##checking for multicollinearity. 
eigen(cor(streamflow))$values

## MODELS AND METHODS 

fitstream<-lm(Y ~X1+X2+X3+X4+X5+X6+X7+X8,data=streamflow)
fitstream
summary(fitstream)

##ANOVA t-test
anova(fitstream)

library(MASS)
## Model Selection
library(olsrr)

#### Print all possible regression models in terms of adjr, Cp, AIC, and BIC.

par(mfrow=c(1,1))
b <- ols_step_all_possible(fitstream)
plot(b)
#### Adjusted R2 ####

b.adjr = data.frame(n=b$n,predictors=b$predictors,adjr=b$adjr)
print(b.adjr)
print(b.adjr[c(93,163,219,247,255),])

#### Cp ####

b.cp = data.frame(n=b$n,predictors=b$predictors,cp=b$cp)
print(b.cp)
print(b.cp[c(93,163,219,247,255),])

#### AIC ####

b.aic = data.frame(n=b$n,predictors=b$predictors,aic=b$aic)
print(b.aic)
print(b.aic[c(93,163,219,247,255),])

#### BIC ####

b.bic = data.frame(n=b$n,predictors=b$predictors,bic=b$sbic)
print(b.bic)
print(b.bic[c(93,163,219,247,255),])

#### PRESS ####

b.press = data.frame(n=b$n,predictors=b$predictors,press=b$msep)
print(b.press)
print(b.press[c(93,163,219,247,255),])

#### Stepwise Regression ####

k <- ols_step_both_p(fitstream,pent=0.10,prem=0.1,details=TRUE)
plot(k)

#### Final Model? ####

reduced.lmfit <- lm(Y ~ X2 + X3+X5, data=streamflow)
summary(reduced.lmfit)

######## Regression Diagnostics ############

res <- rstudent(reduced.lmfit)
fitted.y <- fitted(reduced.lmfit)

######## Residual Plots ##########

par(mfrow=c(2,2))

plot(res ~ streamflow$X2, xlab="X2", ylab="Residual", main="Residuals vs. X2")
abline(h=0)
plot(res ~ streamflow$X3, xlab="X3", ylab="Residual", main="Residuals vs. X3")
abline(h=0)
plot(res ~ streamflow$X5, xlab="X5", ylab="Residual", main="Residuals vs. X5")
abline(h=0)

plot(res ~ fitted.y, xlab="Fitted value", ylab="Residual", main="Residuals vs. Fitted Values")
abline(h=0)

######### Multicollinearity ##########

vif(reduced.lmfit) 


######### Constancy of Error Variances #########
library(lmtest)
bptest(reduced.lmfit)


#Durbin-Watson
#install lmtest
library(lmtest)
dwtest(fitstream, alternative="two.sided")

######### Normality ###########

qqnorm(res);qqline(res)
########Shapiro test########
shapiro.test(res)

#DFFITS values
library(olsrr)

ols_plot_dffits(reduced.lmfit)

#DFBETAS values
ols_plot_dfbetas(reduced.lmfit)

#Cook's distance values
ols_plot_cooksd_chart(reduced.lmfit)

######### Transformation #########
library(EnvStats)

boxcox.summary <- boxcox(reduced.lmfit, optimize=TRUE)
lambda <- boxcox.summary$lambda
lambda 
trans.Y <- streamflow$Y^lambda

streamflow <- cbind(streamflow,trans.Y)
streamflow

######### Re-fitting a model using the transformed response variable. ##########

boxcox.lmfit <- lm(trans.Y ~ X2 + X3 + X5, data=streamflow)
summary(boxcox.lmfit)
anova(boxcox.lmfit)

boxcox.res <- rstudent(boxcox.lmfit)

boxcox.fitted.y <- fitted(boxcox.lmfit)
######## Residual Plots ##########
par(mfrow=c(2,2))

plot(boxcox.res ~ streamflow$X2, xlab="X2", ylab="Residual", main="Residuals vs. X2")
abline(h=0)
plot(boxcox.res ~ streamflow$X3, xlab="X3", ylab="Residual", main="Residuals vs. X3")
abline(h=0)
plot(boxcox.res ~ streamflow$X5, xlab="X5", ylab="Residual", main="Residuals vs. X5")
abline(h=0)

plot(boxcox.res ~ fitted.y, xlab="Fitted value", ylab="Residual", main="Residuals vs. Fitted Values")
abline(h=0)

######### Multicollinearity ##########
library(HH)
vif(boxcox.lmfit)


######### Constancy of Error Variances #########

bptest(boxcox.lmfit)

dwtest(boxcox.lmfit, alternative="two.sided")
######### Normality ###########

qqnorm(boxcox.res);qqline(boxcox.res)
shapiro.test(boxcox.res)

######### Final Model ##########

final.lmfit <- boxcox.lmfit
summary(final.lmfit)
anova(final.lmfit)
##Obtain DFFITS, DFBETAS, and Cook’s distance values
library(olsrr)

#DFFITS values
ols_plot_dffits(final.lmfit)

#DFBETAS values
ols_plot_dfbetas(final.lmfit)

#Cook's distance values
ols_plot_cooksd_chart(final.lmfit)

################ Fit a regression model with interaction terms ##################

streamflow.lmfit <- lm(Y ~ X2 + X3 + X5 + X2*X3 + X2*X5 + X3*X5, data=streamflow)
summary(streamflow.lmfit)
anova(streamflow.lmfit)

################ Fit a regression model with no interaction terms ###############

streamflow.reduced <- lm(Y ~ X2 + X3 + X5, data=streamflow)

################ Test for significance of the interaction terms #################

anova(streamflow.reduced,streamflow.lmfit)

