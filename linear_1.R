library(MASS)
library(ISLR)
fix(Boston)
names(Boston)

#linear
attach(Boston)
lstat2 <- (lstat-mean(lstat))
medv2 <- medv-mean(medv)
lm.fit = lm(medv~(lstat-mean(lstat)),data=Boston)

lm.fit = lm(medv2~lstat2)
lm.fit
summary(lm.fit)
names(lm.fit)
lm.fit$coefficients
coef(lm.fit)
confint(lm.fit)
predict(lm.fit,data.frame(lstat=c(5,10,15))
        ,interval = "confidence")
predict(lm.fit,data.frame(lstat=c(5,10,15))
        ,interval = "prediction")
plot(lstat,medv)
abline(lm.fit)
plot(hatvalues(lm.fit))
which.max(hatvalues(lm.fit))

my_linear <- function(y,x){
  x_bar <- mean(x)
  y_bar <- mean(y)
  beta_1 <- (sum((x-x_bar)*(y-y_bar)))/(sum((x-x_bar)^2))
  beta_0 <- y_bar - beta_1*x_bar
  return(paste("Intercept:", toString(beta_0),", Slope:", toString(beta_1),sep=" "))
}

#multiple linear
lm.fit = lm(medv~lstat+age, data=Boston)
summary(lm.fit)
lm.fit = lm(medv~., data=Boston)
summary(lm.fit)
summary(lm.fit)$sigma 
lm.fit1 = lm(medv~.-age,data=Boston)
summary(lm.fit1)

#non-linear 
lm.fit2=lm(medv~lstat+I(lstat^2))
summary (lm.fit2)
lm.fit=lm(medv~lstat)
anova(lm.fit,lm.fit2)
par(mfrow=c(2,2))
plot(lm.fit2)
