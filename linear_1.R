library(MASS)
library(ISLR)
fix(Boston)
names(Boston)

lm.fit = lm(medv~lstat,data=Boston)
attach(Boston)
lm.fit = lm(medv~lstat)
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


my_linear <- function(y,x){
  x_bar <- mean(x)
  y_bar <- mean(y)
  beta_1 <- (sum((x-x_bar)*(y-y_bar)))/(sum((x-x_bar)^2))
  beta_0 <- y_bar - beta_1*x_bar
  return(paste("Intercept:", toString(beta_0),", Slope:", toString(beta_1),sep=" "))
}


