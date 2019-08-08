#Testing my BIC and AIC functions
my_rss = sum((y - beta %*% x)^2)
my_BIC(df4,my_rss) #510346.3
my_AIC(df4,my_rss) #510154.4
BIC(lm.fit2) #532201.3
AIC(lm.fit2) #532009.4
# Explanation for the slight difference is that we assume the data is Gaussian distribution and 
# calculte log-likelihood accordingly,
# While the package calculate log-likelihood using logLik, and I assume it does not use Gaussian as
# straightforward as I do

# DRAFT AND NOTE, NOT IN MAIN SCRIPT, DONT RUN----------------------------------
# Check my p-value loop, compare with built-in function
y_test <- t(data.matrix(df2[1]))
x_test <- t(data.matrix(df1[1]))
beta_value <- table_beta_geo[1,1]
my_rss = sum((y_test - beta_value %*% x_test)^2)
s <- sqrt(my_rss/(dim(df1)[1]-no_geo))
test_statistics <- beta_value/(s*sqrt(pinv(x_test%*%t(x_test)))) 
# around -33.6, correct when compared with what we got from lm() down there 
p_value <- 2*(1 - pt(abs(test_statistics),dim(df1)[1]-no_geo))
# some values of p will be too small that R will recognise it as 0
my_AIC(df2[1],my_rss)
my_BIC(df2[1],my_rss)

df4 <- df1[1]
df4['y'] <-df2[1]
lm.fit<-lm(y~.,data=df4)
beta_values = lm.fit2$coefficients[2:(no_conf+1)]
table_all_beta_conf[my_row,] <- beta_values
lm.fit(summary)
AIC(lm.fit)
BIC(lm.fit)


# Check Quadratic model
y_test <- t(data.matrix(df2[1]))
x_test_2 <- t(data.matrix(df1[1]))
df_x_test <- df1[1]
df_x_test['squared'] <- I(df1[1]^2)
df_x_test <- center_colmeans(df_x_test)
x_test <- t(data.matrix(df_x_test))
beta_values = y_test %*% pinv(x_test)
# correct beta_values, 0.01132584 -5.536971e-08 
my_rss = sum((y_test - beta_values %*% x_test)^2)
s <- sqrt(my_rss/(dim(df1)[1]-2))
test_statistics <- beta_values[2]/(s*sqrt(solve(x_test%*%t(x_test))[2,2])) 
#correct t_value 0.01132584 
my_AIC(df_x_test,my_rss)
my_BIC(df_x_test,my_rss)

df4 <- df1[1]
df4['y'] <-df2[1]
colnames(df4)[1] <- 'x'
lm.fit<-lm(y~x+I(x^2),data=df4)
summary(lm.fit)
#-------------------------------------------------------------------------        