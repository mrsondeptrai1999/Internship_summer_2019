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

'            X24015.0.0  X24013.0.0  X24011.0.0    X24009.0.0    X24506.0.0    
X24015.0.0  1.000000000  0.84149453 -0.0602832063  0.5099210286 -0.038242630
X24013.0.0  0.841494527  1.00000000 -0.0289006885  0.5312831977 -0.068648301
X24011.0.0 -0.060283206 -0.02890069  1.0000000000  0.0003216369 -0.042148250
X24009.0.0  0.509921029  0.53128320  0.0003216369  1.0000000000  0.001323209
X24506.0.0 -0.038242630 -0.06864830 -0.0421482502  0.0013232085  1.000000000
X24507.0.0 -0.028191497 -0.05612891 -0.0319871393  0.0240725872  0.786546012
X24500.0.0 -0.018243502 -0.04911744 -0.0289425383  0.0103795532  0.942872061
X24503.0.0  0.007761503 -0.02337043 -0.0204448992  0.0389725536  0.804622271
X24501.0.0 -0.045020139 -0.01553550  0.0123599832 -0.0389719545 -0.735139406
X24504.0.0 -0.089159845 -0.05925977  0.0157445464 -0.0680838334 -0.537146455
X24507.0.0  X24500.0.0   X24503.0.0  X24501.0.0  X24504.0.0
X24015.0.0 -0.02819150 -0.01824350  0.007761503 -0.04502014 -0.08915984
X24013.0.0 -0.05612891 -0.04911744 -0.023370432 -0.01553550 -0.05925977
X24011.0.0 -0.03198714 -0.02894254 -0.020444899  0.01235998  0.01574455
X24009.0.0  0.02407259  0.01037955  0.038972554 -0.03897195 -0.06808383
X24506.0.0  0.78654601  0.94287206  0.804622271 -0.73513941 -0.53714646
X24507.0.0  1.00000000  0.74578585  0.866706381 -0.62389019 -0.67078370
X24500.0.0  0.74578585  1.00000000  0.854245744 -0.81907863 -0.60386455
X24503.0.0  0.86670638  0.85424574  1.000000000 -0.74522195 -0.81110220
X24501.0.0 -0.62389019 -0.81907863 -0.745221947  1.00000000  0.81270933
X24504.0.0 -0.67078370 -0.60386455 -0.811102199  0.81270933  1.00000000
'