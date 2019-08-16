# Get biggest beta estimates and p-value
df1 <- read.table('my_output/my_beta_quad_geo_2_Z_scored.csv',sep=',',header=FALSE)
df2 <- read.table('my_output/my_p_quad_geo_2_Z_scored.csv',sep=',',header=FALSE)
mm <- data.matrix(df1)
p_val <- data.matrix(df2)
biggest <- tail(order(mm),3)
x <- read.table('my_data/geo_table_ver_2.csv',sep=',',header=TRUE)
x <- x[-1]
y <- read.table('my_data/value_T1MRI.csv',sep=',',header=TRUE)

x <- mean_impute(x) # Mean impute
x <- center_colmeans(x) # Mean center
y <- mean_impute(y) # Mean impute
y <- center_colmeans(y) # Mean center
x <- data.frame(scale(x))
y <- data.frame(scale(y))

row = c()
col = c()
counter = 1
for (num in biggest){
  row[counter] <- num%%137
  col[counter] <- num%/%137 + 1
  counter <- counter+1
}

'for (i in 1:3){
  print(df2[row[i],col[i]])
} # pick first one (1,9)
'
new.data <- data.frame(dist = seq(from = min(data.matrix(x[col[1]])),
                                  to = max(data.matrix(x[col[1]])), length.out = 2000))
plot(data.matrix(y[row[1]])~data.matrix(x[col[1]]),
     main = paste('p-value: ',formatC(df2[row[1],col[1]], format = "e", digits = 5),
                  ', beta estimate:',formatC(df1[row[1],col[1]], format = "e", digits = 5)),
     xlab = 'Distance to coast',
     ylab = 'Volume of grey matter in Amygdala (left)')

lines(new.data$dist*df1[row[1],col[1]]~new.data$dist,col='blue')


# plot 25794(84) and 25914(123)
df1_2 <- read.table('my_data/geo_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('my_data/value_T1MRI.csv',sep = ',',header=TRUE)
df_beta <- read.table('my_output/my_beta_geo_Z_scored.csv',sep = ',',header=FALSE)
df_p <- read.table('my_output/my_p_geo_adjusted_Z_scored.csv',sep = ',',header=FALSE)
df1 <- mean_impute(df1) # Mean impute
df1 <- center_colmeans(df1) # Mean center
df2 <- mean_impute(df2) # Mean impute
df2 <- center_colmeans(df2) # Mean center

df1 <- data.frame(scale(df1))
df2 <- data.frame(scale(df2))


for (num in 1:length(df1)){
  y <- data.matrix(df2[84])
  x <- data.matrix(df1[num])
  new.data <- data.frame(dist = seq(from = min(x),to = max(x), length.out = 200))
  png(paste("plots_new/Precentral_Gyrus_left_",toString(num),".png",sep='')) 
  plot(y~x,
       main = paste('p-value linear: ',formatC(df_p[84,num], format = "e", digits = 5),
                    ', beta estimate linear:',formatC(df_beta[84,num],format = "e", digits = 5),
                    sep=''),
       xlab = colnames(df1)[num],ylab = 'Volume of grey matter in Precentral Gyrus (left)')
  lines(data.matrix(new.data*df_beta[84,num]) ~ data.matrix(new.data), col = "red")
  legend('topleft', legend=c("Linear", "Quadratic"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  dev.off()
}

for (num in 1:length(df1)){
  y <- data.matrix(df2[123])
  x <- data.matrix(df1[num])
  new.data <- data.frame(dist = seq(from = min(x),to = max(x), length.out = 200))
  png(paste("plots_new/VIIIb_Cerebellum_right_",toString(num),".png",sep='')) 
  plot(y~x,
       main = paste('p-value linear: ',formatC(df_p[123,num], format = "e", digits = 5),
                    ', beta estimate linear:',formatC(df_beta[123,num],format = "e", digits = 5),
                    sep=''),
       xlab = colnames(df1)[num],ylab = 'Volume of grey matter in VIIIb Cerebellum (right)')
  lines(data.matrix(new.data*df_beta[123,num]) ~ data.matrix(new.data), col = "red")
  legend('topleft', legend=c("Linear", "Quadratic"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  dev.off()
}