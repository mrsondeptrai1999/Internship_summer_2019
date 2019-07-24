df1 <- read.csv('value_conf.csv') 
df2 <- read.csv('value_T1MRI.csv')
df3 <- read.csv('value_geo.csv') 


center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
} # Center function

mean_impute <- function(data){
  for(i in 1:ncol(data)) {
    data[ , i][is.na(data[ , i])] <- mean(data[ , i], na.rm = TRUE)  
    
  }
  return (data)
} # for loop to impute


no_T1 <- dim(df2)[2]
no_conf <- dim(df1)[2]
no_geo <- dim(df3)[2]


table_beta_conf <- matrix(0,nrow = no_T1,ncol = no_conf)
table_p_conf <- matrix(0,nrow = no_T1,ncol = no_conf)
table_beta_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_p_geo <- matrix(0,nrow = no_T1,ncol = no_geo)


for (my_row in 1:no_T1){
  for (my_col in 1:no_conf){
    df <- data.frame(matrix(ncol=0,nrow=21855))
    df['y'] <- df2[my_row]
    df['x'] <- df1[my_col]
    df <- mean_impute(df) # Mean impute
    df <- center_colmeans(df) # Mean center
    lm.fit <- lm(y~x,data=df)
    beta_value = lm.fit$coefficients[2]
    p_value = summary(lm.fit)$coefficients[2,4]
    table_beta_conf[my_row,my_col] = beta_value
    table_p_conf[my_row,my_col] = p_value
  }
} # fill in the tables for conf

for (my_row in 1:no_T1){
  for (my_col in 1:no_geo){
    df <- data.frame(matrix(ncol=0,nrow=21855))
    df['y'] <- df2[my_row]
    df['x'] <- df3[my_col]
    df <- mean_impute(df) # Mean impute
    df <- center_colmeans(df) # Mean center
    lm.fit <- lm(y~x,data=df)
    beta_value = lm.fit$coefficients[2]
    p_value = summary(lm.fit)$coefficients[2,4]
    table_beta_geo[my_row,my_col] = beta_value
    table_p_geo[my_row,my_col] = p_value
  }
} # fill in the tables for geo





# Draft
#write.table(table_beta_conf,'table_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)
#write.table(table_p_conf,'table_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)
#write.table(table_beta_geo,'table_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)
#write.table(table_p_geo,'table_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)

'
df <- data.frame(matrix(ncol=0,nrow=20))
df[] <- df2[1]
df[] <- df1[2]
df <- mean_impute(df) # Mean impute
df <- center_colmeans(df) # Mean center
lm.fit <- lm(y~x,data=df)
summary(lm.fit)
beta_value = lm.fit$coefficients[2]
table_beta[1,2] = beta_value
table_beta
'
#df5[colnames(df2[1])] <-  df2[1]
#df5 <- mean_impute(df5) # Mean impute
#df4 <- center_colmeans(df5) # Mean center
#colnames(df2[1])
#X25010.2.0


#df5
#df5 <- mean_impute(df5)
#df5 <- center_colmeans(df5)
#df5
#lm.fit2<-lm(X22704.0.0~X25009.2.0,data=df5)
#summary(lm.fit2)
#plot(df5$X22704.0.0,df5$X25009.2.0)

