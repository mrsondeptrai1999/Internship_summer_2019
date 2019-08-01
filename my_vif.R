(library(car))
(library(carData))

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

my_vif<-function(T1_id){
  df4 <- df1
  df4['y'] <-df2[T1_id]
  df4 <- mean_impute(df4) # Mean impute
  df4 <- center_colmeans(df4) # Mean center
  lm.fit2<-lm(y~.,data=df4)
  return(vif(lm.fit2))
}
#my_vif(1)
