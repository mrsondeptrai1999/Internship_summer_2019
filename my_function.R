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
  lm.fit2<-lm(y~.,data=df4)
  return(vif(lm.fit2))
}

my_AIC<-function(T1_id){
  df4 <- df1
  df4['y'] <-df2[T1_id]
  lm.fit2<-lm(y~.,data=df4)
  return(AIC(lm.fit2))
}

my_BIC<-function(T1_id){
  df4 <- df1
  df4['y'] <-df2[T1_id]
  lm.fit2<-lm(y~.,data=df4)
  return(BIC(lm.fit2))
}

my_read <- function(directory) {
  
  df = read.table(file = directory, sep = '\t', header = FALSE)
  write.table(df[1,-1],'field_id_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  write.table(df[,1],'subject_id_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  write.table(df[,-1],'value_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  
}

my_write <- function(field_id,subject_id,value,directory) {
  df1 = read.table(file = field_id, sep = ',', header = FALSE)
  df2 = read.table(file = subject_id, sep = ',', header = FALSE)
  df3 = read.table(file = value, sep = ',', header = FALSE)
  df3[] <- sapply(df3,as.character) 
  #df4 = rbind(df1,df3)
  df5 = cbind(df2, df3)
  write.table(df5,directory,sep=',',row.names = FALSE,col.names = FALSE)
}
