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

my_vif <- function(df){
  ans = diag(solve(cor(df)))
  return(ans)
} # vif is the diag of inverse of cor matrix

my_AIC <- function(df,vector,epsilon){
  ans = 2*(length(df)+1) + dim(df)[1]*log(sum(epsilon)/(dim(df)[1]))  
  # length + 1 since number of beta estimates and error
  return(ans)
}

my_BIC <- function(df,vector,epsilon){
  ans = log(dim(df)[1])*(length(df)+1) + dim(df)[1]*log(sum(epsilon)/(dim(df)[1]))  
  # length + 1 since number of beta estimates and error
  return(ans)
}

my_read <- function(directory) {
  
  df = read.table(file = directory, sep = '\t', header = FALSE)
  write.table(df[1,-1],'my_data/field_id_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  write.table(df[,1],'my_data/subject_id_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  write.table(df[,-1],'my_data/value_T1MRI.csv',sep=',',row.names = FALSE,col.names = FALSE)
  
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