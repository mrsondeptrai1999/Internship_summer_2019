center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
} # input: dataframe, output: dataframe with column mean_centered 
#Son Ha, 07/08/19

mean_impute <- function(data){
  for(i in 1:ncol(data)) {
    data[ , i][is.na(data[ , i])] <- mean(data[ , i], na.rm = TRUE)
    
  }
  return (data)
} # input: vector, output: vector with all NAs mean_imputed  
#Son Ha, 07/08/19

my_vif <- function(df){
  ans = diag(solve(cor(df)))
  return(ans)
} # input: dataframe, output: vif of dataframe (the diag of inverse of cor matrix)
#Son Ha, 07/08/19

my_AIC <- function(df,epsilon){
  ans = 2*(length(df)+1) + dim(df)[1]*log(2*pi*epsilon/(dim(df)[1]))  
  # length + 1 since number of beta estimates and error
  return(ans)
}# input: dataframe for x and model's rss; output: AIC
#Son Ha, 07/08/19

my_BIC <- function(df,epsilon){
  ans = log(dim(df)[1])*(length(df)+1) + dim(df)[1]*log((2*pi*epsilon)/(dim(df)[1]))
  # length + 1 since number of beta estimates and error
  return(ans)
}# input: dataframe for x and model's rss; output: BIC
#Son Ha, 07/08/19

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

prob_imputation <- function(df){
  for (i in 1:length(df)){
    while (is.na(df[i])){
      df[i] <- df[floor(runif(1,1,length(df)))]
    }
  }
  return (df)
}# input: dataframe, output: dataframe with NAs imputed using freq of each element
#Son Ha, 07/08/19