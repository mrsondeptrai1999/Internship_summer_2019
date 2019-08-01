install.packages('dummies')
df = read.csv('dummy.tsv',sep='\t')
head(df)
library(dummies)
unique(df$stuff)
b <- dummy(df$stuff)

my_grouping_21000<-function(a){
  for (i in 1:length(a)){
    if(! is.na(a[i])){
      if (a[i]-1000 >=0){
        a[i] = a[i]%%1000
    }
    } 
  }
  return(a)
}

prob_imputation <- function(df){
  for (i in 1:length(df)){
    while (is.na(df[i])){
      df[i] <- df[floor(runif(1,1,length(df)))]
    }
  }
  return (df)
}

mean_impute <- function(data){
  for(i in 1:ncol(data)) {
    data[ , i][is.na(data[ , i])] <- mean(data[ , i], na.rm = TRUE)  
    
  }
  return (data)
} # for loop to impute

# cor() to determine the correlation

'
            X24003.0.0 X24004.0.0
X24003.0.0  1.0000000  0.9141271
X24004.0.0  0.9141271  1.000000
'
'
           X24005.0.0 X24007.0.0 X24006.0.0 X24008.0.0
X24005.0.0  1.0000000  0.5081112  0.4921478  0.7924319
X24007.0.0  0.5081112  1.0000000  0.5274889  0.3675723
X24006.0.0  0.4921478  0.5274889  1.0000000  0.1416952
X24008.0.0  0.7924319  0.3675723  0.1416952  1.0000000
'

# 24006 separate, average out 24005,24007,24006
# Average 24003,24004
# remove 24016,24017,24018
# 24506 - 24505