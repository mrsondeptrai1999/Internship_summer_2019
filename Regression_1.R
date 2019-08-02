source('my_function.R')
library(carData)
library(dummies)
library(car)

# tsv files must be in same folder with this


# DATA CLEANING FOR CONFOUNDING VARIABLES--------------------------------

# Clean data in file Conf_Variable_Visit_0_0.tsv

df1 <- read.table('Conf_Variable_Visit_0_0.tsv',sep='\t',header=TRUE)
df2 <- read.table('Conf_Variable.tsv',sep='\t',header=TRUE)

# Take instance 2 of field_id 54 and create dummy variables
df2_54 <- df2[4] 
df2_54_2 <- dummy(df2_54$X54.2.0,sep='_')
df2_54_2 <- data.frame(df2_54_2)
df2_54_2[1] <- (df2_54_2[1]-df2_54_2[2])
df2_54_2 <- df2_54_2[-2] 

# Create dummy variabble for field id 21000
df2_21000 <- df2[37]

#function change small subgroup to large subgroup(1001 to 1 for example)
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

# Imputation
prob_imputation <- function(df){
  for (i in 1:length(df)){
    while (is.na(df[i])){
      df[i] <- df[floor(runif(1,1,length(df)))]
    }
  }
  return (df)
} #specify which column

df2_21000_2 <- prob_imputation(df2_21000$X21000.0.0)
df2_21000_2 <- my_grouping_21000(df2_21000_2) 
df2_21000_3 <- dummy(df2_21000_2,sep='_')
df2_21000_3 <- data.frame(df2_21000_3)
df2_21000_3 <- df2_21000_3[-1] #Dont need 1 column, choose column -3

# Add the above to data frame
df1_2 <- df1[-c(2,9)]
df1_2 <- cbind(df1_2,df2_54_2)
df1_2 <- cbind(df1_2,df2_21000_3)

# Add brother and sister to create sibling column
df1_3 <- df1_2
df1_3[8] <- df1_3[8] + df1_3[9]
df1_3 <- df1_3[-9]

#head(df1_3)
write.table(df1_3,'conf_table_ver_2.csv',sep=',',row.names = FALSE)
print('conf')
#----------------------------------------------------------------
 
  
  
# DATA CLEANING FOR GEOSPATIAL VARIABLES--------------------------------

#----------------------------------------------------------------

  
  
# EXPORT T1MRI TABLES--------------------------------

my_read('T1MRI_Variable.tsv')

#----------------------------------------------------------------
  

  
# ADJUSTING DATA FOR CONFOUNDING FACTORS--------------------------------

df1_2 <- read.table('conf_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('value_T1MRI.csv',sep = ',',header=TRUE)
df2 <- df2[1]

df1 <- mean_impute(df1) # Mean impute
df1 <- center_colmeans(df1) # Mean center
df2 <- mean_impute(df2) # Mean impute
df2 <- center_colmeans(df2) # Mean center

no_T1 <- dim(df2)[2]
no_conf <- dim(df1)[2]

table_all_beta_conf <- matrix(0,nrow = no_T1,ncol = no_conf)

for (my_row in 1:no_T1){
  df4 <- df1
  df4['y'] <-df2[my_row]
  lm.fit2<-lm(y~.,data=df4)
  beta_values = lm.fit2$coefficients[2:(no_conf+1)]
  table_all_beta_conf[my_row,] <- beta_values
}

# change to numeric matrix to apply matrix mult
y = t(data.matrix(df2)) # 1 x no_subject  
x = t(data.matrix(df1)) # no_conf x no_subject
beta = table_all_beta_conf # 1 x no_conf
epsilon <- y - beta %*% x

print('epsilon')
# Sanity check that the dimension is correct
#dim(df1)
#dim(df2)
#dim(table_beta_conf)

#my_vif(1)

# Non-linear model
table_beta_conf_quad <- matrix(0,nrow = no_T1,ncol = no_conf)
for (my_row in 1:no_T1){
  for (my_col in 1:no_conf){
    df <- data.frame(matrix(ncol=0,nrow=21855))
    df["y"] <- df2[my_row]
    df["x"] <- df1[my_col]
    lm.fit <- lm(y~x+I(x^2),data=df)
    beta_value = lm.fit$coefficients[2]
    table_beta_conf_quad[my_row,my_col] = beta_value
  }
} # fill in the tables for beta
my_AIC<-function(T1_id){
  df4 <- df1
  df4['y'] <-df2[T1_id]
  lm.fit2<-lm(y~X21000_1+I(X21000_1^2),data=df4)
  return(AIC(lm.fit2))
}

