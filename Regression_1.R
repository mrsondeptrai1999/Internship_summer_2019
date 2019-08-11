#! /apps/well/R/3.4.3/bin/Rscript
#$ -cwd 
#$ -o $HOME/logs
#$ -e $HOME/logs

source('ASSP/my_function.R')
library(dummies)
library(pracma)

# DATA CLEANING FOR CONFOUNDING VARIABLES--------------------------------

# Clean data in file Conf_Variable_Visit_0_0.tsv

df1 <- read.table('my_data/Conf_Variable_Visit_0_0.tsv',sep='\t',header=TRUE)
df2 <- read.table('my_data/Conf_Variable.tsv',sep='\t',header=TRUE)

# Count number of NA
na_count <-sapply(df1, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count
# All columns has around 500 NA or less, insignificant
# Only field Id 845 has 9838 NAs (almost 50%). Also we have another column for age, so I suggest to cut this 

# Take instance 2 of field_id 54(Assesment center) and create dummy variables
df2_54 <- df2[4] 
df2_54_2 <- dummy.data.frame(data = df2_54,dummy.classes='ALL',sep='_')
df2_54_2[1] <- (df2_54_2[1]-df2_54_2[2])
df2_54_2 <- df2_54_2[-2] 

# Create dummy variabble for field id 21000 (Ethnicity)
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

sum(is.na(df2_21000)) # 6 na, insignificant compared to 22000 subjects
df2_21000_2 <- prob_imputation(df2_21000$X21000.0.0) #imputation
df2_21000_2 <- my_grouping_21000(df2_21000_2) 
df2_21000_2 <- data.frame(X21000.0.0 = df2_21000_2)
df2_21000_3 <- dummy.data.frame(data = df2_21000_2,dummy.classes='ALL',sep='_')
df2_21000_3 <- df2_21000_3[-1] #Dont need 1 column, omit column -3

# Add the above to data frame
df1_2 <- df1[-c(2,9)]
df1_2 <- cbind(df1_2,df2_54_2)
df1_2 <- cbind(df1_2,df2_21000_3)

# Add brother and sister to create sibling column
df1_3 <- df1_2
df1_3[8] <- df1_3[8] + df1_3[9]
df1_3 <- df1_3[-9]

colnames(df1_3)[9] <- 'ASCNT'
colnames(df1_3)[2] <- 'Age'
colnames(df1_3)[3] <- 'ngSex'
colnames(df1_3)[c(10:16)] <- c('dEthn_.1', 'dEthn_1','dEthn_2','dEthn_3','dEthn_4','dEthn_5','dEthn_6')
colnames(df1_3)[8] <- 'NumSibs'
colnames(df1_3)[c(4:7)] <- c('gSex','Age_complete_edu','dIncm','BMI')

# Take field_id 738(Income) and create dummy variables
df2_738 <- df1[7] 
sum(is.na(df2_738)) # 380 na, insignificant compared to 22000 subjects
df2_738_2 <- prob_imputation(df2_738$X738.0.0)
df2_738_2 <- data.frame(dIncm = df2_738_2)
df2_738_2 <- dummy.data.frame(data = df2_738_2,dummy.classes='ALL',sep='_')
df2_738_3 <- df2_738_2[-1] #Dont need 1 column, omit column -3
colnames(df2_738_3)[1] <- 'dIncm_.1' # change name to avoid dataframe errors with minus sign
# Add this to data frame
df1_4 <- cbind(df1_3,df2_738_3)
df1_4 <- df1_4[-6] # remove old Income

#Extracting head size
df_T1 <- read.table('my_data/T1MRI_25000.tsv',sep='\t',header = TRUE)
colnames(df_T1)[2] <- 'HeadSize'
df1_5 <- cbind(df1_4,df_T1[2])

# Remove field ID 845 (Age after finishing full time edu) since data is sparse
df1_5 <- df1_5[-5]

# Change dummy of Sex columns to [1 -1] instead of [1 0],get rid of gsex since we dont use it
# Not use to impute since no NA in ngSex
df1_5$ngSex[df1_5$ngSex==0] <- -1
df1_5 <- df1_5[-4]

# Add quadratic columns
df1_5['Age_squared'] <- (df1_5$Age)^2
df1_5['Age_ngSex'] <- df1_5$Age * df1_5$ngSex
df1_5['Age_squared_ngSex'] <- (df1_5$Age)^2 * df1_5$ngSex

write.table(df1_5,'my_data/conf_table_ver_2.csv',sep=',',row.names = FALSE)
print('conf')

#----------------------------------------------------------------



# EXPORT T1MRI TABLES--------------------------------

my_read('my_data/T1MRI_Variables_full.tsv')
print('T1MRI')
#----------------------------------------------------------------



# ADJUSTING DATA FOR CONFOUNDING FACTORS--------------------------------

df1_2 <- read.table('my_data/conf_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('my_data/value_T1MRI.csv',sep = ',',header=TRUE)

df1 <- mean_impute(df1) # Mean impute
df1 <- center_colmeans(df1) # Mean center
df2 <- mean_impute(df2) # Mean impute
df2 <- center_colmeans(df2) # Mean center

no_T1 <- dim(df2)[2]
no_conf <- dim(df1)[2]

table_all_beta_conf_2 <- matrix(0,nrow = no_T1,ncol = no_conf)

# My own beta_estimate loop
for (my_row in 1:no_T1){
  y = t(data.matrix(df2[my_row]))  
  X = t(data.matrix(df1))
  beta_values = y %*% pinv(X)
  table_all_beta_conf_2[my_row,] <- beta_values
}

write.table(table_all_beta_conf_2,'my_output/my_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)

# change to numeric matrix to apply matrix mult
y = t(data.matrix(df2)) # 1 x no_subject  
x = t(data.matrix(df1)) # no_conf x no_subject
beta = table_all_beta_conf_2 # 1 x no_conf
epsilon <- y - beta %*% x

#print('epsilon')
write.table(epsilon,'my_output/my_epsilon.csv',sep=',',row.names = FALSE,col.names = FALSE)

#-------------------------------------



# DATA CLEANING FOR GEOSPATIAL VARIABLES--------------------------------

df1 <- read.table('my_data/Geo_Variable_Visit_0_0.tsv',sep='\t',header=TRUE)
df2 <- read.table('my_data/Geo_Variable_2.tsv',sep='\t',header=TRUE)
# did not download everything in first file so have to download the missing data is second file

df2_2 <- df2[-c(13:15)] # these field _id unpacked in first file so delete
df2_2 <- df2_2[-c(6:8)] # remove 24016,24017,24108 since we have data in 2010

# air pollution-----------------------------------------------
df2_air_pol <- df2_2[c(2,3,4,5,6,7)] # Deal with air pollution
sum(is.na(df2_air_pol))
#first 4 cols have 752 NA, last 2 have 176 NA, which are all insignificant
df2_air_pol_2 <- mean_impute(df2_air_pol) # mean impute to calculate corr
cor(df2_air_pol_2)

# get this
'           X24005.0.0 X24007.0.0 X24006.0.0 X24008.0.0 X24003.0.0 X24004.0.0
X24005.0.0  1.0000000  0.5081112  0.4921478  0.7924319  0.4415832  0.4625774
X24007.0.0  0.5081112  1.0000000  0.5274889  0.3675723  0.6390328  0.5617424
X24006.0.0  0.4921478  0.5274889  1.0000000  0.1416952  0.8552657  0.8341973
X24008.0.0  0.7924319  0.3675723  0.1416952  1.0000000  0.1059970  0.1554351
X24003.0.0  0.4415832  0.6390328  0.8552657  0.1059970  1.0000000  0.9141271
X24004.0.0  0.4625774  0.5617424  0.8341973  0.1554351  0.9141271  1.0000000'
# 05,08 seems to have high correlation to we average them out (pm10,2.5-10um)
# 06,03,04,07 seems to have high correlation to we average them out ((pm2.5) absorbance,(pm2.5),...)

df2_3 <- df2_2
df2_3[2] <- (df2_3[2]+df2_3[5])/2
df2_3[3] <- (df2_3[3]+df2_3[4]+df2_3[6]+df2_3[7])/4
df2_3 <- df2_3[-c(4,5,6,7)] # average out 24005,07,08
#-----------------------------------------------

df2_3 <- df2_3[-4] #delete 24019

# Traffic and Greenspace -----------------------------------
df2_4 <- df2_3[-c(6,8,10,12,14,16,18,20)] #delete second visit since I forgot to flag -v 0 in funpack
df2_4 <- df2_4[-c(11,12)] # not need water
df2_5 <- cbind(df1,df2_4[-1])
df2_misc <- df2_5[c(20,13,14,15,21,22,23,24,25,26)] 
sum(is.na(df2_misc))
# 176 NAs for first 4 cols, 185 NAs for next 2 cols, 886 for last 4 cols
# All are relatively small compared to 20000 so imputation is safe
df2_misc_2 <- mean_impute(df2_misc)
cor(df2_misc_2)

df2_6 <- df2_5

# Traffic: 24015,13,09 has high enough correlation so we average them out, 
# 11 has low cor with the other 3 so it is seperated
df2_6[20] <- (df2_6[20]+df2_6[13]+df2_6[15])/3
# Natural environment percentage(24506,24507): cor is around 0.79, high enough to average
df2_6[21] <- (df2_6[21]+df2_6[22])/2
# Greenspace percentage(24500,24503): cor is around 0.85, high enough to avr
df2_6[23] <- (df2_6[23]+df2_6[24])/2
# Domestic garden percentage(24501,24504): cor is around 0.81, high enough to average
df2_6[25] <- (df2_6[25]+df2_6[26])/2

df2_6 <- df2_6[-c(13,15,22,24,26)]
df2_6 <- df2_6[-c(7,8,9)] # Omit noise pollution
#-------------------------------------------

# count number of NAs in each col
na_count <-sapply(df2_6, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count
# Highest NA in oe column is 1263, which is insignificant compared to around 22000 subjects
# So no column is too sparse

#--------------------------------------------

# Rename
colnames(df2_6)<- c('eid','home_loc_east','home_loc_north','PoB_east','PoB_north','Townsend',
                    'close_maj_road','Distance_nearest_road_1','Distance_nearest_road_2','Traffic_1',
                    'dist_to_coast','pop_density','air_pol_1','air_pol_2','Traffic_2','Greenspace_1',
                    'Greenspace','Greenspace_2')
# Greensapce and natural env have high pos correlation so we can avr
# but above and domestic garden high neg cor so not avr
df2_6[16] <- (df2_6[16]+df2_6[17])/2
df2_6 <- df2_6[-17]

# Using inverse of nearest major road and nearest road instead
df2_6[8] <- mean_impute(df2_6[8])
df2_6[9] <- mean_impute(df2_6[9])
cor(df2_6$Distance_nearest_road_1, df2_6$Distance_nearest_road_2)
#around 0.1251037 so low correlation

a <- mean_impute(df2_6[2]) 
b <- mean_impute(df2_6[4])
cor(a,b) # around 0.3
a <- mean_impute(df2_6[3]) 
b <- mean_impute(df2_6[5])
cor(a,b) # around 0.4
#low correltion betwen location variables

# Delete 'close to major road' since we have inverse distance to road 
df2_6<-df2_6[-7]

df2_7 <- df2_6
df_pop <- df2_7[11]
df_pop <- prob_imputation(df_pop$pop_density)
df_pop_2 <- data.frame(pop_density = df_pop)
df_pop_2 <- dummy.data.frame(data = df_pop_2,dummy.classes='ALL',sep='_')
df_pop_3 <- df_pop_2[2]+df_pop_2[3]+df_pop_2[7]+df_pop_2[8]
colnames(df_pop_3)<-'pop_density_urban'

df2_7 <- cbind(df2_7,df_pop_3)
df2_7 <- df2_7[-11]
#--------------------------------------------

write.table(df2_7,'my_data/geo_table_ver_2.csv',sep=',',row.names = FALSE)
print('geo')

#----------------------------------------------------------------


# ESTIMATE BETA VALUES FOR T1MRI AGAINST GEO AND TEST THE HYPOTHESIS--------
df1_2 <- read.table('my_data/geo_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('my_data/value_T1MRI.csv',sep = ',',header=TRUE)

df1 <- mean_impute(df1) # Mean impute
df1 <- center_colmeans(df1) # Mean center
df2 <- mean_impute(df2) # Mean impute
df2 <- center_colmeans(df2) # Mean center

no_T1 <- dim(df2)[2]
no_geo <- dim(df1)[2]

table_beta_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_p_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_AIC_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_BIC_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_t_stats_geo <- matrix(0,nrow = no_T1,ncol = no_geo)

# My Beta_estimate loop
for (my_row in 1:no_T1){
  for(my_col in 1:no_geo){
  y = t(data.matrix(df2[my_row]))  
  X = t(data.matrix(df1[my_col]))
  beta_values = y %*% pinv(X)
  table_beta_geo[my_row,my_col] <- beta_values
  }
}
write.table(table_beta_geo,'my_output/my_beta_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)

# My t_statistics,p_value loop for linear model
for (my_row in 1:no_T1){
  for(my_col in 1:no_geo){
    y <- t(data.matrix(df2[my_row]))
    x <- t(data.matrix(df1[my_col]))
    beta_value <- table_beta_geo[my_row,my_col]
    my_rss = sum((y - beta_value %*% x)^2)
    s <- sqrt(my_rss/(dim(df1)[1]-no_geo))
    test_statistics <- beta_value/(s*sqrt(pinv(x%*%t(x))))
    table_t_stats_geo[my_row,my_col] <- test_statistics
    p_values <- 2*(1 - pt(abs(test_statistics),dim(df1)[1]-no_geo))
    table_p_geo[my_row,my_col] <- p_values
  }
}
write.table(table_p_geo,'my_output/my_p_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_t_stats_geo,'my_output/my_t_stats_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)

# Create table for residual
y = t(data.matrix(df2)) 
x = t(data.matrix(df1)) 
beta = data.matrix(read.table('my_output/my_beta_geo.csv',sep = ',',header=FALSE))
table_residual_geo <- y - beta %*% x
table_residual_geo_2 <- t(table_residual_geo)

write.table(table_residual_geo_2,'my_output/my_residual_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)

# My AIC,BIC loop
for (my_row in 1:no_T1){
  for(my_col in 1:no_geo){
    y <- t(data.matrix(df2[my_row]))
    x <- t(data.matrix(df1[my_col]))
    beta_value <- table_beta_geo[my_row,my_col]
    my_rss = sum((y - beta_value %*% x)^2)
    AIC_value <- my_AIC(df2[my_row],my_rss)
    BIC_value <- my_BIC(df2[my_row],my_rss)
    table_AIC_geo[my_row,my_col] <- AIC_value
    table_BIC_geo[my_row,my_col] <- BIC_value
  }
}
write.table(table_AIC_geo,'my_output/my_AIC_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_BIC_geo,'my_output/my_BIC_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)

table_beta_quad_geo_1 <- matrix(0,nrow = no_T1,ncol = no_geo)
table_beta_quad_geo_2 <- matrix(0,nrow = no_T1,ncol = no_geo)
table_p_quad_geo_1 <- matrix(0,nrow = no_T1,ncol = no_geo)
table_p_quad_geo_2 <- matrix(0,nrow = no_T1,ncol = no_geo)
table_AIC_quad_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_BIC_quad_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_t_stats_quad_geo_1 <- matrix(0,nrow = no_T1,ncol = no_geo)
table_t_stats_quad_geo_2 <- matrix(0,nrow = no_T1,ncol = no_geo)

# My beta, p-value,AIC, BIC loop for quadratic model
for (my_row in 1:no_T1){
  for(my_col in 1:no_geo){
    y <- t(data.matrix(df2[my_row]))
    df_x <- df1[my_col]
    df_x['squared'] <- I(df1[my_col]^2)
    df_x <- center_colmeans(df_x)
    x <- t(data.matrix(df_x))
    beta_values = y %*% pinv(x)
    table_beta_quad_geo_1[my_row,my_col] <- beta_values[1]
    table_beta_quad_geo_2[my_row,my_col] <- beta_values[2]
    my_rss = sum((y - beta_values %*% x)^2)
    s <- sqrt(my_rss/(dim(df1)[1]-2))
    test_statistics_1 <- beta_values[1]/(s*sqrt(solve(x%*%t(x))[1,1]))
    table_t_stats_quad_geo_1[my_row,my_col] <- test_statistics_1
    p_value_1 <- 2*(1 - pt(abs(test_statistics_1),dim(df1)[1]-2))
    table_p_quad_geo_1[my_row,my_col] <- p_value_1
    test_statistics_2 <- beta_values[2]/(s*sqrt(solve(x%*%t(x))[2,2]))
    table_t_stats_quad_geo_2[my_row,my_col] <- test_statistics_2
    p_value_2 <- 2*(1 - pt(abs(test_statistics_2),dim(df1)[1]-2))
    table_p_quad_geo_2[my_row,my_col] <- p_value_2
    AIC_value <- my_AIC(df2[my_row],my_rss)
    BIC_value <- my_BIC(df2[my_row],my_rss)
    table_AIC_quad_geo[my_row,my_col] <- AIC_value
    table_BIC_quad_geo[my_row,my_col] <- BIC_value
  }
}
write.table(table_beta_quad_geo_1,'my_output/my_beta_quad_geo_1.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_beta_quad_geo_2,'my_output/my_beta_quad_geo_2.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_p_quad_geo_1,'my_output/my_p_quad_geo_1.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_p_quad_geo_2,'my_output/my_p_quad_geo_2.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_AIC_quad_geo,'my_output/my_AIC_quad_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_BIC_quad_geo,'my_output/my_BIC_quad_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_t_stats_quad_geo_1,'my_output/my_t_stats_quad_geo_1.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_t_stats_quad_geo_2,'my_output/my_t_stats_quad_geo_2.csv',sep=',',row.names = FALSE,col.names = FALSE)

# Create table for residual
y = t(data.matrix(df2)) 
x_1 = t(data.matrix(df1))
x_2 = t(data.matrix(df1))^2
beta_1 = data.matrix(read.table('my_output/my_beta_quad_geo_1.csv',sep = ',',header=FALSE))
beta_2 = data.matrix(read.table('my_output/my_beta_quad_geo_2.csv',sep = ',',header=FALSE))
table_residual_geo <- y - beta_1 %*% x_1 - beta_2 %*% x_2
table_residual_geo_2 <- t(table_residual_geo)

write.table(table_residual_geo_2,'my_output/my_residual_quad_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)


# Adjust p-values using Benjamini-Hochberg
p_values_linear <- read.table('my_output/my_p_geo.csv',sep=',')
p_values_linear_adjusted <- matrix(0,nrow = dim(p_values_linear)[1],ncol = dim(p_values_linear)[2])
for (my_col in 1:dim(p_values_linear)[2]){
  p_adjusted <- p.adjust(data.matrix(p_values_linear[,my_col]),method='BH')
  p_values_linear_adjusted[,my_col] <- p_adjusted
}
write.table(p_values_linear_adjusted,'my_output/my_p_geo_adjusted.csv',sep=',',row.names = FALSE,col.names = FALSE)

p_values_quad_1 <- read.table('my_output/my_p_quad_geo_1.csv',sep=',')
p_values_quad_1_adjusted <- matrix(0,nrow = dim(p_values_quad_1)[1],ncol = dim(p_values_quad_1)[2])
for (my_col in 1:dim(p_values_quad_1)[2]){
  p_adjusted <- p.adjust(data.matrix(p_values_quad_1[,my_col]),method='BH')
  p_values_quad_1_adjusted[,my_col] <- p_adjusted
}
write.table(p_values_quad_1_adjusted,'my_output/my_p_quad_geo_1_adjusted.csv',sep=',',row.names = FALSE,col.names = FALSE)

p_values_quad_2 <- read.table('my_output/my_p_quad_geo_2.csv',sep=',')
p_values_quad_2_adjusted <- matrix(0,nrow = dim(p_values_quad_2)[1],ncol = dim(p_values_quad_2)[2])
for (my_col in 1:dim(p_values_quad_2)[2]){
  p_adjusted <- p.adjust(data.matrix(p_values_quad_2[,my_col]),method='BH')
  p_values_quad_2_adjusted[,my_col] <- p_adjusted
}
write.table(p_values_quad_2_adjusted,'my_output/my_p_quad_geo_2_adjusted.csv',sep=',',row.names = FALSE,col.names = FALSE)

#----------------------------------------------------

# DETERMINE WHICH MODEL WE CHOOSE USING AIC AND BIC

df1 <- read.table('my_output/my_AIC_geo.csv',sep=',',header=FALSE)
df2 <- read.table('my_output/my_AIC_quad_geo.csv',sep=',',header=FALSE)
table_compare_AIC_linear <- matrix(0,nrow = dim(df1)[1],ncol = dim(df1)[2])
# 0 is linear model, 1 is quadratic model
for (my_row in 1:dim(table_compare)[1]){
  for (my_col in 1:dim(table_compare)[2]){
    if(df2[my_row,my_col] < df1[my_row,my_col]){
      table_compare_AIC_linear[my_row,my_col] <- 1
    }
    if(df2[my_row,my_col] >= df1[my_row,my_col]){
      table_compare_AIC_linear[my_row,my_col] <- 0
    }
  }
}

df3 <- read.table('my_output/my_BIC_geo.csv',sep=',',header=FALSE)
df4 <- read.table('my_output/my_BIC_quad_geo.csv',sep=',',header=FALSE)
table_compare_BIC_linear <- matrix(0,nrow = dim(df3)[1],ncol = dim(df3)[2])
# 0 is linear model, 1 is quadratic model
for (my_row in 1:dim(table_compare)[1]){
  for (my_col in 1:dim(table_compare)[2]){
    if(df4[my_row,my_col] < df3[my_row,my_col]){
      table_compare_BIC_linear[my_row,my_col] <- 1
    }
    if(df4[my_row,my_col] >= df3[my_row,my_col]){
      table_compare_BIC_linear[my_row,my_col] <- 0
    }
  }
}

# Create scatter plot of hippocampus vs geo

df1_2 <- read.table('my_data/geo_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('my_data/value_T1MRI.csv',sep = ',',header=TRUE)
df1 <- mean_impute(df1) # Mean impute
df2 <- mean_impute(df2) 
for (num in 1:length(df1)){
  y <- data.matrix(df2[47])
  x <- data.matrix(df1[num])
  lm.fit <- lm(y~x)
  #abline(lm.fit,col="red", lwd=3, lty=2)
  lm.fit2 <- lm(y~x + I(x^2))
  new.data <- data.frame(dist = seq(from = min(x),to = max(x)), length.out = 200)
  pred_lm <- predict(lm.fit, df2[47])
  pred_lm2 <- predict(lm.fit2, df2[47])
  pdf(paste("plots/rplot_left_hippo_",toString(num),".pdf",sep='')) 
  plot(y~x)
  lines(pred_lm ~ x, col = "red")
  lines(pred_lm2 ~ x, col = "blue")
  dev.off()
}
for (num in 1:length(df1)){
  y <- data.matrix(df2[48])
  x <- data.matrix(df1[num])
  lm.fit <- lm(y~x)
  #abline(lm.fit,col="red", lwd=3, lty=2)
  lm.fit2 <- lm(y~x + I(x^2))
  new.data <- data.frame(dist = seq(from = min(x),to = max(x)), length.out = 200)
  pred_lm <- predict(lm.fit, df2[48])
  pred_lm2 <- predict(lm.fit2, df2[48])
  pdf(paste("plots/rplot_right_hippo_",toString(num),".pdf",sep='')) 
  plot(y~x)
  lines(pred_lm ~ x, col = "red")
  lines(pred_lm2 ~ x, col = "blue")
  dev.off()
}

# Create scatter plot of hippocampus vs geo z-scored
df1_2 <- read.table('my_data/geo_table_ver_2.csv',sep = ',',header=TRUE)
df1 <- df1_2[-1]
df2 <- read.table('my_data/value_T1MRI.csv',sep = ',',header=TRUE)
df_beta <- read.table('my_output/my_beta_geo_Z_scored.csv',sep = ',',header=FALSE)
df_p <- read.table('my_output/my_p_geo_Z_scored.csv',sep = ',',header=FALSE)
df_beta_2 <- read.table('my_output/my_beta_quad_geo_2_Z_scored.csv',sep = ',',header=FALSE)
df_beta_1 <- read.table('my_output/my_beta_quad_geo_1_Z_scored.csv',sep = ',',header=FALSE)
df_p_2 <- read.table('my_output/my_p_quad_geo_2_Z_scored.csv',sep = ',',header=FALSE)

df1 <- mean_impute(df1) # Mean impute
df1 <- center_colmeans(df1) # Mean center
df2 <- mean_impute(df2) # Mean impute
df2 <- center_colmeans(df2) # Mean center

df1 <- data.frame(scale(df1))
df2 <- data.frame(scale(df2))


for (num in 1:length(df1)){
  y <- data.matrix(df2[47])
  x <- data.matrix(df1[num])
  lm.fit <- lm(y~x)
  lm.fit2 <- lm(y~x + I(x^2))
  new.data <- data.frame(dist = seq(from = min(x),to = max(x), length.out = 200))
  #pdf(paste("plots/rplot_left_hippo_",toString(num),".pdf",sep='')) 
  plot(y~x,
       main = paste('Scatter plot, p-value linear: ',formatC(df_p[47,num], format = "e", digits = 5),
                    ', beta estimate linear:',formatC(df_beta[47,num],format = "e", digits = 5),
                    ', \n p-value estimate quad:',formatC(df_p_2[47,num],format = "e", digits = 5),
                    ', beta estimate quad:',formatC(df_beta_2[47,num],format = "e", digits = 5),
                    sep=''),
       xlab = colnames(df1)[num],ylab = colnames(df2)[47])
  lines(data.matrix(new.data*df_beta[47,num]) ~ data.matrix(new.data), col = "red")
  lines(data.matrix(df_beta_1[47,num]*new.data+df_beta_2[47,num]*new.data^2) ~ data.matrix(new.data), col = "blue")
  legend('topleft', legend=c("Linear", "Quadratic"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  #dev.off()
}

for (num in 1:length(df1)){
  y <- data.matrix(df2[48])
  x <- data.matrix(df1[num])
  lm.fit <- lm(y~x)
  lm.fit2 <- lm(y~x + I(x^2))
  new.data <- data.frame(dist = seq(from = min(x),to = max(x), length.out = 200))
  #pdf(paste("plots/rplot_right_hippo_",toString(num),".pdf",sep='')) 
  plot(y~x,
       main = paste('Scatter plot, p-value: ',formatC(df_p[48,num], format = "e", digits = 5),
                    ', beta estimate:',formatC(df_beta[48,num], format = "e", digits = 5),
                    ', \n p-value estimate quad:',formatC(df_p_2[47,num],format = "e", digits = 5),
                    ', beta estimate quad:',formatC(df_beta_2[47,num],format = "e", digits = 5),
                    sep=''),
       xlab = colnames(df1)[num],ylab = colnames(df2)[48])
  lines(data.matrix(new.data*df_beta[48,num]) ~ data.matrix(new.data), col = "red")
  lines(data.matrix(df_beta_1[48,num]*new.data+df_beta_2[48,num]*new.data^2) ~ data.matrix(new.data), col = "blue")
  legend('topleft', legend=c("Linear", "Quadratic"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  #dev.off()
}  


