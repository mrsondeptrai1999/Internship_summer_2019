table_all_beta_conf = matrix(0,nrow = no_T1,ncol = no_conf)
table_all_p_conf = matrix(0,nrow = no_T1,ncol = no_conf)
table_all_beta_geo <- matrix(0,nrow = no_T1,ncol = no_geo)
table_all_p_geo <- matrix(0,nrow = no_T1,ncol = no_geo)

for (my_row in 1:no_T1){
  df4 <- df1
  df4['y'] <-df2[my_row]
  df4 <- mean_impute(df4) # Mean impute
  df4 <- center_colmeans(df4) # Mean center
  lm.fit2<-lm(y~.,data=df4)
  beta_values = lm.fit2$coefficients[2:(no_conf+1)]
  p_value = summary(lm.fit2)$coefficients[2:(no_conf+1),4]
  table_all_beta_conf[my_row,] <- beta_values
  table_all_p_conf[my_row,] <- p_value
}

for (my_row in 1:no_T1){
  df4 <- df3
  df4['y'] <-df2[my_row]
  df4 <- mean_impute(df4) # Mean impute
  df4 <- center_colmeans(df4) # Mean center
  lm.fit2<-lm(y~.,data=df4)
  beta_values = lm.fit2$coefficients[2:(no_geo+1)]
  p_value = summary(lm.fit2)$coefficients[2:(no_geo+1),4]
  table_all_beta_geo[my_row,] <- beta_values
  table_all_p_geo[my_row,] <- p_value
}

write.table(table_all_beta_conf,'table_all_beta_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_all_p_conf,'table_all_p_conf.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_all_beta_geo,'table_all_beta_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)
write.table(table_all_p_geo,'table_all_p_geo.csv',sep=',',row.names = FALSE,col.names = FALSE)

