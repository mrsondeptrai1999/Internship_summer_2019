import numpy as np
import pandas as pd
import csv

def my_read(directory):
    df = pd.read_csv(directory, sep="\t")
    
    subject_id = df['eid']
    print(subject_id)
    
    field_id = list(df)
    print(field_id)
    
    df2 = df.drop('eid',axis=1) 
    df2.to_csv('value2.csv', header=False, index=False)
    
    subject_id.to_csv('subject_id.csv',index=False)
  
    a= open('field_id.csv',"w+")
    for i in field_id:
        #if i != 'eid':
        a.write(str(i) + '\n')   
 
    
    return(None)
      

def my_write(field_id,value,subject_id,directory):
    df1 = pd.read_csv(field_id, sep="\t")

    df2 = pd.read_csv(value, sep=",",names = df1['eid'].tolist())

    df3 = pd.read_csv(subject_id, sep=",",names=['eid'])

    df2.insert(loc=0, column='eid', value=df3)

    df2.to_csv(directory,index=False)


