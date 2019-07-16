import numpy as np
import pandas as pd
'UKB/Variables/Conf_VariableList.txt'

def csv_to_column (directory, name):
    f = open(directory,'r')
    #print(f.read())
    ans = []
    data = f.readlines()
    for line in data:
         new_string = ''.join(ch for ch in line if ch.isdigit())
         #if new_string != '':
         ans.append(new_string)
    print(ans)
    
    a= open(name,"w+")
    for i in ans:
        a.write(i + '\n')
    return(None)

csv_to_column ('UKB/Variables/Conf_VariableList.txt', 'Conf_VariableList_ID.txt')