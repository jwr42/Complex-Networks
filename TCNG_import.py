"""
This program creates adjacency matrices of network files downloaded from
The Cancer Network Galaxy (TCNG) website: 
http://tcng.hgc.jp/

@author: Jonathan Roberts (https://github.com/jwr42)
"""

# =============================================================================
# Please provide the address of the TCNG file and a filename for the output:

data = ''
#e.g. '../Downloads/GSE11135_TCNG.txt'
name = ''
#e.g. 'Leukemia_Adjacency_Matrix'

# =============================================================================

# Importing relevant modules
import pandas as pd
import numpy as np
import time
# Starting timer
Start = time.clock()
# Importing TCNG data from file
df = pd.read_csv(data,sep='\t')
# Select the two relevant columns
col_A = df['Parent']
col_B = df['Child']
# List all datapoint (but this includes duplicates)
List = pd.concat([col_A,col_B])
# Number of edges
N_edges = len(col_A)
# Number of unique datapoints
N = List.nunique() 
# List of all unique datapoints
UniqueList = List.drop_duplicates()
# Correcting the index of unique datapoints
new_index = np.arange(0,N)
UniqueList.index = new_index
# Creating the adjacency matrix
A_matrix = np.zeros((N,N),dtype=int)
# Populating the adjacency matrix with the datapoints
for i in range(N_edges):    
    i_head = 0
    for j in range(N):
        if UniqueList[j] == col_A[i]:
            i_head = j    
    i_foot = 0
    for k in range(N):
        if UniqueList[k] == col_B[i]:
            i_foot = k
    A_matrix[i_head,i_foot]=1
    print(i+1)
# Saving the adjacency matrix
np.save(name,A_matrix)
# Printing the runtime
Stop = time.clock()
print("Program Run Time: %.2f mins"%((Stop-Start)/60))