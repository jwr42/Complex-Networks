"""
This program creates adjacency matrices of network files downloaded from
the High-quality INTeractomes (HINT) website: 
http://hint.yulab.org/

@author: Jonathan Roberts
"""

# =============================================================================
# Please provide the address of the HINT file and a filename for the output:

file = ''
# e.g. '../Downloads/HomoSapiens_htb_hq.txt'
filename = ''
# e.g.'Human_Adjacency_Matrix'

# =============================================================================

# Importing relevant modules
import pandas as pd
import numpy as np
import time
# Starting timer
Start = time.clock()
# Importing HINT data from file
df = pd.read_csv(file,sep='\t')
# Select the two relevant columns
col_A = df['Uniprot_A']
col_B = df['Uniprot_B']
# List all datapoints (but this includes duplicates)
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
np.save(filename,A_matrix)
# Printing the runtime
Stop = time.clock()
print("Program Run Time: %.2f mins"%((Stop-Start)/60))