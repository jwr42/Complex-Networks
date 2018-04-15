"""
This file contains a number of functions designed to create model networks, 
help in the analysis of networks and help create plots of the behaviours of 
different networks. All functions only work on network edges that are 
unweighted and undirected.

@author: Jonathan Roberts (https://github.com/jwr42)
"""

import numpy as np
from numpy.random import random as rand

# =============================================================================
# Model Networks
# =============================================================================

def Random(N,P):
    """
    Random creates a model random network.
    Inputs: 
        N = number of vertices
        P = probability of connection
    Output:
        C = adjacency matrix
    """
    A = np.zeros((N,N))
    c=0
    for i in np.arange(0,N,1):
        c+=1
        for j in np.arange(c,N,1):
            random_number = rand()
            if random_number <= P:
                A[i,j] = 1
    B = A.T #transpose of matrix A
    C = np.add(A,B)
    return C

def SmallWorld(N,P): 
    """
    SmallWorld creates a model small world network.
    Inputs:
        N = number of vertices
        P = probability of rewiring edges
    Output:
        A = adjacency matrix
    """
    A = np.zeros((N,N))
    #establishing connections with 2 nearest neighbours
    for i in range(0,N):
        if i == 0:
            A[i,1:3] = A[i,N-2:N] = 1
        if i == 1:
            A[i,0] = A[i,2:4] = A[i,N-1:N] = 1
        if i > 1 and i < N-2:
            A[i,i-2:i+3] = 1
            A[i,i] = 0
        if i == N-2:
            A[i,0] = A[i,N-1] = A[i,N-4:N-2] = 1
        if i == N-1:
            A[i,0:2] = A[i,N-3:N-1] = 1
    #randomly rewiring the network
    for i in range(0,N):
        for j in range(0,N):
            r = rand()
            if r <= P and A[i,j] == 1:
                A[i,j] = 0 #delete old wire
                A[j,i] = 0
                count = 0
                while (count == 0):
                    loc = np.random.randint(N)
                    if A[i,loc] == 0 and loc != i:
                        A[i,loc] = 1
                        A[loc,i] = 1
                        count += 1
    return A

def ScaleFree(N): 
    """
    ScaleFree creates a model Scale-Free network.
    Input:
        N = number of vertices
    Output:
        A = adjacency matrix
    """
    A = np.zeros((N,N))
    #creating a "seed network"
    n=3 #size of the seed network
    A[0:n,0:n] = Random(n,1)
    #adding vertices to the network
    for i in range(n,N):
        for j in range(0,N):
            k_sum = sum(sum(A[0:i,0:N-1]))
            k_j = sum(A[j,:])
            random = rand()
            if k_j and k_sum != 0:
                connection = k_j/k_sum
                if random <= connection:
                    A[i,j] = 1
                    A[j,i] = 1
            A[i,i] = 0
    return A

# =============================================================================
# Network Analysis
# =============================================================================

def Degree(A):
    """
    Degree calculates the degree of each vertex in the network.
    Input:
        A = adjacency matrix
    Output:
        K = array of degree values 
    """
    N=len(A[0,:])
    K = np.zeros(N)
    for i in np.arange(0,N,1):
        K[i] = sum(A[i,:])
    return K

def Cluster(A):
    """
    Cluster calculates the local clustering coefficient of each vertex in the
    network.
    Input:
        A = adjacency matrix
    Output:
        C = array of clustering coefficient values
    """
    N=len(A[0,:])
    K=Degree(A)
    C = np.zeros(N)
    for i in range(0,N):
        count = 0
        for j in range(0,N):
            for k in range(0,N):
                var1 = A[i,j]
                var2 = A[i,k]
                var3 = A[j,k]
                if var1*var2*var3 == 1: 
                    count += 1
        if K[i] > 1:
            # no (k-1)/2 term to account for double counting
            C[i] = (1/(K[i]*(K[i]-1)))*count 
        else:
            C[i] = 0
    return C

def MultiRandom(N,PA): 
    """
    MultiRandom creates multiple random graph networks and calculates the
    degree of each vertex in each network. Note the length of the PA array 
    determines the number of networks created.
    Input:
        N = number of vertices per network
        PA = 1D array of probability of connection values
    Output:
        MR = 2D array of degree values (each column is a network)
    """
    Points = len(PA)
    MR = np.zeros((N,Points))
    for i in range(0,Points):
        A = Random(N,PA[i])
        MR[:,i] = Degree(A)
    return MR 

def Edges(A,K): 
    """
    Edges creates a 2D array of all the edges in a network.
    Input:
        A = adjacency matrix
        K = array of degree values
    Output:
        edges = 2D array of edges (each row is an edge)
    """
    N = len(A[0,:])
    N_edges = int(sum(K)/2)
    edges = np.zeros((int(N_edges),2))
    count1=0 #only count a triangle portion of the matrix
    count2=0 #location in the edges array
    for i in np.arange(0,N,1):
        count1+=1
        for j in np.arange(count1,N,1):
            if A[i,j] == 1:
                edges[count2,0]=i
                edges[count2,1]=j
                count2+=1
    return edges

def TotalEdges(K):
    """
    TotalEdges calculates the total number of edges in a network.
    Input:
        K = array of degree values
    Output:
        m = total number of edges
    """
    m = int(sum(K)/2.0)
    return m

def TotalVertices(A):
    """
    TotalVertices calculates the total number of vertices in a network.
    Input:
        A = adjacency matrix
    Output:
        n = total number of vertices
    """
    n = int(len(A[0,:]))
    return n

def AverageDegree(K):
    """
    AverageDegree calculates the average degree value of a network.
    Input:
        K = array of degree values
    Output:
        z = average degree value
    """
    n = len(K)
    z = sum(K)/n
    return z

def AverageCluster(C):
    """
    AverageCluster calculates the average local clustering coefficient value
    of a network.
    Input:
        C = array of clustering coefficient values
    Output:
        c = average clustering coefficient value
    """
    n = len(C)
    c = sum(C)/n
    return c

# =============================================================================
# Network Plots
# =============================================================================

def Histogram(K): 
    """
    Histogram calculates the frequency of all of the different degree values 
    in a network.
    Input:
        K = array of degree values
    Output:
        histogram = 2D array with degree values in the first column and 
        frequency values in the second column
    """
    max_k = max(K)
    deg = np.linspace(0,max_k,int(max_k+1))
    freq = np.zeros(int(max_k+1))
    for i in range(0,len(K)):
        for j in range(0,len(deg)):    
            if K[i] == deg[j]:
                freq[j]+=1
    N = 0 #count non-zero values
    for i in range(0,len(freq)):
        if freq[i] != 0:
            N += 1
    histogram = np.zeros((N,2))
    count = 0
    for i in range(0,len(freq)):
        if freq[i] != 0:
            histogram[count,1] = freq[i]
            histogram[count,0] = deg[i]
            count+=1
    return histogram

def DVP(Pa,Nt,N,Method): 
    """
    DVP calculates how degree varies with regards to probability in random 
    graph and small world networks. Note that the Pa array corresponds to 
    probability of connection in random graph networks and probability of
    rewiring in small world networks.
    Input:
        Pa = 1D array of probability values
        Nt = number of trials per probability value
        N = number of vertices per network
        Method = either 'Random' or 'SmallWorld'
    Output:
        k_av = 1D array of average degree values of all trial networks at 
        each probability value
        k_minmax = 2D array of the difference between the minimum and maximum 
        average degree values and the average degree value at each probability 
        value (difference to minimum in first row and difference to maximum in
        the second row)
    """
    pts = len(Pa)
    data = np.zeros((Nt,pts))
    for i in range(0,pts):
        for j in range(0,Nt):
            if Method == 'Random':
                A = Random(N,Pa[i])
            if Method == 'SmallWorld':
                A = SmallWorld(N,Pa[i])
            K = Degree(A)
            data[j,i] = AverageDegree(K)
    k_av = np.zeros(pts)
    k_minmax = np.zeros((2,pts))
    for i in range(pts):
        k_av[i] = np.average(data[:,i])
        k_minmax[0,i] = k_av[i]-min(data[:,i])
        k_minmax[1,i] = max(data[:,i])-k_av[i]
    return k_av, k_minmax

def CVP(Pa,Nt,N,Method): 
    """
    CVP calculates how clustering coefficient varies with regards to 
    probability in random graph and small world networks. Note that the Pa 
    array corresponds to probability of connection in random graph networks 
    and probability of rewiring in small world networks.
    Input:
        Pa = 1D array of probability values
        Nt = number of trials per probability value
        N = number of vertices per network
        Method = either 'Random' or 'SmallWorld'
    Output:
        c_av = 1D array of average clustering coefficient values of all trial 
        networks at each probability value
        c_minmax = 2D array of the difference between the minimum and maximum 
        average clustering coefficient values and the average clustering 
        coefficient value at each probability value (difference to minimum in 
        first row and difference to maximum in the second row)
    """
    pts = len(Pa)
    data = np.zeros((Nt,pts))
    for i in range(0,pts):
        for j in range(0,Nt):
            if Method == 'Random':
                A = Random(N,Pa[i])
            if Method == 'SmallWorld':
                A = SmallWorld(N,Pa[i])
            C = Cluster(A)
            data[j,i] = np.average(C)
    c_av = np.zeros(pts)
    c_minmax = np.zeros((2,pts))
    for i in range(pts):
        c_av[i] = np.average(data[:,i])
        c_minmax[0,i] = c_av[i]-min(data[:,i])
        c_minmax[1,i] = max(data[:,i])-c_av[i]
    return c_av, c_minmax

def ScaleFreeDVN(Na,Nt): 
    """
    ScaleFreeDVN calculates how degree varies with regards to number of 
    vertices in scale free networks.
    Input:
        Na = 1D array of number of vertices
        Nt = number of trials per number of vertices value
    Output:
        k_av = 1D array of average degree values of all trial networks at 
        each number of vertices value
        k_minmax = 2D array of the difference between the minimum and maximum 
        average degree values and the average degree value at each number of 
        vertices value (difference to minimum in first row and difference to 
        maximum in the second row)
    """
    pts = len(Na)
    data = np.zeros((Nt,pts))
    for i in range(0,pts):
        for j in range(0,Nt):
            A = ScaleFree(Na[i])
            K = Degree(A)
            data[j,i] = AverageDegree(K)
    k_av = np.zeros(pts)
    k_minmax = np.zeros((2,pts))
    for i in range(pts):
        k_av[i] = np.average(data[:,i])
        k_minmax[0,i] = k_av[i]-min(data[:,i])
        k_minmax[1,i] = max(data[:,i])-k_av[i]
    return k_av, k_minmax

def ScaleFreeCVN(Na,Nt): 
    """
    ScaleFreeCVN calculates how clustering coefficient varies with regards to 
    number of vertices in scale free networks.
    Input:
        Na = 1D array of number of vertices
        Nt = number of trials per number of vertices value
    Output:
        c_av = 1D array of average clustering coefficient values of all trial 
        networks at each number of vertices value
        c_minmax = 2D array of the difference between the minimum and maximum 
        average clustering coefficient values and the average clustering 
        coefficient value at each number of vertices value (difference to 
        minimum in first row and difference to maximum in the second row)
    """
    pts = len(Na)
    data = np.zeros((Nt,pts))
    for i in range(0,pts):
        for j in range(0,Nt):
            A = ScaleFree(int(Na[i]))
            C = Cluster(A)
            data[j,i] = np.average(C)
            print('network %s, trial %s'%(i+1,j+1))
    c_av = np.zeros(pts)
    c_minmax = np.zeros((2,pts))
    for i in range(pts):
        c_av[i] = np.average(data[:,i])
        c_minmax[0,i] = c_av[i]-min(data[:,i])
        c_minmax[1,i] = max(data[:,i])-c_av[i]
    return c_av, c_minmax

def DCVP(Pa,Nt,N,Method): 
    '''
    DCVP calculates how both degree and clustering coefficient vary with 
    regards to probability in random graph networks and small world networks.
    Note that the Pa array corresponds to probability of connection in random 
    graph networks and probability of rewiring in small world networks.
    Input:
        Pa = 1D array of probability values
        Nt = number of trials per probability value
        N = number of vertices per network
        Method = either 'Random' or 'SmallWorld'
    Output:
        k_array = 1D array of average degree values at each probability value
        c_array = 1D array of average clustering coefficient valuess at each
        probability value
    '''
    points = len(Pa)
    c_array = np.zeros(points)
    k_array = np.zeros(points)
    for i in range(0,points):
        c_sum = np.zeros(Nt)
        k_sum = np.zeros(Nt)
        for j in range(0,Nt):
            if Method == 'Random':
                A = Random(N,Pa[i])
            if Method == 'SmallWorld':
                A = SmallWorld(N,Pa[i])
            C = Cluster(A)
            K = Degree(A)
            k_sum[j] = AverageDegree(K)
            c_sum[j] = AverageCluster(C)
            print('DCVP network %s, trial %s'%(i+1,j+1))
        c_array[i] = sum(c_sum)/Nt
        k_array[i] = sum(k_sum)/Nt
    return k_array, c_array

def DCVN(Na,Nt,Method,P): 
    '''
    DCVN calculates how both degree and clustering coefficient vary with 
    regards to number of vertices in random graph, small world and scale-free 
    networks. Note that the P input is only used if the random graph or small
    world methods are chosen, however as value must still be inputted for the
    function to work.
    Input:
        Na = 1D array of number of vertices values
        Nt = number of trials per probability value
        Method = either 'Random', 'SmallWorld' or 'ScaleFree'
        P = probability value, which is only used is the Method input is 
        'Random' or 'SmallWorld'.
    Output:
        k_array = 1D array of average degree values at each number of vertices 
        value
        c_array = 1D array of average clustering coefficient valuess at each
        number of vertices value
    '''
    points = len(Na)
    k_array = np.zeros(points)
    c_array = np.zeros(points)
    for i in range(0,points):
        k_sum = np.zeros(Nt)
        c_sum = np.zeros(Nt)
        for j in range(0,Nt):
            N = int(Na[i])
            if Method == 'Random':
                A = Random(N,P)
            if Method == 'SmallWorld':
                A = SmallWorld(N,P)
            if Method == 'ScaleFree':
                A = ScaleFree(N)
            K = Degree(A)
            C = Cluster(A)
            k_sum[j] = AverageDegree(K)
            c_sum[j] = AverageCluster(C)
            print('DCVN network %s, trial %s'%(i+1,j+1))
        k_array[i] = sum(k_sum)/Nt
        c_array[i] = sum(c_sum)/Nt
    return k_array, c_array