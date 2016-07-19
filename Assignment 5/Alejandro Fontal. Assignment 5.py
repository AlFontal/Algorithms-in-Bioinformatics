#!/usr/bin/env python

"""
Author:
Student nr:
Script to:
"""

from __future__ import division
from numpy import sqrt, array
import numpy as np
from scipy import sparse
from collections import defaultdict
from heapq import *

def shortest_path(edges, f, t):
    """Return tuple of cost and path for the shortest path from node f to t

    edges: list of edge-tuples [(node1,node2,edge_weight)]
    f: string, label of start node
    t: string, label of end node

    This is a python implementation of Dijkstra's shortest path algorithm.
    """
    g = defaultdict(list)
    for l,r,c in edges:
        g[l].append((c,r))

    q, seen = [(0,f,())], set()
    while q:
        (cost,v1,path) = heappop(q)
        if v1 not in seen:
            seen.add(v1)
            path = (v1, path)
            if v1 == t: 
                return (cost, path)
            for c, v2 in g.get(v1, ()):
                if v2 not in seen:
                    heappush(q, (cost+c, v2, path))

    return float("inf")

class Clusterset(object):
    """Clusterset object describing a cluster

       """
    def __init__(self,left=None,right=None,dissimilarity=0.0,ident=None):
        """ 
        ident: identifier of a leaf node (-1 for internal nodes)
        left: clusterset, left child of the node
        right: clusterset, right child of the node
        dissimilarity: dissimilarity between left and right children
        """
        self.left = left
        self.right = right
        self.ident = ident
        self.dissimilarity = dissimilarity

def print_tree(clust,labels=None,n=0):
    """Print graphical representation of clusterset object
    """
    for i in range(n): print ' ',
    if clust.ident<0: # negative id means that this is branch
        print '-'
    else: # positive id means that this is an endpoint
        if labels==None: 
            print clust.ident
        else: 
            print labels[clust.ident]
    # now print the right and left branches
    if clust.left!=None: 
        print_tree(clust.left,labels=labels,n=n+1)
    if clust.right!=None: 
        print_tree(clust.right,labels=labels,n=n+1)

def get_ordered_elements_dissimilarity(clust,list_of_elements=None,\
    list_of_dissim=None):
    """Return ordered list of elements and dissimilarity from clusterset object
    
    clust: Clusterset object
    list_of_elements: list with node identifiers
    list_of_dissim: list of dissimilarity values 
    """
    if list_of_elements is None: 
        list_of_elements=[]
        list_of_dissim=[]
    if clust.ident < 0: # negative id means that this is a branch
        list_of_elements.append('-')
        # append dissimilarity between the right and left branches
        list_of_dissim.append(clust.dissimilarity)
    else:  # positive id indicates leaf node (the dissim will be zero)
        list_of_elements.append(clust.ident)
        list_of_dissim.append(clust.dissimilarity)
    if clust.left is not None:
        get_ordered_elements_dissimilarity(clust.left,\
            list_of_elements=list_of_elements,list_of_dissim=list_of_dissim)
    if clust.right is not None: 
        get_ordered_elements_dissimilarity(clust.right,\
            list_of_elements=list_of_elements,list_of_dissim=list_of_dissim)

    return list_of_elements, list_of_dissim

def cut_tree(list_of_elements, list_of_dissim, h=1):
    """Return list of clusters by 'cutting' the dendrogram at height h
    
    list_of_elements: list of node identifiers
    list_of_dissim: list of dissimilarity values
    """
    clustlist = []
    cl = []
    for i in range(len(list_of_elements)):
        if list_of_dissim[i] < h:
            if (list_of_elements[i] != '-'):
                cl.append(list_of_elements[i])

        if list_of_dissim[i] >= h:
            if len(cl) > 0:
                clustlist.append(cl)
            cl=[]
    clustlist.append(cl)  
    return clustlist  

def csv_to_matrix(csv):
    """Creates a matrix from a .csv file""" 

    matrix=[]
    for row in csv:
        matrix.append(row.strip().split(','))
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j]= float(matrix[i][j])
    return matrix
 
def csv_to_edges(lines, weight=1):
    """Return list of edges from csv adjacency matrix
    
    lines: open csv file, which behaves like a list of lines
    
    If node1 and node2 are connected in the adjacency matrix,
    they will be added to the list as (node1,node2,1) and 
    (node2,node1,1). In this way the graph is undirected,
    and all edges have weight 1.
    """
    res = []
    node_i = 0
    for line in lines:
        if not line.strip():
            continue
        parts = map(int,line.strip().split(','))
        for node_j, linked in enumerate(parts):
            if linked:
                res.append((node_i, node_j, weight))
        node_i += 1
    return res


def labels_csv_to_dict(filename):
    """

    :param filename: .csv filename
    :return: Dictionary with integers from 0 to n as keys and name of the nodes
    as values
    """

    labels_file = open(filename)
    labels = labels_file.read()
    labels = labels.replace('"','')
    labels = labels.splitlines()
    labels = labels[1:]
    label_list = []
    for node in labels:
        label_list.append(node.split(","))
    label_dict = {}
    for node in label_list:
        label_dict[int(node[0])] = node[1]

    return label_dict

def shortest_dist_from_adjacency(adjacency, edges):
    """

    :param adjacency:
    :param edges:
    :return:
    """
    SD = []
    for i in range(len(adjacency)):
        SD.append([])
        for j in range(len(adjacency)):
            if i == j:
                SD[i].append(0)

            elif i > j:
                SD[i].append(SD[j][i])

            else:
                SD[i].append(shortest_path(edges, i, j)[0])

    return SD


def association_from_shortest_dist(SD):
    """

    :param SD: Matrix of size NxN nodes containing the shortest distance
    between each pair of nodes
    :return: Association matrix of size NxN containing the association values
    for each pair of nodes (1/d^2)
    """
    A = []
    for i in range(len(SD)):
        A.append([])
        for j in range(len(SD)):
            if i == j:
                A[i].append(1)

            elif i > j:
                A[i].append(A[j][i])

            else:
                A[i].append(round(1/(SD[i][j]*SD[i][j]),4))

    return A


def correlation(x, y):
    """

    :param x: list of integers of length n
    :param y: list of integers of length n
    :return: Pearson's correlation between x and y
    """
    n = len(x)
    cov = 0
    sumx = 0
    sumy = 0
    sumxsq = 0
    sumysq = 0

    for i in range(n):
        cov += x[i]*y[i]
        sumx += x[i]
        sumy += y[i]
        sumxsq += x[i] * x[i]
        sumysq += y[i] * y[i]
    num = n*cov - sumx*sumy
    denom = sqrt(n*sumxsq - sumx * sumx) * sqrt(n*sumysq - sumy * sumy)
    corr = num/denom

    return corr


def dissimilarity_matrix_from_correlation(M):
    """

    :param M: Association matrix of size NxN containing the association values
    for each pair of nodes
    :return: Dissimilarity matrix. Contains the dissimilarity (1-Corr) velue
    for each pair of nodes.
    """
    D = []

    for i in range(len(M)):
        D.append([])
        for j in range(len(M)):
            if i == j:
                D[i].append(0)

            elif i > j:
                    D[i].append(D[j][i])

            else:
                    D[i].append(1-correlation(M[i],M[j]))

    return D
 
def hclust(matrix):
    """Return Clusterset object, representing a hierarchical tree
    """

    # The first step transforms the matrix in a numpy array, it will be useful
    # to deal with operations within the matrix.

    matrix = np.array(matrix)

    # initialize the cluster, initially each element is in its own cluster
    clust = [Clusterset(ident=i) for i in range(len(matrix))]

    while len(clust) > 1:  # Loop until a big unique cluster contains all nodes

        # Find the two least dissimilar elements, called A and B
        # To do this, search all the values distinct from 0 in the matrix
        # and return the minimum.
        mindist = np.min(matrix[np.nonzero(matrix)])

        # Given the simetric nature of the matrix, all values will appear at
        # least twice, so just take the first that appears.
        coords = np.argwhere(matrix == mindist)
        AB = coords[0]
        A = AB[0]
        B = AB[1]

        # Delete the columns containing A and B and calculate the avg
        # distance of the two nodes. Generate new row for this new node A-B
        # containing these distances.
        matrix = np.delete(matrix, [A,B], 1)
        newrow = (matrix[A] + matrix[B])/2

        # Now delete the rows containing the A and B nodes, add a new row with
        # the average of the deleted nodes to the end of the matrix.
        matrix = np.delete(matrix, [A,B], 0)
        matrix = np.vstack([matrix, newrow])

        # Add a 0 to the average row to make it possible to add it as the last
        # column of the distance matrix, making it symmetric again.

        newcol = np.append(newrow, 0)
        matrix = np.insert(matrix, len(matrix)-1, [newcol], axis = 1)

        # Now there is a new cluster formed by elements A and that have
        # dissimilarity mindist and an identifier of -1,because it is an
        # internal node
        newcluster = Clusterset(left=clust[A], right=clust[B],\
            dissimilarity=mindist, ident=-1) 
        # now delete the old nodes  from the cluster
        del clust[A]  
        del clust[B-1]
        # add the new cluster
        clust.append(newcluster)
        # If there is still more than one cluster, the loop will start again.

    return clust[0]


 
if __name__ == "__main__": 

    # QUESTION 1:

    print "\nQUESTION 1\n"

    edges = csv_to_edges(open('toy_example.csv'))
    adj = csv_to_matrix(open('toy_example.csv'))
    SD = shortest_dist_from_adjacency(adj, edges)
    print "Distance matrix:\n"
    print('\n'.join([''.join(['{:6}'.format(item) for item in row])for
                     row in SD]))

    ass = association_from_shortest_dist(SD)
    print "\n Association matrix:\n"
    print('\n\n'.join([''.join(['{:10}'.format(item) for item in row])for
                     row in ass]))

    D = dissimilarity_matrix_from_correlation(ass)
    simple_D = []
    print "\n Dissimilarity matrix:\n"
    for row in D:
        rounded_row = [ round(elem, 4) for elem in row ]
        simple_D.append(rounded_row)
    print('\n\n'.join([''.join(['{:10}'.format(item) for item in row])for
                     row in simple_D]))

    # QUESTION 2

    print "\nQUESTION 2\n"


    clust = hclust(D)

    print "The generated tree is: \n"
    print_tree(clust)

    heights = np.arange(0.01,2,0.001)
    len_count = 0
    list_of_elements, list_of_dissim =  get_ordered_elements_dissimilarity(clust)

    for h in heights:
        groups = cut_tree(list_of_elements, list_of_dissim, h=h)
        if len(groups) != len_count:
            print "\nFor height = ", h, " :\n"
            print groups
            len_count = len(groups)

    # QUESTION 3

    print "\nQUESTION 3\n"

    SD = csv_to_matrix(open("shortest_path_distance_high_risk_network.csv"))
    ass = association_from_shortest_dist(SD)
    D = dissimilarity_matrix_from_correlation(ass)
    clust = hclust(D)
    labels = labels_csv_to_dict("metabolite_nodes.csv")
    print "\nThe generated tree for the High Risk Network is: \n"
    print_tree(clust, labels)
    list_of_elements, list_of_dissim =  get_ordered_elements_dissimilarity(clust)



    heights = np.arange(0.1, 1.3, 0.2)
    for h in heights:
        print "\nIf we cut the dendogram at height", h, ":\n"
        k = cut_tree(list_of_elements, list_of_dissim, h)
        labeled_group = []
        for cluster in k:
            nodes = []
            for node in cluster:
                nodes.append(labels[node])
            labeled_group.append(nodes)
        print labeled_group

