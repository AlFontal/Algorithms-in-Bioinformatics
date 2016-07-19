#!/usr/bin/env python

"""
Author: Alejandro Fontal
Student nr: 920110242090
Script to: REVEAL algorithm
"""
from __future__ import division
import time
from itertools import combinations
from copy import deepcopy
import math
from operator import mul
eps = 1e-100


def log2(x):
    """
    Executes logarithm of base 2 to the input integer

    :param x: Integer, must be bigger than 0.
    :return: y, the logarithm of base 2 of x.
    """

    y = math.log(x)/math.log(2)

    return y



def entropy(*args):
    """
    Returns the entropy of a single node or the joint entropy of the input args

    :args = list or list of lists containing the on/off states of one or more
    nodes.
    """

    if type(args[0][0]) == list:
        args = args[0]

    array = []
    for arg in args:
        array.append(arg)

    # transform list of states to binary values so as to be able to count them
    # easily for as many arguments as wanted.

    bit_transf = []
    for i in range(len(array)):
        bit_transf.append([pow(2,i)] * len(array[0]))

    newarray = []
    for idx, element in enumerate(array):
        new = map(mul, element, bit_transf[idx])
        newarray.append(new)

    binary = []
    for i in range(len(newarray[0])):
        sum = 0
        for j in range(len(newarray)):
            sum += newarray[j][i]
        binary.append(sum)

    probs = []

    # for each of the possible binary values, count how many times it happens
    # and divide by the total number of occurrences to estimate p(i)

    for pos in range(pow(2,len(newarray))):
        freq = binary.count(pos)
        prob = freq/len(binary)
        probs.append(prob + eps)

    # Calculate the entropy with the given formula.
    entropy = 0
    for p in probs:
        entropy += -p * log2(p)

    return round(entropy,3)


def mutual_information(x, *args):
    """
    Return the mutual information M(x,args).

    x = list containing the on/off states of a certain node.
    args = list or list of lists containing the on/off states of one or more
    nodes.
    """

    tot_list = [x]

    if type(args[0][0]) == list:
        args = args[0]

    for arg in args:

        tot_list.append(arg)
    mut = entropy(x) + entropy(tot_list[1:]) - entropy(tot_list)

    return mut

def reveal(input,output,kmax=3):
    """Returns network, a dictionary containing for every target node
       a list of most likely source nodes. kmax determines how large
       the subsets explored are.
    """
    network = {}
    nodes = sorted(output.keys())  # Sort nodes so they follow same order.
    k = 0
    counter = 0
    while len(nodes) > 0:  # Repeat until all nodes have a solution

        current_node = nodes[0]
        if kmax == 1:
            if counter == 1:
                nodes.remove(current_node)
                if len(nodes) > 0:
                    current_node = nodes[0]
                    counter = 0
            k = 1
            subsets = sorted(input.keys())
            for node in subsets:
                out = output[current_node]
                inp = input[node]
                ent_a = entropy(out)
                ent_b = entropy(inp, out)
                if abs(ent_a - ent_b) < 0.001:
                    nodes.remove(current_node)
                    network[current_node] = node
                    if len(nodes) > 0:
                        current_node = nodes[0]
                        counter = 0
                        break

            counter += 1

        else:

            if k < kmax:
                k += 1
            else:  # If k is kmax then no possible solution can be found with
                #    this kmax for the node, so mark it as unsolvable.
                nodes.remove(current_node)
                current_node = nodes[0]
                k = 1
            subsets = list(combinations(sorted(input.keys()), k))
            for j, subset in enumerate(subsets):
                imp_list = []
                for element in subset:
                    imp_list.append(input[element])

                ent_a = entropy(imp_list)
                imp_list.append(output[current_node])
                ent_b = entropy(imp_list)

                if abs(ent_a - ent_b) < 0.001:
                    network[current_node] = subset
                    nodes.remove(current_node)
                    if len(nodes) > 0:
                        k = 0
                        break

                    else:
                        return network
                elif j == len(subsets)-1:
                    break

    return network

    
def read_tsv_file(tsv):
    """Reads a timeseries dataset and returns the dictionaries
       inputs (time 0,...,T-1) and outputs (time 1,...,T).
    """
    
    header = tsv.readline().rstrip()
    nodes  = header.split("\t")
    
    inputs  = {}; outputs = {}
    for i in range(1,len(nodes)):
      inputs[nodes[i]] = []; outputs[nodes[i]] = []

    line = tsv.readline().rstrip()
    while len(line) > 0:
      vals = line.split("\t")
      for i in range(1,len(nodes)):
        inputs[nodes[i]].append(float(vals[i]))
      line = tsv.readline().rstrip()

    outputs = deepcopy(inputs)
    n = len(inputs[nodes[1]])
    for i in range(1,len(nodes)):
      del inputs[nodes[i]][n-1]; del outputs[nodes[i]][0]
          
    return inputs, outputs

def check_table(input,output):
    """Prints all entropy and mutual information values
       in Figure 5 of Liang et al, PSB 1998 for checking
       your implementation.
    """

    print "H(A), H(B), H(C) =",       entropy(input['A']), entropy(input['B']), entropy(input['C'])
    print "H(A,B), H(B,C), H(A,C) =", entropy(input['A'],input['B']), entropy(input['B'],input['C']), entropy(input['A'],input['C'])
    print "H(A,B,C) =",               entropy(input['A'],input['B'],input['C'])

    print
    print "H(A') =", entropy(output['A'])
    print "H(A',A), H(A',B), H(A',C) =", \
        entropy(output['A'],input['A']), \
        entropy(output['A'],input['B']), \
        entropy(output['A'],input['C'])
    print "M(A',A), M(A',B), M(A',C) =", \
        mutual_information(output['A'],input['A']), \
        mutual_information(output['A'],input['B']), \
        mutual_information(output['A'],input['C'])

    print
    print "H(B') =", entropy(output['B'])
    print "H(B',A), H(B',B), H(B',C) =", \
        entropy(output['B'],input['A']), \
        entropy(output['B'],input['B']), \
        entropy(output['B'],input['C'])
    print "H(B',[A,B]), H(B',[B,C]), H(B',[A,C]) =", \
        entropy(output['B'],input['A'],input['B']), \
        entropy(output['B'],input['B'],input['C']), \
        entropy(output['B'],input['A'],input['C'])
    print "M(B',A), M(B',B), M(B',C) =", \
        mutual_information(output['B'],input['A']), \
        mutual_information(output['B'],input['B']), \
        mutual_information(output['B'],input['C'])
    print "M(B',[A,B]), M(B',[B,C]), M(B',[A,C]) =", \
        mutual_information(output['B'],input['A'],input['B']), \
        mutual_information(output['B'],input['B'],input['C']), \
        mutual_information(output['B'],input['A'],input['C'])

    print
    print "H(C') =", entropy(output['C'])
    print "H(C',A), H(C',B), H(C',C) =", \
        entropy(output['C'],input['A']), \
        entropy(output['C'],input['B']), \
        entropy(output['C'],input['C'])
    print "H(C',[A,B]), H(C',[B,C]), H(C',[A,C]) =", \
        entropy(output['C'],input['A'],input['B']), \
        entropy(output['C'],input['B'],input['C']), \
        entropy(output['C'],input['A'],input['C'])
    print "H(C',[A,B,C]) =", entropy(output['C'],input['A'],input['B'],input['C'])
    print "M(C',A), M(C',B), M(C',C) =", \
        mutual_information(output['C'],input['A']), \
        mutual_information(output['C'],input['B']), \
        mutual_information(output['C'],input['C'])
    print "M(C',[A,B]), M(C',[B,C]), M(C',[A,C]) =", \
        mutual_information(output['C'],input['A'],input['B']), \
        mutual_information(output['C'],input['B'],input['C']), \
        mutual_information(output['C'],input['A'],input['C'])
    print "M(C',[A,B,C]) =", mutual_information(output['C'],input['A'],input['B'],input['C'])
                                                                                                                                                                                                                      
if __name__ == "__main__":
    start_time = time.time()

    a = [0,1,1,1,1,1,1,0,0,0]
    b = [0,0,0,1,1,0,0,1,1,1]
    c = [0,0,1,1,0,0,1,1,0,0]



    # Data used in Liang et al, PSB 1998, in Fig. 1 and on pp. 23:
    
    input = {'A':[0,0,0,0,1,1,1,1],'B':[0,0,1,1,0,0,1,1],'C':[0,1,0,1,0,1,0,1]}
    output = {'A':[0,0,1,1,0,0,1,1],'B':[0,1,0,1,1,1,1,1],'C':[0,0,0,1,0,1,1,1]}

    #QUESTION 1

    print "\nQUESTION 1\n"

    check_table(input,output)
    network = reveal(input,output,3)

    for node in input.keys():
        if node in network.keys():
            print "\nThe expression of node", node,\
                "is best explained by node(s)", network[node]
        else:
            print "\nNode ", node, " seems to be a source node, so it can't" \
                                 " be explained by other nodes. At least not " \
                                 "with kmax nodes only."
    # Data simulated from a real network:

    #QUESTION 2

    input, output = read_tsv_file(open('yeast_bin.tsv'))


    #QUESTION 3

    print "\nQUESTION 3\n"


    for k in range(1, 4):
        print "\nWith kmax = ", k
        network = reveal(input,output,k)
    
        for node in input.keys():
            if node in network.keys():
                print "\nThe expression of node", node,\
                    "is best explained by node(s)", network[node]
            else:
                print "\nNode ", node, "'s behaviour can't be predicted by " \
                                       , k, "other node(s) only."

    #QUESTION 4:


    input1000 = {}
    for key in input.keys():
        input1000[key] = input[key][8999:]

    output1000 = {}
    for key in input.keys():
        output1000[key] = output[key][8999:]

    input100 = {}
    for key in input1000.keys():
        input100[key] = input1000[key][900:]

    output100 = {}
    for key in output1000.keys():
        output100[key] = output1000[key][900:]

    print "\nFor a subset of the last 1000 timepoints:\n"

    network = reveal(input1000,output1000,3)

    for node in input1000.keys():
        if node in network.keys():
            print "\nThe expression of node", node,\
                    "is best explained by node(s)", network[node]
        else:
            print "\nNode ", node, "'s behaviour can't be predicted by " \
                                        "3 other node(s) only."

    print "\nFor a subset of the last 100 timepoints:\n"
    network = reveal(input100,output100,3)

    for node in input100.keys():
        if node in network.keys():
            print "\nThe expression of node", node,\
                "is best explained by node(s)", network[node]
        else:
            print "\nNode ", node, "'s behaviour can't be predicted by " \
                                   , "3 other node(s) only."


    print("--- %s seconds ---" % (time.time() - start_time))

