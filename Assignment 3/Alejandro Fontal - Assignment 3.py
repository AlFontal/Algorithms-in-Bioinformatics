#!/usr/bin/env python

"""
Author: Alejandro Fontal
Student nr: 920110242090
Script to: Assignment 3. Finding Eulerian paths and cycles.
"""

# Implement your functions here


def get_degree(node, graph):
    """
    Returns the degree(inbound edges - outbound edges) of a node in a graph

    :param node: Name of the node contained in a graph.
    :param graph: Graph that can be in form of dictionary or list of tuples.
    :return: degree, integer specifying the difference between the amount of
    inbound edges and outbound edges in a certain node of the graph.
    """
    if type(graph) is dict:

        indegree = [i[j] for i in graph.values() for j in range(len(i))]
        outdegree = len(graph[node])
        degree = indegree.count(node) - outdegree

        return degree

    if type(graph) is list:

        outdegree = [edge[0] for edge in graph]
        indegree = [edge[1] for edge in graph]
        degree = indegree.count(node) - outdegree.count(node)

        return degree

def is_eulerian(graph):
    """
    Evaluates whether a graph is Eulerian or not.

    :param graph: Graph in the form of a dictionary or list of tuples.
    :return: Logical value TRUE/FALSE depending on the Eulerian condition of
    the input graph.
    """

    if type(graph) is dict:
        for node in graph:
            if get_degree(node, graph) != 0:
                return False

    if type(graph) is list:
        nodes = []
        for edge in graph:
            for node in edge:
                if node not in nodes:
                    nodes.append(node)
        for node in nodes:
            if get_degree(node,graph) != 0:
                return False

    return True


def has_eulerian_path(graph, getNodes = False):
    """
    Evaluates whether a graph contains an Eulerian path or not

    :param graph: Graph in the form of a dictionary or list of tuples.
    :return: Logical value TRUE/FALSE depending on whether it is possible to
    obtain an Eulerian path within the graph or not. If getNodes = True, first
    node and last node are also returned. They contain the nodes that start and
    end the eulerian path.
    """
    if type(graph) is dict:
        nodes = graph.keys()
        graph = get_edges(graph)

    if type(graph) is list:
        nodes = []
        for edge in graph:
            for node in edge:
                if node not in nodes:
                    nodes.append(node)

    in_counter = 0
    out_counter = 0

    for node in nodes:
        if get_degree(node, graph) != 0:

            if get_degree(node, graph) == 1:
                in_counter += 1
                last_node = node

            elif get_degree(node, graph) == -1:
                out_counter += 1
                first_node = node

            else:
                print "Graph doesn't contain an Eulerian path"
                return False

    if in_counter < 2 and out_counter < 2:
        if not getNodes:
            return True

        if getNodes:
            return first_node, last_node, True
    else:
        print "Graph doesn't contain an Eulerian path"
        return False


def rotate(cycle, element):
    """
    Rotates a cycle to make a certain element the first and last depicted.

    :param element: Element of the cycle that will end up being the first and
    last of it.
    :param cycle: List of elements starting and ending up in the same one.
    :return: Rotated list. If
    """

    idx = cycle.index(element)
    if idx == 0:
        return cycle
    else:
        return cycle[idx:] + cycle[1:idx] + [element]


def join_cycles(cycle_1, cycle_2, element):
    """
    Takes two cycles and joins them by a common element to make one longer cyc.

    :param cycle_1: List containing a cycle of elements
    :param cycle_2: List containing a cycle of elements
    :param element: Common element in the two lists that will be used to join
    them
    :return: New list containing the two joined cycles in one.
    """
    idx = cycle_1.index(element)
    prev = cycle_1[:idx]
    after = cycle_1[idx + 1:]

    new_cycle = prev + cycle_2 + after

    return new_cycle


def get_edges(graph):
    """
    Takes a dictionary graph and outputs a list of edges as tuples

    :param graph: Graph in the form of a dictionary.
    :return: a sorted list of tuples representing all the edges in a graph. The
    first element of every tuple is the outbound node and the second element is
    the inbound node, so 5 -> 4 is (5, 4)
    """

    edges = []
    for tuple in graph.items():
        for i in range(0, len(tuple[1])):
            edges.append((tuple[0],tuple[1][i]))

    return edges


def find_eulerian_cycle(graph):
    """
    Finds an Eulerian cycle within a graph (cycle visiting all edges)

    :param graph: Graph in the form of a dictionary or list of edges as tuples.
    :return: List of the nodes followed by the cycle.
    """

    if not is_eulerian(graph):
        print "The selected graph is not Eulerian"
        return

    if type(graph) is dict:
        graph = get_edges(graph)

    cycles = {}
    starting_node = graph[0][0]
    current_node = starting_node
    cycles[0] = []
    cycles[0].append(starting_node)
    cyc_count = 0
    step_count = 0

    while len(graph) > 0:   #len(graph) > 0 until all edges have been visited.

        for edge in graph:
            if current_node == edge[0]:
                current_node = edge[1]
                graph.remove(edge)
                cycles[cyc_count].append(current_node)
                step_count += 1

        if current_node == starting_node and step_count > 0 and len(graph) > 0:
        #If this happens, a cycle has been completed, so get a new starting
        #node, reset the steps counter and start searching for next cycle
                cyc_count += 1
                current_node = graph[0][0]
                step_count = 0
                starting_node = current_node
                cycles[cyc_count] = []
                cycles[cyc_count].append(current_node)



    cycles_list = []
    for i in range(len(cycles)):
        cycles_list.append(cycles[i])

    while len(cycles_list) > 1:

        for idx, cycle in enumerate(cycles_list):
            if idx < len(cycles_list) - 1:

                for i, node in enumerate(cycle):
                    for j in range(1, len(cycles_list) - idx):
                        if node in cycles_list[idx + j]:
                            cycles_list[idx + j] = rotate(cycles_list[idx + j],
                                                          node)
                            cycles_list[idx] = join_cycles(cycles_list[idx],
                                                           cycles_list[idx + j],
                                                           node)
                            del cycles_list[idx + j]

    cycles_list = cycles_list[0]

    nodes_list = []

    for element in cycles_list:
        if type(element) is list:
            for node in element:
                nodes_list.append(node)
        else:
            nodes_list.append(element)

    return nodes_list


def find_eulerian_path(graph):
    """

    :param graph:
    :return:
    """
    if is_eulerian(graph):
        cycle = find_eulerian_cycle(graph)
        path = cycle[1:]
        return path

    elif not has_eulerian_path(graph):
        print "The graph doesn't contain any Eulerian path"
        return

    else:
        if type(graph) == dict:
            graph = get_edges(graph)

        first_node, last_node, has_path = has_eulerian_path(graph, getNodes = True)

        new_edge = (last_node, first_node)
        graph.append(new_edge)

        cycle = find_eulerian_cycle(graph)
        cycle = rotate(cycle, last_node)
        path = cycle[1:]

    return path


def spectrum_to_graph(spectrum, edges = False):
    """

    :param spectrum: list of lmers
    :return:
    """

    graph = {}
    edges = []
    for lmer in spectrum:
        edges.append((lmer[0:len(lmer)-1], lmer[1:len(lmer)]))

    for edge in edges:
        if edge[0] in graph:
            graph[edge[0]].append(edge[1])
        else:
            graph[edge[0]] = [edge[1]]

    if edges == True:
        return edges, graph

    else:
        return graph

def path_to_seq(path):

    seq = []
    for idx,node in enumerate(path):
        if idx == 0:
            seq.append(node)
        else:
            seq.append(node[-1])

    seq = "".join(seq)
    return seq




if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {1:[3],3:[7],10:[9],9:[4],4:[2,5],7:[10,8],2:[1],5:[8],\
        8:[4,6],6:[7]}
    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {'A':['C'],'C':['E'],'F':['G'],'G':['H'],'H':['A','J'],\
        'E':['F','K'],'J':['K'],'K':['H','D'],'D':['E','I'],\
        'I':['B'], 'B':['D']}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 


# Question 1
print "\nQUESTION 1 \n"
if is_eulerian(graph_822):
    print "The graph is Eulerian since all of its vertices are balanced"
else:
    print "The graph is not Eulerian since there is at least one unbalanced" \
          " vertex"



# Question 2
print "\nQUESTION 2 \n"

if has_eulerian_path(graph_822):
    print "The graph contains an Eulerian path because it has, if any, 2 " \
          "unbalanced vertices at most with degree 1 or -1"
else:
    print "The graph doesn't contain an Eulerian path"

# Question 3
print "\nQUESTION 3 \n"

cycle = find_eulerian_cycle(graph_822)
print "The Eulerian cycle found by my algorithm for graph_822 is: \n"
print cycle

# Question 4
print "\nQUESTION 4 \n"

print "First attempt:"
cycle = find_eulerian_cycle(graph_822)
print cycle

print "\nSecond attempt:"
cycle = find_eulerian_cycle(graph_822)
print cycle
print "\nThird attempt:"
cycle = find_eulerian_cycle(graph_822)
print cycle

# Question 5
print "\nQUESTION 5 \n"
graph = spectrum_to_graph(s)
print "The graph obtained from the spectrum s is:"

for k,v in graph.items():
    print k, v

# Question 6
print "\nQUESTION 6 \n"
print is_eulerian(graph)
print has_eulerian_path(graph)

# Question 7
print "\nQUESTION 7 \n"
path = find_eulerian_path(graph)
print "The path found for spectrum s is:"
print path
seq = path_to_seq(path)
print "\nThe DNA sequence to which that path corresponds is:"
print seq
# Question 8
print "\nQUESTION 8 \n"
if is_eulerian(bigger_graph):
    cycle = find_eulerian_cycle(bigger_graph)
    print "bigger_graph is Eulerian so an Eulerian cycle could be found:"
    print cycle
else:
    if has_eulerian_path(bigger_graph):
        path = find_eulerian_path(bigger_graph)
        print "bigger_graph isn't Eulerian but it contains an Eulerian path:"
        print path

