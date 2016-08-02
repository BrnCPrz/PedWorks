"""
    PedWorks.py | v0.80
    Modeling genealogical structured data with networks theory.

    Input:
        - .ini file: ".txt" file containing setup/options
        - Pedigree: ".txt" file
                    - standard format > animal | sire | dam

    Output:
        - Network/Pedigree statistics
        - Network/Pedigree drawing
        - Pedigree Reordering/Renumbering
        - Community detection and plot

    To be implemented:
        - relationship matrix as Fruchterman & Reingold matrix entry
        - plotting by centrality (highlight n higher valued nodes (size/color))
        - Adjacency matrix plot (improve output)

"""

__author__ = """\n""".join(['Bruno C. Perez (brunocpvet@gmail.com)',
                            'Ricardo V. Ventura (rventura@uoguelph.ca)'])

import ConfigParser as CoP
import os
import sys
import math
from collections import deque
from operator import itemgetter
import itertools

import community #http://perso.crans.org/aynaud/communities/
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import patches


def main():
    """
    In continuous development.
    """
    #print matplotlib.get_backend()

    # PedWork's greeter
    print "\n\n\tStarting \033[0;34mPedWork\033[0;m.py (version - 0.1.0)"
    # print "\n\n\tStarting \033[0;31mPedWork\033[0;m.py (version - 0.1.0)"

    # check command line entries
    print "\n\tChecking command line entries..."
    if len(sys.argv) < 2:
        sys.stderr.write("\t> ERROR: Incorrect entries format!")
        print "\n\tCorrect format: [python] [pedWork.py] [iniFile.ini]\n\n"
        print "\n\tPedWorks.py was built under Python2 language syntax!\n\n"
        sys.stderr.flush()
        sys.exit()
    else:
        print "\tEntries > OK!"

    config = CoP.RawConfigParser()

    # check pedigree file in command line
    print "\n\tChecking", sys.argv[1], "file..."
    if os.path.exists(sys.argv[1]):
        print "\t", sys.argv[1], "> OK!"
        config.read(sys.argv[1])
    else:
        print "\t", sys.argv[1], "> ERROR - File does not exist!"
        print "\tCheck [SETUP].pedfile option in", sys.argv[1]
        sys.exit()


    # get pedfile name
    pedfile = config.get('SETUP', 'pedfile')
    # transform the input pedigree file into a networkX graph object
    PedigNetwork = input_diGraph(pedfile)


    #aMatrix(PedigNetwork)

    # PedigNetwork must be a nx.DiGraph object
    # check if there is any cycles in the pedigree
    cycle_list = list(nx.simple_cycles(PedigNetwork))
    if len(cycle_list) > 0:
        print '\n\tCycles detected for individuals...'
        print '\n\tThe following individuals are involved in cycles:'
        listdc = []
        for i in cycle_list:
            print '\t -', i[0]
            listdc.append(i[0])

    # pedigree cycle detection
    rmCyNode = ''  # declaring (any) value
    if len(cycle_list) > 0:
        print '\tPedigree cycles detected.'
        rmCyNode = raw_input('\n\tRemove individuals to proceed? (y/n): ')
        if rmCyNode == 'y':
            for i in cycle_list:
                PedigNetwork.remove_node(i[0])
            print '\n\tNodes removed with success.'
            print '\n\tProceeding with the algorithm.'
            # creates another file to store the corrected pedigree
            with open(pedfile) as oldfile, open('newpedfile.txt', 'w') as newfile:
                for line in oldfile:
                    elements = line.split()
                    first = elements[0].strip()
                    if first not in listdc and first != 'TP0':
                        newfile.write(line)
        elif rmCyNode == 'n':
            print '\n\tCorrect pedigree errors before proceeding..'
            print '\n\tClosing PedWork.py session...'
            print '\n\tExit > OK!'
            sys.exit()
        else:
            print '\n\tPlease enter y/n.'
            sys.exit()

    # use the corrected file (without cycles)
    if rmCyNode == 'y':
        pedfile = 'newpedfile.txt'

    # grab function [SETUP] from ini. file
    function = config.get('SETUP', 'function')

    # proceed to the function [SETUP] defined
    if function == "reord":
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Pedigree Reordering/Renumbering"
        # gets the [OPTION] values from .ini file
        outfile = config.get('OPTIONS', 'outfile')
        format = config.get('OPTIONS', 'format')
        outheader = config.getboolean('OPTIONS', 'outheader')
        write_reord_ped(pedfile, outfile=outfile, frmt=format, outheader=outheader)
    elif function == "analysis":
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Pedigree Network Analysis"
        cPlot = config.getboolean('OPTIONS', 'ctplot')
        run_analysis(PedigNetwork, cent_plot=cPlot)
    elif function == "draw_simple":
        PedigNetwork = PedigNetwork.to_undirected()
        nscale = config.getfloat('OPTIONS', 'nscale')
        nalpha = config.getfloat('OPTIONS', 'nalpha')
        nsize = config.getfloat('OPTIONS', 'nsize')
        ncolor = config.get('OPTIONS', 'ncolor')
        ealpha = config.getfloat('OPTIONS', 'ealpha')
        ewidth = config.getfloat('OPTIONS', 'ewidth')
        ecolor = config.get('OPTIONS', 'ecolor')
        print "\n\tFUNCTION", ">", function
        print "\tStarting... Simple Pedigree Drawing"
        ped_draw(pedgraph=PedigNetwork, nscale=nscale, nalpha=nalpha, nsize=nsize, ncolor=ncolor,
                 ealpha=ealpha, ewidth=ewidth, ecolor=ecolor)
    elif function == "draw_cluster":
        PedigNetwork = PedigNetwork.to_undirected()
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Cluster Identification Pedigree Drawing"
        nscale = config.getfloat('OPTIONS', 'nscale')
        nalpha = config.getfloat('OPTIONS', 'nalpha')
        nsize = config.getfloat('OPTIONS', 'nsize')
        ealpha = config.getfloat('OPTIONS', 'ealpha')
        ewidth = config.getfloat('OPTIONS', 'ewidth')
        ecolor = config.get('OPTIONS', 'ecolor')
        cSize = config.getint('OPTIONS', 'cSize')
        ped_report(PedigNetwork, ComSize=cSize)
        ped_clus(pedgraph=PedigNetwork, ComSize=cSize, nscale=nscale, nalpha=nalpha, nsize=nsize,
                 ealpha=ealpha, ewidth=ewidth, ecolor=ecolor)
    elif function == "draw_group":
        PedigNetwork = PedigNetwork.to_undirected()
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Group Highlighting Pedigree Drawing"
        nscale = config.getfloat('OPTIONS', 'nscale')
        nalpha = config.getfloat('OPTIONS', 'nalpha')
        nsize = config.getfloat('OPTIONS', 'nsize')
        ealpha = config.getfloat('OPTIONS', 'ealpha')
        ewidth = config.getfloat('OPTIONS', 'ewidth')
        ecolor = config.get('OPTIONS', 'ecolor')
        nlist = config.get('OPTIONS', 'nlist')
        atName = config.get('OPTIONS', 'atName')
        atCol = config.getint('OPTIONS', 'atCol')
        ped_group(pedgraph=PedigNetwork, nlist=nlist, nscale=nscale, nalpha=nalpha, nsize=nsize,
                  ealpha=ealpha, ewidth=ewidth, ecolor=ecolor,
                  atName=atName, atCol=atCol)
    elif function == "draw_multigroup":
        PedigNetwork = PedigNetwork.to_undirected()
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Group Highlighting Pedigree Drawing"
        nscale = config.getfloat('OPTIONS', 'nscale')
        nalpha = config.getfloat('OPTIONS', 'nalpha')
        nsize = config.getfloat('OPTIONS', 'nsize')
        ealpha = config.getfloat('OPTIONS', 'ealpha')
        ewidth = config.getfloat('OPTIONS', 'ewidth')
        ecolor = config.get('OPTIONS', 'ecolor')
        group_list = [e.strip() for e in config.get('OPTIONS', 'group_list').split(',')]
        color_list = [e.strip() for e in config.get('OPTIONS', 'color_list').split(',')]
        draw_multigroup(pedgraph=PedigNetwork, nscale=nscale, nalpha=nalpha, nsize=nsize,
                        ealpha=ealpha, ewidth=ewidth, ecolor=ecolor,
                        group_list=group_list, color_list=color_list)
    elif function == "draw_induced":
        PedigNetwork = PedigNetwork.to_undirected()
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Induced Community Drawing"
        nscale = config.getfloat('OPTIONS', 'nscale')
        nalpha = config.getfloat('OPTIONS', 'nalpha')
        nsize = config.getfloat('OPTIONS', 'nsize')
        ealpha = config.getfloat('OPTIONS', 'ealpha')
        ewidth = config.getfloat('OPTIONS', 'ewidth')
        ecolor = config.get('OPTIONS', 'ecolor')
        draw_induced(pedgraph=PedigNetwork, nscale=nscale, nalpha=nalpha, nsize=nsize,
                     ealpha=ealpha, ewidth=ewidth, ecolor=ecolor)
    elif function == "breed_comp":
        print "\n\tFUNCTION", ">", function
        print "\tStarting...\n\t> Breed Composition Calculation"
        infile = config.get('OPTIONS', 'infile')
        nbreed = config.getint('OPTIONS', 'nbreed')
        breed_comp(pedgraph=PedigNetwork, inFile=infile, nbreed=nbreed)
    else:
        print "\n\t> ERROR:", function, "is not a valid function."
        sys.exit()


def calculate_density(graph):
    print "\n\tCalculating density..."
    g = graph
    dens = nx.density(g)
    print "\t   >  Graph density:", dens
    return g, dens


def calculate_degree(graph):
    print "\n\tCalculating node degree...\n"
    g = graph
    deg = nx.degree(g)
    nx.set_node_attributes(g, 'degree', deg)
    return g, deg


def calculate_indegree(graph):
    # will only work on DiGraph (directed graph)
    print "Calculating Indegree..."
    g = graph
    indeg = g.in_degree()
    nx.set_node_attributes(g, 'indegree', indeg)
    indeg_sorted = sorted(indeg.items(), key=itemgetter(1), reverse=True)
    for key, value in indeg_sorted[0:10]:
        print "   > ", key, value
    return g, indeg


def out_degree_centrality(G):
    """Compute the out-degree centrality for nodes.

    The out-degree centrality for a node v is the fraction of nodes its
    outgoing edges are connected to.

    Parameters
    ----------
    G : graph
        A NetworkX graph

    Returns
    -------
    nodes : dictionary
        Dictionary of nodes with out-degree centrality as values.

    See Also
    --------
    degree_centrality, in_degree_centrality

    Notes
    -----
    The degree centrality values are normalized by dividing by the maximum
    possible degree in a simple graph n-1 where n is the number of nodes in G.

    For multigraphs or graphs with self loops the maximum degree might
    be higher than n-1 and values of degree centrality greater than 1
    are possible.
    """
    if not G.is_directed():
        raise nx.NetworkXError( \
            "out_degree_centrality() not defined for undirected graphs.")
    centrality = {}
    s = 1.0 / (len(G) - 1.0)
    centrality = dict((n, d * s) for n, d in G.out_degree_iter())
    return centrality


def calculate_outdegree(graph):
    # will only work on DiGraph (directed graph)
    print "\tCalculating Outdegree Centrality..."
    g = graph
    outdeg = out_degree_centrality(g)
    nx.set_node_attributes(g, 'outdegree', outdeg)
    outdeg_sorted = sorted(outdeg.items(), key=itemgetter(1), reverse=True)
    for key, value in outdeg_sorted[0:10]:
        print "\t   > ", key, value
    return g, outdeg


def calculate_closeness(graph):
    print "\n\tCalculating Closeness Centrality..."
    g = graph
    clo = nx.closeness_centrality(g)
    nx.set_node_attributes(g, 'closeness', clo)
    degclos_sorted = sorted(clo.items(), key=itemgetter(1), reverse=True)
    for key, value in degclos_sorted[0:10]:
        print "\t   > ", key, value
    return g, clo


def calculate_betweenness(graph):
    print "\n\tCalculating Betweenness Centrality..."
    g = graph
    bc = nx.betweenness_centrality(g.to_undirected())
    nx.set_node_attributes(g, 'betweenness', bc)
    degbetw_sorted = sorted(bc.items(), key=itemgetter(1), reverse=True)
    for key, value in degbetw_sorted[0:10]:
        print "\t   > ", key, value
    return g, bc


def calculate_eigenvector_centrality(graph):
    print "\n\tCalculating Eigenvector Centrality..."
    g = graph
    ec = nx.eigenvector_centrality_numpy(g)
    nx.set_node_attributes(g, 'eigenvector', ec)
    degeign_sorted = sorted(ec.items(), key=itemgetter(1), reverse=True)
    for key, value in degeign_sorted[0:10]:
        print "\t   > ", key, value

    return g, ec


def calculate_katz_centrality(graph):
    """
    Compute the katz centrality for nodes.
    """
    #if not graph.is_directed():
    #    raise nx.NetworkXError( \
    #       "katz_centrality() not defined for undirected graphs.")
    print "\n\tCalculating Katz Centrality..."
    g = graph
    kt = nx.katz_centrality(g)
    nx.set_node_attributes(g, 'katz', kt)
    katz_sorted = sorted(kt.items(), key=itemgetter(1), reverse=True)
    for key, value in katz_sorted[0:10]:
        print "\t   > ", key, value
    return g, kt


def calculate_degree_centrality(graph):
    print "\nCalculating Degree Centrality..."
    g = graph
    dc = nx.degree_centrality(g)
    nx.set_node_attributes(g, 'degree_cent', dc)
    degcent_sorted = sorted(dc.items(), key=itemgetter(1), reverse=True)
    for key, value in degcent_sorted[0:10]:
        print "   > ", key, value

    return graph, dc


def write_node_attributes(graph, filename):
    # utility function to let you print the node + various attributes in a csv format
    for node in graph.nodes(data=True):
        #    print graph.report_node_data(undir_g)
        node_idx, node_dict = node
        attrs = ','.join(str(v) for v in node_dict.values)
        print node  # nx.get_node_attributes(graph, node_idx) #, ",", ",".join(vals)


def write_edge_attributes(graph, filepath, format='networkx', with_data=False):
    """
     Utility function to let you write an edgelist
    """
    print "Writing edgelist to file..."
    if format == 'networkx':
        nx.write_edgelist(graph, filepath, data=with_data)
    else:
        print "generate csv"
    print "Done"


def report_node_data(graph, node=""):
    g = graph
    if len(node) == 0:
        print "Attributes found for the nodes:"
        print g.nodes(data=True)
    else:
        print "Values for node " + node
        print [d for n, d in g.nodes_iter(data=True) if n == node]


def getKey(item):
    return int(item[0])


def run_analysis(dg, cent_plot=False):
    '''
        Include nrank argument!!!
    '''
    print "\n\tBasic graph aspects: \n"
    # on screen basic reports
    print "\t   >  Number of nodes/animals:", len(dg.nodes())
    print "\t   >  Number of connections:", len(dg.edges())
    print "\t   >  Average (out)degree:", (sum(dg.out_degree().values()) / float(len(dg.nodes())))

    #print dg.out_degree().values()

    print "\n\tBuilding Degree Histogram plot."
    degree_histogram(dg)

    # draw_induced(dg)
    dg, dens = calculate_density(dg)
    dg, deg = calculate_degree(dg)

    if type(dg).__name__ == 'DiGraph':
        # dg, indeg = calculate_indegree(dg)
        # dg, indegcent = out_degree_centrality(dg)
        dg, outdeg = calculate_outdegree(dg)
    else:
        dg, degcent = calculate_degree_centrality(dg)

    dg, closn = calculate_closeness(dg)
    dg, bet = calculate_betweenness(dg)
    dg, eigen = calculate_eigenvector_centrality(dg.to_undirected())
    dg, katz = calculate_katz_centrality(dg)

    # create and write pedStats.txt
    with open('pedStats.txt', 'w') as stats:
        for n, d in sorted(dg.nodes_iter(data=True), key=getKey):
            stats.writelines('{:5s} {:4d} {:4f} {:4f} {:4f} {:4f} {:4f} {:4f} \n'.format(str(n), d['degree'],
                                                                                   d['outdegree'],
                                                                                   d['closeness'],
                                                                                   d['betweenness'],
                                                                                   d['outdegree'],
                                                                                   d['eigenvector'],
                                                                                   d['katz']))
    stats.close()
    print "\n\t> Pedigree/Network statistics written to pedStats.txt"

    # routine to print centrality plots
    if cent_plot == True:
        print "\n\t   Preparing to create centrality plots:\n"

        values = [outdeg.get(node) for node in dg.nodes()]
        pos = nx.spring_layout(dg, scale=500)

        nx.draw_networkx_nodes(dg, pos,
                               alpha=1.0, node_color=values, node_size=250, linewidths=0.3,
                               cmap=plt.get_cmap('PuBu'))
        nx.draw_networkx_edges(dg, pos=pos, width=0.3, arrows=False)
        nx.draw_networkx_labels(dg, pos)

        plt.axis("off")
        plt.savefig("outdegPlot.png")
        print "\t   outdegPlot.png > Created"
        plt.show()

        values = [closn.get(node) for node in dg.nodes()]
        nx.draw_networkx_nodes(dg, pos,
                               alpha=1.0, node_color=values, node_size=250, linewidths=0.3,
                               cmap=plt.get_cmap('PuBu'))
        nx.draw_networkx_edges(dg, pos=pos, width=0.3, arrows=False)
        nx.draw_networkx_labels(dg, pos)

        plt.axis("off")
        plt.savefig("closPlot.png")
        print "\t   closPlot.png > Created"
        plt.show()

        values = [bet.get(node) for node in dg.nodes()]
        nx.draw_networkx_nodes(dg, pos,
                               alpha=1.0, node_color=values, node_size=250, linewidths=0.3,
                               cmap=plt.get_cmap('PuBu'))
        nx.draw_networkx_edges(dg, pos=pos, width=0.3, arrows=False)
        nx.draw_networkx_labels(dg, pos)

        plt.axis("off")
        plt.savefig("betwPlot.png")
        print "\t   betwPlot.png > Created"
        plt.show()

        values = [eigen.get(node) for node in dg.nodes()]
        nx.draw_networkx_nodes(dg, pos,
                               alpha=1.0, node_color=values, node_size=250, linewidths=0.3,
                               cmap=plt.get_cmap('PuBu'))
        nx.draw_networkx_edges(dg, pos=pos, width=0.3, arrows=False)
        nx.draw_networkx_labels(dg, pos)

        plt.axis("off")
        plt.savefig("eingenPlot.png")
        plt.show()
        print "\t   eingenPlot.png > Created"

        values = [katz.get(node) for node in dg.nodes()]
        nx.draw_networkx_nodes(dg, pos,
                               alpha=1.0, node_color=values, node_size=250, linewidths=0.3,
                               cmap=plt.get_cmap('PuBu'))
        nx.draw_networkx_edges(dg, pos=pos, width=0.3, arrows=False)
        nx.draw_networkx_labels(dg, pos)

        plt.axis("off")
        plt.savefig("katzplot.png")
        plt.show()
        print "\t   katzPlot.png > Created"

    return dg


def find_partition(graph):
    # code and lib from http://perso.crans.org/aynaud/communities/
    # must be an undirected graph
    g = graph
    partition = community.best_partition(g)
    print "Partitions found: ", len(set(partition.values()))
    # to show members of each partition:
    for i in set(partition.values()):
        members = [nodes for nodes in partition.keys() if partition[nodes] == i]
        print i, len(members)

        # if i==0:
        #      # write out the subgraph
        #      community_graph = graph.subgraph(members)
        #      #draw_graph(community_graph)
        #      #nx.write_edgelist(community_graph, "community.edgelist", data=False)
        #      #for member in members:
        # #         print member, i

        # print "Partition for node johncoogan: ", partition[node_id]
    nx.set_node_attributes(g, 'partition', partition)
    return g, partition


def draw_partition(graph, partition):
    # requires matplotlib.pyplot, uncomment above
    # uses community code and sample from http://perso.crans.org/aynaud/communities/ to draw matplotlib graph in shades of gray
    g = graph
    count = 0
    size = float(len(set(partition.values())))
    pos = nx.spring_layout(g)
    for com in set(partition.values()):
        count = count + 1
        list_nodes = [nodes for nodes in partition.keys()
                      if partition[nodes] == com]
        nx.draw_networkx_nodes(g, pos, list_nodes, node_size=20,
                               node_color=str(count / size))
    nx.draw_networkx_edges(g, pos, alpha=0.5)
    plt.show()


def directed_modularity_matrix(G, nodelist=None):
    """Return the directed modularity matrix of G.
    The modularity matrix is the matrix B = A - <A>, where A is the adjacency
    matrix and <A> is the expected adjacency matrix, assuming that the graph
    is described by the configuration model.
    More specifically, the element B_ij of B is defined as
        B_ij = A_ij - k_i(out) k_j(in)/m
    where k_i(in) is the in degree of node i, and k_j(out) is the out degree
    of node j, with m the number of edges in the graph.
    Parameters
    ----------
    G : DiGraph
       A NetworkX DiGraph
    nodelist : list, optional
       The rows and columns are ordered according to the nodes in nodelist.
       If nodelist is None, then the ordering is produced by G.nodes().
    Returns
    -------
    B : Numpy matrix
      The modularity matrix of G.
    Notes
    -----
    NetworkX defines the element A_ij of the adjacency matrix as 1 if there
    is a link going from node i to node j. Leicht and Newman use the opposite
    definition. This explains the different expression for B_ij.
    See Also
    --------
    to_numpy_matrix
    adjacency_matrix
    laplacian_matrix
    modularity_matrix
    References
    ----------
    .. [1] E. A. Leicht, M. E. J. Newman,
       "Community structure in directed networks",
        Phys. Rev Lett., vol. 100, no. 11, p. 118703, 2008.
    """
    if nodelist is None:
        nodelist = G.nodes()
    A = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, format='csr')
    k_in = A.sum(axis=0)
    k_out = A.sum(axis=1)
    m = G.number_of_edges()
    # Expected adjacency matrix
    X = k_out * k_in / m
    return A - X


def getColumns(inFile, delim="\t", header=True):
    '''
        Get columns of data from inFile. The order of the rows is respected

        :param inFile: column file separated by delim
        :param header: if True the first line will be considered a header line
        :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that
        are headings in the inFile, and values are a list of all the entries in that
        column. indexToName dict maps column index to names that are used as keys in
        the cols dict. The names are the same as the headings used in inFile. If
        header is False, then column indices (starting from 0) are used for the
        heading names (i.e. the keys in the cols dict)
        '''
    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1
    return cols, indexToName


def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[]):
    """
    - G is a networkx graph
    - node_order (optional) is a list of nodes, where each node in G
          appears exactly once
    - partitions is a list of node lists, where each node in G appears
          in exactly one node list
    - colors is a list of strings indicating what color each
          partition should be
    If partitions is specified, the same number of colors needs to be
    specified.
    """
    adjacency_matrix = nx.to_numpy_matrix(G, dtype=np.bool, nodelist=node_order)

    # Plot adjacency matrix in toned-down black and white
    fig = plt.figure(figsize=(8, 8))  # in inches
    plt.imshow(adjacency_matrix,
               cmap="Paired",
               interpolation="none")

    # The rest is just if you have sorted nodes by a partition and want to
    # highlight the module boundaries
    assert len(partitions) == len(colors)
    ax = plt.gca()
    for partition, color in zip(partitions, colors):
        current_idx = 0
        for module in partition:
            ax.add_patch(patches.Rectangle((current_idx, current_idx),
                                           len(module),  # Width
                                           len(module),  # Height
                                           facecolor="none",
                                           edgecolor=color,
                                           linewidth="1"))
            current_idx += len(module)
    plt.show()


def input_ped(inFile, animal=0, sire=1, dam=2, rm=True):
    '''
        inFile - pedigree as .txt file.
        animal - column for the animal ID
        sire - column for the sire ID
        dam - column for the dam ID
        rm - True for removing node = "0"

        > "rm" argument eliminates "0" (missing) individuals from graph
        > if not removed, network topological organization will be impaired

        Output:
            pedForGraph: list of tuples - (animal1/animal2, parent)
            PedGraph: networkX graph
    '''
    ped_file = open(inFile, "r")
    cols, indexToName = getColumns(ped_file, delim=' ', header=False)
    ped_file.close()
    # separate animal / sire / dam in lists
    animal_in = cols[animal]
    sire_in = cols[sire]
    dam_in = cols[dam]

    # creates anim_parent(1/2) tuples
    anim_parent1 = zip(sire_in, animal_in)
    anim_parent2 = zip(dam_in, animal_in)

    # create single object containing all parent_child relationship
    pedForGraph = anim_parent1 + anim_parent2

    PedGraph = nx.from_edgelist(pedForGraph)
    if rm == True:
        PedGraph.remove_node("0")
        return PedGraph

    else:
        return PedGraph


def input_diGraph(inFile, animal=0, sire=1, dam=2, rm=True):
    '''
        inFile - pedigree as .txt file.
        animal - column for the animal ID
        sire - column for the sire ID
        dam - column for the dam ID
        rm - True for removing node = "0"

        > "rm" argument eliminates "0" (missing) individuals from graph
        > if not removed, network topological organization will be impaired

        Output:
            pedForGraph: list of tuples - (animal1/animal2, parent)
            PedGraph: networkX graph
    '''
    ped_file = open(inFile, "r")
    cols, indexToName = getColumns(ped_file, delim=' ', header=False)
    ped_file.close()
    # separate animal / sire / dam in lists
    animal_in = cols[animal]
    sire_in = cols[sire]
    dam_in = cols[dam]

    # creates anim_parent(1/2) tuples
    anim_parent1 = zip(sire_in, animal_in)
    anim_parent2 = zip(dam_in, animal_in)

    # create single object containing all parent_child relationship
    pedForGraph = anim_parent1 + anim_parent2

    PedGraph = nx.DiGraph(pedForGraph)
    if rm == True:
        PedGraph.remove_node("0")
        return PedGraph

    else:
        return PedGraph


def add_node_attribute(inFile, pedgraph, animal=1, atCol=4, atName="attr1"):
    """
    inFile - pedigree as .txt file
    pedgraph - Pedigree as a networkX graph object
    animal - column for the animal ID
    atCol - column for the attribute
    atName - name for the attribute
    """
    ped_df = pd.read_table(inFile, header=None, delim_whitespace=True)
    dic_ped = dict(zip(ped_df[animal - 1], ped_df[atCol - 1]))
    # print dic_ped

    o = nx.set_node_attributes(pedgraph, atName, dic_ped)
    return dic_ped


def ped_draw(pedgraph, nscale=200, nalpha=0.9, nsize=15, ncolor='blue', ealpha=0.3, ewidth=0.5, ecolor="#000000"):
    '''
    Receives a networkx graph and plots.
        - needs matplotlib package

    :param pedgraph: networkX graph object (pedigree)
    :param nscale: scale of the plot
    :param nalpha: node transparency
    :param nsize:  node size
    :param ncolor: node color
    :param ealpha: edge transparency
    :param ewidth: edge width
    :param ecolor: edge color
    :return:
    '''

    pos = nx.spring_layout(pedgraph, scale=nscale)

    # table= [[0 for node in range(0,len(pedgraph.nodes()))] for node in range(0,len(pedgraph.nodes()))]

    # generate D matrix (nodes x nodes)
    Dmatrix = []
    # writes the linear distance between nodes in the graph (varies with POS)
    with open('calculatedD.txt', 'w') as dist:
        for k, v in sorted(pos.iteritems(), key=getKey):
            for l, o in sorted(pos.iteritems(), key=getKey):
                #print k, 'x', l, node_distance(pos[k][0], pos[k][1], pos[l][0], pos[l][1])
                Dmatrix.append([node_distance(pos[k][0], pos[k][1], pos[l][0], pos[l][1])])
                dist.writelines('{:7s} {:7s} {:4f}\n'.format(k, l, node_distance(pos[k][0], pos[k][1], pos[l][0], pos[l][1])))

    #Dmatrixnp = np.array(Dmatrix)
    #Dmatrixnp.shape = (len(pedgraph.nodes()),len(pedgraph.nodes()))


    nx.draw_networkx_nodes(pedgraph, alpha=nalpha, pos=pos,
                           node_color=ncolor, node_size=nsize, linewidths=0.6)

    # label plot not feasible for larger networks/pedigrees
    # nx.draw_networkx_labels(pedgraph, pos)

    nx.draw_networkx_edges(pedgraph, pos, alpha=ealpha, width=ewidth, edge_color=ecolor)

    plt.axis("off")
    plt.show()


def ped_clus(pedgraph, ComSize, nscale=200, nalpha=0.95, nsize=15, ealpha=0.2, ewidth=0.3, ecolor="#000000"):
    '''
    Receives a networkx graph and plots.
        - needs matplotlib package

    :param pedgraph: networkX graph object (pedigree)
    :param nscale: scale of the plot
    :param nalpha: node transparency
    :param nsize:  node size
    :param ncolor: node color
    :param ealpha: edge transparency
    :param ewidth: edge width
    :param ecolor: edge color
    :return:
    '''
    # function level variable definition
    #config = CoP.RawConfigParser()
    #config.read(sys.argv[1])
    #pedfile = config.get('SETUP', 'pedfile')
    #PedigNetwork = input_ped(pedfile)

    part = community.best_partition(pedgraph,resolution=1.0)
    pos = nx.spring_layout(pedgraph, scale=nscale)

    lessOneList = []
    for i in set(part.values()):
            n_member = [nodes for nodes in part.keys() if part[nodes] == i]
            if len(n_member) < ComSize:
                lessOneList.append(n_member)

    #print lessOneList
    adjList = list(itertools.chain(*lessOneList))

    #print adjList
    drawNodes = []
    for i in pedgraph.nodes():
        if i not in adjList:
            drawNodes.append(i)

    #print adjList
    #print len(drawNodes)
    #print len(part)
    for k, v in part.items():
        if k not in drawNodes:
            part.pop(k, None)

    #print len(part)
    values = []
    for i in pedgraph.nodes():
        if i in drawNodes:
            values.append(part.get(i))

    #values = [part.get(k)+1 for k in pedgraph.nodes()]

    nx.draw_networkx_nodes(pedgraph, pos,
                           alpha=nalpha, nodelist=drawNodes, node_color=values, node_size=nsize, linewidths=0.2,
                           cmap=plt.get_cmap('Paired'))

    nx.draw_networkx_nodes(pedgraph, pos,
                           alpha=nalpha, nodelist=adjList, node_color="white", node_size=nsize, linewidths=0.4)

    # label plot not feasible for larger networks
    #nx.draw_networkx_labels(pedgraph, pos)

    nx.draw_networkx_edges(pedgraph, pos, alpha=ealpha, width=ewidth, edge_color=ecolor)
    # ped_report(pedgraph)
    plt.axis("off")
    plt.show()


def ped_report(pedgraph, ComSize):
    # number of detected groups
    part = community.best_partition(pedgraph, resolution=1.0)
    values = [part.get(k)+1 for k in pedgraph.nodes()]
    #print values
    print "\n\tTotal Detected Groups:", max(values)
    # number of groups containing less then 5 individuals
    xlist = []
    for i in set(part.values()):
        n_member = [nodes for nodes in part.keys() if part[nodes] == i]
        if len(n_member) < ComSize:
            xlist.append(i)

    print "\tConflicting Groups Detected: ", len(xlist)
    # number of real communitites detected (>=5 individuals)
    print "\tConsistent Groups Detected: ", max(values) - len(xlist)
    print "\tFinal modularity: ", community.modularity(part, pedgraph)

    # creates Report_01 file containing group_number and number of individuals within group
    with open('Report_01.txt', 'w') as map_report1:
        for i in set(part.values()):
            n_member = [nodes for nodes in part.keys() if part[nodes] == i]
            map_report1.writelines('{:5s} {:4d}\n'.format(str(i + 1), len(n_member)))

    map_report1.close()

    # creates Report_02 file containing animal_id and community_id
    with open('Report_02.txt', 'w') as map_report2:
        map_report2.writelines('{:10s} {:4d}\n'.format(k, v + 1) for k, v in part.items())

    map_report2.close()


def cycle_detect(G):
    # G = nx.DiGraph object
    cycle_list = list(nx.simple_cycles(G))
    print '\n\tCycles detected for individuals: '
    for i in cycle_list:
        print '\t', i[0]

    rmCyNode = raw_input('\n\t> Remove nodes? (y/n): ')
    print rmCyNode
    if rmCyNode == 'y':
        for i in cycle_list:
            G.remove_node(i[0])
        print '\n\tNodes removed with success.'
        print '\n\tProceeding with the algorithm.'
    elif rmCyNode == 'n':
        print '\n\tCorrect pedigree errors before proceeding..'
        print '\n\tClosing PedWork.py session...'
        print '\n\tExit > OK!'
        sys.exit()
    else:
        print '\n\tPlease enter y/n.'
        sys.exit()

    return G


def ped_sort(file):
    """
        - Reorders a pedigree (dict) by the Kahn's Algorithm.

            TO IMPLEMENT:
                - create map
                - get values
                - write .txt reordered and renumbered
        """
    pedgraph = input_diGraph(inFile=file)
    print "\n\tApplying Kahn's Algorithm... "
    in_degree = {u: 0 for u in pedgraph}  # determine in-degree
    for u in pedgraph:  # of each node
        for v in pedgraph[u]:
            in_degree[v] += 1

    Q = deque()  # collect nodes with zero in-degree
    for u in in_degree:
        if in_degree[u] == 0:
            Q.appendleft(u)

    order_list = []  # list for order of nodes

    while Q:
        u = Q.pop()  # choose node of zero in-degree
        order_list.append(u)  # and 'remove' it from graph
        for v in pedgraph[u]:
            in_degree[v] -= 1
            if in_degree[v] == 0:
                Q.appendleft(v)

    if len(order_list) == len(pedgraph):
        return order_list
    else:  # if there is a cycle,
        print "Error: At least one cycle detected!\n"
        return []  # return an empty list


def write_reord_ped(inFile, outfile='output.txt', frmt="fwf", outheader=False):
    """
    :param inFile: pedigree file (animal | sire | dam)
    :param frmt: output format (spc, csv1, csv2, fwf)
    :param outheader: header in the renumbered output
    """
    print '\n\tPreparing ID maps...'
    PedOrig = pd.read_csv(inFile, delim_whitespace=True, header=None)
    PedOrig.columns = ['anID', 'siID', 'daID']
    PedOrig['anID'] = PedOrig['anID'].astype(str)
    PedOrig['siID'] = PedOrig['siID'].astype(str)
    PedOrig['daID'] = PedOrig['daID'].astype(str)

    order_array = np.array(ped_sort(inFile))
    PedData = pd.DataFrame(order_array)
    PedData.columns = ['anID']

    PedMerge01 = pd.merge(PedData, PedOrig, on='anID', how='left')

    if len(PedMerge01['anID'].unique()) < len(PedMerge01['anID']):
        print '\n\tDuplicates were detected!'

    PedMerge02 = PedMerge01.drop_duplicates(subset='anID')
    PedMerge03 = PedMerge02.fillna('0')

    PedData.columns = ['map_id']

    row_count = np.array(range(1, (len(PedData['map_id']) + 1)))

    PedData['renum_id'] = row_count
    PedData['map_id'] = PedData['map_id'].astype(str)

    d = PedData.set_index('map_id').to_dict()

    df3 = PedMerge03.replace(d['renum_id'])
    # print df3
    df3.columns = ['aRenum', 'sRenum', 'dRenum']
    # df3.columns = [col + '/' for col in df3.columns]

    ped1 = pd.concat([PedMerge03, df3], axis=1)
    ped1['aRenum'] = ped1['aRenum'].astype('int')
    ped2 = ped1.sort_values(by=['aRenum'], ascending=[True])

    print '\n\tGenerating', outfile, '...'
    ped3 = ped2[['aRenum', 'sRenum', 'dRenum', 'anID', 'siID', 'daID']]
    # print ped3

    # prepare output file format
    if frmt == "spc":  # simple space delimited
        ped3.to_csv(outfile, index=False, header=outheader, sep=" ")
        print '\t', outfile, 'created with success!'
        return 0
    elif frmt == "csv1":  # comma delimited
        ped3.to_csv(outfile, index=False, header=outheader, sep=",")
        print '\t', outfile, 'created with success!'
        return 0
    elif frmt == "csv2":  # delimited
        ped3.to_csv(outfile, index=False, header=outheader, sep=";")
        print '\t', outfile, 'created with success!'
        return 0
    elif frmt == "fwf":
        with open(outfile, mode='w') as f:
            ped3.to_string(f, header=outheader, index=False)
        print '\t', outfile, 'created with success!'
        return 0
    else:
        print "\tFormat Error!"
        print "\tPlease use a valid 'format' argument -> [spc, csv1, csv2, fwf]"
        return 1


def ped_group(pedgraph, nlist, nscale=600, nalpha=0.95, nsize=15, ealpha=0.2, ewidth=0.3, ecolor="#000000",
              atName="attr1", atCol=4):
    '''
    Receives a networkx graph and plots.
        - needs matplotlib package

    :param pedgraph: networkX graph object (pedigree)
    :param nscale: scale of the plot
    :param nalpha: node transparency
    :param nsize:  node size
    :param ncolor: node color
    :param ealpha: edge transparency
    :param ewidth: edge width
    :param ecolor: edge color
    :return:
    '''
    grpList = [line.strip() for line in open(nlist, 'r')]

    part = add_node_attribute("ped_testherd.in", pedgraph, atName=atName, atCol=atCol)
    values = [part.get(node) for node in grpList]
    pos = nx.spring_layout(pedgraph, scale=nscale)

    nx.draw_networkx_nodes(pedgraph, pos, nodelist=grpList,
                           alpha=nalpha, node_color=values, node_size=nsize, linewidths=0.1,
                           cmap=plt.get_cmap('Paired'))

    # label plot not feasable for larger networks
    # nx.draw_networkx_labels(pedgraph, pos)

    # nx.draw_networkx_edges(pedgraph, pos, alpha=ealpha, width=ewidth, edge_color=ecolor)

    plt.axis("off")
    plt.show()


def draw_multigroup(pedgraph, group_list=[], color_list=[], nscale=200,
                    nalpha=0.95, nsize=15, ealpha=0.2, ewidth=0.3, ecolor="#000000"):
    '''
    Receives a networkx graph and plots.
        - needs matplotlib package

    :param pedgraph: networkX graph object (pedigree)
    :param nscale: scale of the plot
    :param nalpha: node transparency
    :param nsize:  node size
    :param ncolor: node color
    :param ealpha: edge transparency
    :param ewidth: edge width
    :param ecolor: edge color
    :param group_list: list containing file names with the groups
    :return:
    '''

    # list containing the file name for n groups
    grp_list = []
    color = []
    for i in range(len(group_list)):
        grp_list.append([line.strip() for line in open(group_list[i], 'r')])
    for i in range(len(color_list)):
        color.append(color_list[i])

    pos = nx.spring_layout(pedgraph, scale=nscale)

    nx.draw_networkx_nodes(pedgraph, pos, nodelist=None,
                           alpha=nalpha, node_color='#1E90FF', node_size=nsize, linewidths=0.1)

    if len(grp_list) == len(color_list):
        for i in range(0, len(group_list)):
            nx.draw_networkx_nodes(pedgraph, pos, nodelist=grp_list[i],
                                   alpha=nalpha, node_color=color_list[i], node_size=nsize, linewidths=0.1)
    else:
        print "\n\tNumber of groups don't match number of colors!"
        print "\n\tPlease check the .ini file."
        sys.exit()


        # part = add_node_attribute("ped_testherd.in", pedgraph, atName=atName, atCol=atCol)
        # values = [part.get(node) for node in grpList]
        # pos = nx.spring_layout(pedgraph, scale=nscale)

        # nx.draw_networkx_nodes(pedgraph, pos, nodelist=grpList,
        # alpha=nalpha, node_color=values, node_size=nsize, linewidths=0.1,
        # cmap=plt.get_cmap('Paired'))

    # label plot not feasable for larger networks
    # nx.draw_networkx_labels(pedgraph, pos)

    # nx.draw_networkx_edges(pedgraph, pos, alpha=ealpha, width=ewidth, edge_color=ecolor)

    plt.axis("off")
    plt.show()


def degree_histogram(pedgraph):
    """Return a list of the frequency of each degree value.

    Parameters
    ----------
    pedgraph : Networkx graph
       A graph

    Notes
    -----
    Note: the bins are width one, hence len(list) can be large
    (Order(number_of_edges))
    """
    degree_sequence = sorted(nx.degree(pedgraph).values(), reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    dmax = max(degree_sequence)

    plt.loglog(degree_sequence, 'b-', marker='o', markersize=5, markerfacecolor='#FF8C00', antialiased=True,
               color='#000000')
    plt.title("(out)Degree Rank Plot")
    plt.ylabel("(out)Degree")
    plt.xlabel("Rank")

    print "\t   >  Degree histogram plot created in ~/dgRankPlot.png"
    plt.savefig("dgRankPlot.png")
    plt.show()


def ped_dendogram(pedgraph):
    """
    ... under implementation process.
    Now working on:
        - plotting strategy
                > SciPy? MatplotLib?

    """
    dendo = community.generate_dendrogram(pedgraph)
    for level in range(len(dendo) - 1):
        print("partition at level", level,
              "is", community.partition_at_level(dendo, level))


def draw_induced(pedgraph, nscale=10, nalpha=0.95, nsize=80, ealpha=0.2, ewidth=0.3, ecolor="#000000"):
    '''
    Receives a networkx graph and plots.
        - needs matplotlib package

    :param pedgraph: networkX graph object (pedigree)
    :param nscale: scale of the plot
    :param nalpha: node transparency
    :param nsize:  node size
    :param ncolor: node color
    :param ealpha: edge transparency
    :param ewidth: edge width
    :param ecolor: edge color
    :return:
    '''
    # function level variable definition
    config = CoP.RawConfigParser()
    config.read(sys.argv[1])
    pedfile = config.get('SETUP', 'pedfile')
    PedigNetwork = input_ped(pedfile)

    part = community.best_partition(PedigNetwork)
    ind = community.induced_graph(part, PedigNetwork)

    values = []
    for i in ind.nodes():
        values.append(i)

    pos = nx.spring_layout(ind, scale=nscale)

    nx.draw_networkx_nodes(ind, pos,
                           alpha=nalpha, node_color=values, node_size=nsize, linewidths=0.4,
                           cmap=plt.get_cmap('Paired'))

    # label plot not feasable for larger networks
    nx.draw_networkx_labels(ind, pos)

    nx.draw_networkx_edges(ind, pos, alpha=ealpha, width=ewidth, edge_color=ecolor)

    # print list(ind.edges_iter(data=True))
    plt.axis("off")
    plt.savefig("indCommPlot.png")
    plt.show()


def node_distance(xi, yi, xii, yii):
    sq1 = (xi - xii) * (xi - xii)
    sq2 = (yi - yii) * (yi - yii)
    return math.sqrt(sq1 + sq2)


def breed_comp(pedgraph, inFile, nbreed):
    print "\n\tChecking", inFile, "file..."
    if os.path.exists(inFile):
        print "\t", inFile, "> OK!"
    else:
        print "\t", inFile, "> ERROR - File does not exist!"
        print "\tCheck [OPTIONS].infile option"
        sys.exit()

    print "\n\tVerifying breed composition data..."
    with open(inFile, 'r') as document:
        baselist = []
        breedcount = []
        for line in document:
            line = line.split()
            # print line
            if len(line) == 5:
                pedgraph.node[line[0]]['Br1'] = float(line[3])
                pedgraph.node[line[0]]['Br2'] = float(line[4])
                baselist.append(line[0])
                breedcount.append(2)
            elif len(line) == 6:
                pedgraph.node[line[0]]['Br1'] = float(line[3])
                pedgraph.node[line[0]]['Br2'] = float(line[4])
                pedgraph.node[line[0]]['Br3'] = float(line[5])
                baselist.append(line[0])
                breedcount.append(3)
            elif len(line) == 7:
                pedgraph.node[line[0]]['Br1'] = float(line[3])
                pedgraph.node[line[0]]['Br2'] = float(line[4])
                pedgraph.node[line[0]]['Br3'] = float(line[5])
                pedgraph.node[line[0]]['Br4'] = float(line[6])
                baselist.append(line[0])
                breedcount.append(4)
            else:
                pass

    if sum(breedcount) / len(breedcount) == nbreed:
        print "\n\tNumber of Breeds:", nbreed, " > OK!"
        pass
    else:
        print "\n\t[ERROR] Number of breeds vary among individuals"
        print "\n\tPlease check", inFile, "file!"
        sys.exit()

    print "\tbreed_comp data > OK!"
    print "\tbreed_comp baselist > OK"
    # print baselist

    v = []
    for i in sorted(pedgraph.nodes(), key=int):
        v.append(i)

    otherlist = []
    for i in v:
        if i not in baselist:
            otherlist.append(i)
    # print otherlist

    print "\n\tStarting breed_comp calculations."
    for i in sorted(otherlist, key=int):
        b = [int(x) for x, y in pedgraph.in_edges(i)]
        # print i, b
        if not b:
            if nbreed == 2:
                b.append(0)
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
            # print i, b[0], b[1]
            elif nbreed == 3:
                b.append(0)
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
                pedgraph.node[i]['Br3'] = 0
            elif nbreed == 4:
                b.append(0)
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
                pedgraph.node[i]['Br3'] = 0
                pedgraph.node[i]['Br4'] = 0
        elif len(b) == 1:
            if nbreed == 2:
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
            # print i, b[0], b[1]
            elif nbreed == 3:
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
                pedgraph.node[i]['Br3'] = 0
            # print i, b[0], b[1]
            elif nbreed == 4:
                b.append(0)
                pedgraph.node[i]['Br1'] = 0
                pedgraph.node[i]['Br2'] = 0
                pedgraph.node[i]['Br3'] = 0
                pedgraph.node[i]['Br4'] = 0
                # print i, b[0], b[1]
        else:
            if nbreed == 2:
                # print float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br1'] = float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br2'] = float((pedgraph.node[str(b[0])]['Br2'] + pedgraph.node[str(b[1])]['Br2'])) / 2
            # print i, b[0], b[1]
            elif nbreed == 3:
                # print float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br1'] = float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br2'] = float((pedgraph.node[str(b[0])]['Br2'] + pedgraph.node[str(b[1])]['Br2'])) / 2
                pedgraph.node[i]['Br3'] = float((pedgraph.node[str(b[0])]['Br3'] + pedgraph.node[str(b[1])]['Br3'])) / 2
            elif nbreed == 4:
                # print float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br1'] = float((pedgraph.node[str(b[0])]['Br1'] + pedgraph.node[str(b[1])]['Br1'])) / 2
                pedgraph.node[i]['Br2'] = float((pedgraph.node[str(b[0])]['Br2'] + pedgraph.node[str(b[1])]['Br2'])) / 2
                pedgraph.node[i]['Br3'] = float((pedgraph.node[str(b[0])]['Br3'] + pedgraph.node[str(b[1])]['Br3'])) / 2
                pedgraph.node[i]['Br4'] = float((pedgraph.node[str(b[0])]['Br4'] + pedgraph.node[str(b[1])]['Br4'])) / 2
                # print i, b[0], b[1]

    # report_node_data(PedigNetwork)

    with open('breedComp.txt', 'w') as breedComp:
        for n, d in sorted(pedgraph.nodes_iter(data=True), key=getKey):
            if nbreed == 2:
                breedComp.writelines('{:5s} {:4f} {:4f} \n'.format(str(n), d['Br1'], d['Br2']))
            if nbreed == 3:
                breedComp.writelines('{:5s} {:4f} {:4f} {:4f} \n'.format(str(n), d['Br1'], d['Br2'], d['Br3']))
            if nbreed == 4:
                breedComp.writelines('{:5s} {:4f} {:4f} {:4f} {:4f} \n'.format(str(n),
                                                                               d['Br1'],
                                                                               d['Br2'],
                                                                               d['Br3'],
                                                                               d['Br4']))
    breedComp.close()


def aMatrix(pedgraph):
    print "\n\tStarting A-Matrix algorithm."
    sireList = []
    damList = []
    print "\n\t\tObtaining parentage information."
    # ordering nodes (type must be int)
    for i in sorted(pedgraph.nodes(), key=int):
        b = [int(x) for x, y in pedgraph.in_edges(i)]
        if not b:
            b.append(0)
            b.append(0)
            sireList.append(b[0])
            damList.append(b[1])
            # print i, b[0], b[1]
        elif len(b) == 1:
            b.append(0)
            sireList.append(b[0])
            damList.append(b[1])
        else:
           sireList.append(b[0])
           damList.append(b[1])
           # print i, b[0], b[1]

    sireArray = np.asarray(sireList, dtype=int)
    damArray = np.asarray(damList, dtype=int)

    #print sireArray
    #print damArray

    def BuildAMatrix(sirev, damv):
        n = len(sirev)
        N = n + 1
        A = np.zeros((N, N))
        #x = len(sirev[sirev == 0])

        sirev[sirev == 0] = N
        damv[damv == 0] = N
        sirev[sirev !=0] -= 1
        damv[damv !=0] -= 1

        #print sirev
        #print damv

        for a in range(0, n):
            #print n
            A[a, a] = 1 + (A[sirev[a], damv[a]]) / 2
            for j in range(a + 1, n):
                if j > n:
                    pass
                A[a, j] = (A[a, sirev[j]] + A[a, damv[j]]) / 2
                A[j, a] = A[a, j]

        return A[0:n,0:n]

    print "\n\t\tBuilding A-Matrix..."
    A = BuildAMatrix(sireArray, damArray)
    #print A

    print "\n\t\tAttaching kinship coefficients to edge weights"
    for i, u in pedgraph.in_edges():
        #print i,"parent of ", u
        pedgraph.add_edge(str(i), str(u), weight=A[int(i)-1, int(u)-1])

    # generate D matrix (nodes x nodes)
    Amatrix = []
    # writes the linear distance between nodes in the graph (varies with POS)
    with open('calculatedA.txt', 'w') as dist:
        for k in range(0,(len(A))):
            for l in range(0,(len(A))):
                if k==l:
                    Amatrix.append(1.0)
                    dist.writelines('{:7d} {:7d} {:4f}\n'.format(k+1, l+1, 1.0))
                else:
                #print k, 'x', l, node_distance(pos[k][0], pos[k][1], pos[l][0], pos[l][1])
                    Amatrix.append(A[k][l])
                    dist.writelines('{:7d} {:7d} {:4f}\n'.format(k+1, l+1, A[k][l]))

    #for i,u in pedgraph.in_edges():
        #print i, "to", u, pedgraph.get_edge_data(str(i),str(u))
    print "\n\t\tA-Matrix successfully generated."


if __name__ == "__main__":

    main()

    # standard exiting
    print "\n\tPreparing to exit..."
    print "\n\tExit > OK!\n"
    sys.exit(0)
