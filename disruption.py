#!/usr/bin/env python3
########################################################################
# File: disruption.py
# Executable: disruption.py
#
# Usage: python3 disruption.py pathway.sif upstream.input downstream.input 
#               > output.txt
#
# Purpose: This program aims to find the most influential genes in a 
#   network based on how much they disrupt the flow of the network
#   when they are removed.
#
# Author: Marina Haukness
#   Based on work by Tulio Torezan Silingardi Del Claro in Summer 2016
# Date: November 8, 2016
########################################################################

import sys
import math
import os
import random
import timeit
import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path

########################################################################
# Disruption
########################################################################

def disruption(G, originalMaxFlow, upstreamNodes, downstreamNodes):
    """
    Calculate the disruption for each node in the graph.
    This is done by removing the node, calculating the new flow, and
    seeing how different that was from the original.

    input:  G, the graph
            originalMaxFlow, the original maximum flow value
            upstreamNodes, a list of the upstream nodes
            downstreamNodes, a list of the downstream nodes
    output: none, but prints the table of results to stdout 
    """
    #print("Total number of nodes in graph = ", nx.number_of_nodes(G))
    print("Calculating Disruption... ")
   
    disruptionDict = {'source': float('inf'), 'sink': float('inf')}
    sortedGraph = sorted(G.nodes())
    for node in sortedGraph:
        if node != 'source' and node != 'sink':
            outEdges = G.out_edges(node, data=True)
            inEdges = G.in_edges(node, data=True) 
            G.remove_node(node)
            try:
                flow_value = nx.maximum_flow_value(G, 'source', 'sink')
                flowDifference = originalMaxFlow - flow_value
                disruptionDict[node] = flowDifference   
            except:
                print(sys.exc_info()[0])
                print("Error evaluating flow when removing node ", node)
                sys.exit(1)
            G.add_node(node)
            G.add_edges_from(inEdges)
            G.add_edges_from(outEdges)

    nodes = G.nodes()
    sortedDisruption = sorted(nodes, key=disruptionDict.__getitem__, reverse=True)
    sortedDisruption.remove('source')
    sortedDisruption.remove('sink')
    zeroEffectGenes = ""

    print("Flow Difference {0:2} % Change {1:8} Node type {2:2} Node ".format("","",""))
    print('------------------------------------------------------------------')
    for node in sortedDisruption:
        flowDifference = disruptionDict[node]
        percentChange = flowDifference / originalMaxFlow * 100
        nodeType = "linker"
        if node in upstreamNodes:
            nodeType = "UPSTREAM"
        elif node in downstreamNodes:
            nodeType = "DOWNSTREAM"
        if flowDifference <= 0.005:
            zeroEffectGenes += node + "\t "
        else:
            print("{0:20} {1:15} {2:12} {3}".format('%.2f' % flowDifference, 
                '%.2f' % percentChange, nodeType, node))
    # These were clogging up the output, so I condensed them a bit
    print("\nGenes that had zero effect on flow: ")
    print(zeroEffectGenes)


def removeEmptyEdges(G, flow_dict):
    """ 
    Remove all edges from the graph that did not have any flow going through
    them after running a max flow algorithm.

    input: G, the graph
            flow_dict, the dictionary containing the flow through each edge
                resulting from finding the maximum flow
    """
    for node, edges in flow_dict.items():
        for neighbor, flow in edges.items():
            if flow == 0:
                G.remove_edge(node, neighbor)


########################################################################
# Parser
########################################################################
# Note: These functions are based on the assumption that the input files
# are all formatted correctly. A more robust program would do more 
# error checking when reading in these files.

def readSifFile(filename, G, baseEdgeWeight):
    """
    Read in a .sif file.

    input: filename, the name of the .sif file
            G, the graph that is being created
            baseEdgeWeight, the weight for the edges that do not have scores
    output: none, but modifies the graph G
    """
    f = open(filename)
    for line in f:
        line = line.rstrip()
        splitLine = line.split('\t')
        gene1 = splitLine[0]
        gene2 = splitLine[2]
        G.add_node(gene1)
        G.add_node(gene2)
        G.add_edge(gene1, gene2, capacity=baseEdgeWeight)
    f.close()

def readInputFile(filename, scoreDictionary):
    """
    Read in a file containing info about the upstream or downstream genes. 

    input: filename, the name of the file
            scoreDictionary, a dictionary used to keep track of the scores 
    output: a list of all of the genes in the file
            scoreDictionary is also modified.
    """
    f = open(filename)
    genes = []
    for line in f:
        line = line.rstrip()
        splitLine = line.split('\t')
        gene = splitLine[0]
        score = splitLine[1]
        scoreDictionary[gene] =  float(score)
        genes.append(gene)
    f.close()
    return genes

########################################################################
# Main
########################################################################


def main():
    """ 
    Run the main program. Parse the input files from the command line,
    create a graph, calculate the flow of the graph, and calculate the 
    disruption when removing each node in the graph.
    """
    sifFile = sys.argv[1]
    upstreamFile = sys.argv[2]
    downstreamFile = sys.argv[3]

    # Make graph and add the global source and sink
    G = nx.DiGraph()
    G.add_node('source')
    G.add_node('sink')

    # Read in upstream / downstream genes and their scores, add to graph
    scores = {}
    upstreamNodes = readInputFile(upstreamFile, scores)
    downstreamNodes = readInputFile(downstreamFile, scores)
    for gene in upstreamNodes:
        G.add_node(gene)
        G.add_edge('source', gene, capacity=float('inf'))
    for gene in downstreamNodes:
         G.add_node(gene)
         if not G.has_edge('source', gene):
            G.add_edge(gene, 'sink', capacity=float('inf'))
    
    # Read .sif file to add edges between genes that influence each other
    # Use base weight for edges that don't have explicit scores from up/downstream files
    # The average score seen is used here as the base, but this could be changed
    scoreList = list(scores.values())
    averageScore = sum(scoreList) / len(scoreList)
    readSifFile(sifFile, G, averageScore)
    print("Base edge weight: ", averageScore)
    for (u,v) in G.edges():
        if u in scores and v in scores:
            averagedScore = (scores[u] + scores[v]) / 2.0
            G[u][v]['capacity'] = averagedScore
 
    # Find the maximum flow for the constructed graph
    (flow_value, flow_dict) = nx.maximum_flow(G, 'source', 'sink')
    print("Flow value: ", "%.2f" % flow_value)

    # Remove any edges that didn't have flow going through them
    removeEmptyEdges(G, flow_dict)

    # Calculate disruption!
    disruption(G, flow_value, upstreamNodes, downstreamNodes)


if __name__ == "__main__":
    main()
    raise SystemExit
