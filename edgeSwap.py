#!/usr/bin/env python3
########################################################################
# File: edgeSwap.py
# Executable: edgeSwap.py
#
# Usage: python3 edgeSwap.py -i 2 -p pathway.sif -u upstream.input -d downstream.input > output.txt
#
# Purpose: Run the disruption function on a null model graph with edges
#           randomly swapped among the source and sink nodes and their
#           connecting linker nodes. 
#          Do this as many times as the user specifies on the command line.
#          Outputs disruption scores for linker nodes only.
#          
#
# Author: Lauren Sanders 
#   Built on work by Marina Haukness in Fall 2016
#   (Based on work by Tulio Torezan Silingardi Del Claro in Summer 2016)
# Date: 12/21/16 
########################################################################

import argparse
import sys
import math
import os
import random
import timeit
import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path

########################################################################
# Command Line
########################################################################

class CommandLine() :
    '''
    to do 
    '''
    def __init__(self) :
        self.parser = argparse.ArgumentParser(description = "This program creates a null model in which upstream-linker and\n"
                                                            "downstream-linker edges have been randomly swapped. Then it runs\n"
                                                            "the disruption function on that null model, as many times as specified\n"
                                                            "in the --iterations option.\n"
                                                            "\n"
                                                            "The pathway interactions file should be a .sif file specified by the --pathway option.\n"
                                                            "The file containing upstream genes and scores should be specified by the --upstream option.\n"
                                                            "The file containing downstream genes and scores should be specified by the --downstream option.\n"
                                                            "Output writes to filename specified by sys.stdout.\n"
                                                            "\n"
                                                            "Output will be same format as the output from the disruption.py program.\n"
                                                            
                                                            "This program has a runtime of ~10 minutes per null model run.\n",
                                              add_help = True, 
                                              usage = "python3 edgeSwap.py -i 2 -p pathway.sif -u upstream.input -d downstream.input > output.txt"
                                                )
       
        self.parser.add_argument('-i', '--iterations', type=int, action='store', help='Number of different null models to create and disrupt')
        self.parser.add_argument('-p', '--pathway', action='store', help='Filename of pathway interactions (.sif) file')
        self.parser.add_argument('-u', '--upstream',action='store', help='Filename of file containing upstream genes and scores')
        self.parser.add_argument('-d', '--downstream', action='store', help='Filename of file containing downstream genes and scores')
        
        self.args = self.parser.parse_args()

########################################################################
# Disruption
########################################################################

def disruption(G, originalMaxFlow, upstreamNodes, downstreamNodes):
    """
    Code by Marina Haukness. 
    
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
    #print("Calculating Disruption... ")
   
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
                print("Error evaluating flow when removing node ", node, file=sys.stderr)
                sys.exit(1)
            G.add_node(node)
            G.add_edges_from(inEdges)
            G.add_edges_from(outEdges)

    nodes = G.nodes()
    sortedDisruption = sorted(nodes, key=disruptionDict.__getitem__, reverse=True)
    sortedDisruption.remove('source')
    sortedDisruption.remove('sink')
    zeroEffectGenes = ""

    #print("Flow Difference {0:2} % Change {1:8} Node type {2:2} Node ".format("","",""), file=sys.stdout)
    #print("Flow Difference" + "\t" + "Change" + "\t" + "Node type" + "Node", file=sys.stdout)
    #print('------------------------------------------------------------------', file=sys.stdout)
    for node in sortedDisruption:
        flowDifference = disruptionDict[node]
        percentChange = flowDifference / originalMaxFlow * 100
        nodeType = "linker"
        if node in upstreamNodes:
            #nodeType = "UPSTREAM" 
            pass
        elif node in downstreamNodes:
            #nodeType = "DOWNSTREAM"
            pass
        if flowDifference <= 0.005:
            pass
            #zeroEffectGenes += node + "\t "
        else:
            print(str(flowDifference) + '\t' + str(percentChange) + '\t' + str(nodeType) + '\t' + str(node), file=sys.stdout)

def removeEmptyEdges(G, flow_dict):
    """ 
    Code by Marina Haukness. 
    
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
# Edge Swap
########################################################################
                
def edgeSwap(filename, G, baseEdgeWeight, upstreamNodes, downstreamNodes, scores):
    """
    Code by Lauren Sanders.
    
    Read in a pathway (.sif) file.

    input: filename, the name of the .sif file
            G, the graph that is being created
            baseEdgeWeight, the weight for the edges that do not have scores
            upstreamNodes, the source nodes
            downstreamNodes, the sink nodes
            scores, a list of all source and sink genes
    output: none, but modifies the graph G
    """
    
    upstreamLinkers = []
    upstreamUsed = []
    downstreamLinkers = []
    downstreamUsed = []
    
    f = open(filename)
    for line in f:
        line = line.rstrip()
        splitLine = line.split('\t')
        gene1 = splitLine[0]
        gene2 = splitLine[2]
        
        # If neither gene is source or sink, 
        # add an edge connecting them normally.
        if not gene1 in scores and not gene2 in scores: 
            G.add_node(gene1)
            G.add_node(gene2)
            G.add_edge(gene1, gene2, capacity=baseEdgeWeight)

        # If either gene is a source or sink, save both genes 
        # for edge-swapping later.
        elif gene1 in upstreamNodes:
            upstreamLinkers.append(gene2)
            upstreamUsed.append(gene1)
        elif gene2 in upstreamNodes:
            upstreamLinkers.append(gene1)
            upstreamUsed.append(gene2)
        elif gene1 in downstreamNodes:
            downstreamLinkers.append(gene2)
            downstreamUsed.append(gene1)
        elif gene2 in downstreamNodes:
            downstreamLinkers.append(gene1)
            downstreamUsed.append(gene2)
    f.close()

    # Swap edges, finish building graph
    for node in upstreamUsed:
        randomUpNode = random.choice(upstreamLinkers)
        G.add_edge(node, randomUpNode, capacity=baseEdgeWeight)
        upstreamLinkers.remove(randomUpNode)
    for node in downstreamUsed:
        randomDownNode = random.choice(downstreamLinkers)
        G.add_edge(randomDownNode, node, capacity=baseEdgeWeight)
        downstreamLinkers.remove(randomDownNode)


########################################################################
# Input File Parser 
########################################################################

def readInputFile(filename, scoreDictionary):
    """
    Code by Marina Haukness. 
    
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
    Code by Marina Haukness and Lauren Sanders.
        
    Run the main program. Parse the input files from the command line,
    create a graph, and calculate the original flow of the graph. 
    
    Create a new graph with randomly swapped edges between linker nodes and 
    upstream and downstream nodes (but not between linker-linker nodes.) Run
    the disruption program through this graph. Repeat this process multiple 
    times according to number of iterations entered on the command line.
        
    """
    # Get variables from user command line input
    cL = CommandLine()
    iterations = cL.args.iterations
    sifFile = cL.args.pathway
    upstreamFile = cL.args.upstream
    downstreamFile = cL.args.downstream
    
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
            
    # Get average score of all upstream / downstream nodes
    # This average will be used as the score for all linker nodes
    scoreList = list(scores.values())
    averageScore = sum(scoreList) / len(scoreList)
    
    # Print header to file 
    print("Flow Difference" + "\t" + "Change" + "\t" + "Node type" + "\t" + "Node", file=sys.stdout)
    
    ###########################################################################
    # Iterate through edge swapping and disruption as many times as specified #
    ###########################################################################
    while True: 
        iterations -= 1
        
        # Read .sif file to add edges between genes that influence each other
        # Swap edges between source/sink genes and their linker nodes 
        # Use base weight for edges that don't have explicit scores from up/downstream files
        # The average score seen is used here as the base, but this could be changed
        edgeSwap(sifFile, G, averageScore, upstreamNodes, downstreamNodes, scores)
        #print("Base edge weight: ", averageScore)
        for (u,v) in G.edges():
            if u in scores and v in scores:
                averagedScore = (scores[u] + scores[v]) / 2.0
                G[u][v]['capacity'] = averagedScore
                
        # Find the maximum flow for the constructed graph
        (flow_value, flow_dict) = nx.maximum_flow(G, 'source', 'sink')
        #print("Flow value: ", "%.2f" % flow_value)
        
        # Remove any edges that didn't have flow going through them
        removeEmptyEdges(G, flow_dict)
        
        # Calculate disruption!
        disruption(G, flow_value, upstreamNodes, downstreamNodes)
        
        for (u,v) in G.edges():
            if G[u][v]['capacity'] == averageScore or G[u][v]['capacity'] == averagedScore:
                G.remove_edge(u,v)
        
        if iterations <= 0 : 
            break
        
if __name__ == "__main__":
    main()
    raise SystemExit
