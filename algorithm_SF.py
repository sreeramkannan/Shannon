#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Kayvon
#
# Created:     17/09/2013
# Copyright:   (c) Kayvon 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from __future__ import division
from cvxopt import matrix,spmatrix, solvers, spdiag, sparse
import copy,random, numpy
import pdb
import sys, time
from path_decompose_sparse import path_decompose
import os
from time import sleep

sys.setrecursionlimit(100000)
solvers.options['show_progress'] = False
solvers.options['msg_lev'] = 'GLP_MSG_OFF'
worry_abt_unique =0  ## Determines whether the algorithm will make effort to find unique sparsest flow when decomposing a node.
use_norm = 'l1' #Options 'l1' or 'l2'
unit_normalization = False
restored_normalization = True
use_Y_paths = True
use_smoothing = False
use_GLPK = False
path_sparsity = 10

burden_factor = 100

run_penalized = 0

debug_mode = 0
overwrite_normalization = 0 # Set this to 1 if you want to overload the variable normalization to the copy_count of the true
#if overwrite_normalization=1: equivalent to using original Copy counts for thresholding
n_inp = len(sys.argv)
comp = ''
sample_name = ''
sample_name_out = ''
if n_inp>1:
    comp = sys.argv[1]
    if n_inp>2:
        sample_name = sys.argv[2]
        sample_name_out = sample_name
        if n_inp>3:
            sample_name_out = sys.argv[3]

def run_cmd(s1):
        print(s1); os.system(s1)
    
edges_file = sample_name+'intermediate/edges' + comp + '.txt'
nodes_file = sample_name+'intermediate/nodes' + comp + '.txt'
single_nodes_file = sample_name+'intermediate/single_nodes.txt'
KnownPathsFile = sample_name+'intermediate/paths' + comp + '.txt'
reconstr_file = sample_name_out+'algo_output/reconstructed_comp_' +str(comp) + '.fasta'
reconstr_Y_file = sample_name_out+'algo_output/reconstructed' + '.fasta'
use_file =1 #whether you want to read from files or not

if not os.path.exists(sample_name_out+'algo_output'):
    run_cmd('mkdir ' + sample_name_out+'algo_output')



hash_to_node = {}      ## Map: hash value => corresponding node.
node_to_hash = {}      ## Map: node => node's hash value.
known_paths = []       ## A list of known paths.
known_paths_str = []
paths_for_node = {}    ## Map: node => known_paths that the node is involved in.


def single_nodes_to_fasta(): 
    ## The function outputs individual nodes without any edges as reconstructed transcripts.
    sname_head,sname_tail = os.path.split(sample_name)
    try:
            spath,sname = sname_tail.split('__')
    except ValueError:
            sname = sname_tail

    with open(reconstr_file, 'a') as reconstFile:
        i = 0
        for lines in open(single_nodes_file):
            fields = lines.strip().split()
            reconstFile.write('>Shannon_'+sname + '_single_'+str(i)+'\t Copycount:' + fields[2])
            reconstFile.write('\n'+fields[1] +'\n')
            i+=1


def ParseKnownPathsFile(KnownPathsFile, graph): 
    ## This function builds the known_paths and paths_for_node data structures.
    f = open(KnownPathsFile, 'r')
    lines = f.readlines()
    i = 0
    for node1 in graph.nodes:
        paths_for_node[node1] = []
    for (i,line) in enumerate(lines):
        if i != 0:
            tokens = line.split()
            nodes_in_path = []
            tmp_string = ""
            #prev_node = None
            for (j,hashcode) in enumerate(tokens):
                node = hash_to_node[hashcode]
                nodes_in_path.append(node)
                if paths_for_node.get(node) == None:
                    paths_for_node[node]=[i-1]
                else:
                    if len(paths_for_node[node]) < path_sparsity: #Append only if the no of known paths is smaller than path_sparsity
                        paths_for_node[node].append(i-1)
                #prev_node = node
            known_paths.append(nodes_in_path)
    f.close()


# must be called first
def ParseNodeFile(NodeFile, graph):
    ## Builds node_to_hash and hash_to_node.
    ## NodeFile: The file with all the nodes.
    ## graph: the graph object we are using.
    f=open(NodeFile,'r')
    lines = f.readlines()
    i = 0
    for line in lines:
        if i != 0:
            tokens = line.split()
            try:
                t2=float(tokens[2])
            except ValueError:
                t2 = 0
            try:
                t3=float(tokens[3])
            except ValueError:
                t3 = 0
            new_node = Node(tokens[1], t2,t3,tokens[0])
            hash_to_node[tokens[0]] = new_node
            node_to_hash[new_node] = tokens[0]
            graph.add_node(new_node)
        i += 1
    f.close()

def ParseEdgeFile(EdgeFile, graph):
    ## Adds each edge to list of connections for both nodes involved.  *EHC
    ## EdgeFile: The file with all the edges.
    ## graph: the graph object we are using.
    f = open(EdgeFile, 'r')
    lines = f.readlines()
    i = 0
    for line in lines:
        if i != 0:
            tokens = line.split()
            start_node = hash_to_node[tokens[0]]
            end_node = hash_to_node[tokens[1]]
        
            start_node.out_edges.append([end_node, int(tokens[2]), float(tokens[3]), float(tokens[4])])
            end_node.in_edges.append([start_node, int(tokens[2]), float(tokens[3]), float(tokens[4])])
        i += 1
    f.close()


def intersect(a, b, c ):
    ## Returns the intersection of the three lists.  
     return list(set(a) & set(b) & set(c))
     
def intersect5(a, b, c, d, e):
    ## Returns the intersection of the five lists.  
     return list(set(a) & set(b) & set(c) & set(d) & set(e))
     
class Edge(object):  ## Edge object (used for building copycount filter matrices).
    def __init__(self, start_node, end_node, overlap_weight, weight, L):
        self.start = start_node  ## The starting node in the edge.
        self.end = end_node  ## The ending node in the edge.
        self.overlap_weight = overlap_weight  ## How much the connected nodes overlap.
        self.weight = weight  ## The copycount for the edge.
        self.L = L  ## The normalization used for error term in copycount filtering. 

class Node(object):  ## Node object (used universally)
    def __init__(self, node_string, node_weight, L,name):
        self.string = node_string  ## The sequence of bases that the node represents.
        self.in_edges = []  ## A list of edges in which this node object is the end node.
        self.out_edges = []  ## A list of edges in which this node object is the start node.
        self.name = name  ## Hash value for the node.
        self.weight = (node_weight)  ## The copycount for the node.
        self.L = (L)  ## The normalization used for the error term in copycount filtering.  
        self.DNA_start_pos = None  ## Start position in the reference DNA at which this node has bases from 
        self.DNA_end_pos = None  ## Start position in the reference DNA at which this node has bases from.
    def set_string(self, node_string):  ## This is the sequence on the node.
        self.string = node_string
    def add_in_edge(self, node, overlap_weight, weight, L):
        self.in_edges.append([node, overlap_weight, weight, L])
    def add_out_edge(self, node, overlap_weight, weight, L):
        self.out_edges.append([node, overlap_weight, weight, L])

class Graph(object):    ## Graph object (used universally)
    def __init__(self):
        self.start = None    ## Node with no in-edges
        self.end = None      ## Node with no out-edges
        self.nodes = []
        self.tobereduced = []  ## list of nodes with more than one in-edge
        self.paths = []
        #for matrix
        self.edges = []    ## This is used for building the matrix with the node's in/out edge information.
        #for edge filter
        self.edges2 = []
        self.edge_weights = []
        self.node_weights = []
        self.normalization = []
        self.penalization = []
        #for unique solution determination
        self.no_unique_solution = False      ## Is there a unique sparsest flow for each decomposed node in the graph?
        self.paths=[]       ## known paths in graph
        self.paths_Y = []   
        self.og_nodes = {}  ## original nodes in graph before sparse flow is run
        self.constituent_nodes = {}  ## dictionary with nodes as keys and nodes that were condensed to form node as values


    def add_node(self, node):
        self.nodes.append(node)
    def remove_node(self, node):
        self.nodes.remove(node)
    def add_edge(self, start_node, end_node, overlap_weight, weight, L):
        start_node.out_edges.append([end_node, overlap_weight, weight, L])
        end_node.in_edges.append([start_node, overlap_weight, weight, L])

    def findStartAndEnd2(self):
        ## Adds a dummy node labeled start_node that has an out edge to all original nodes with in degree 0.
        ## Adds a dummy node labeled end_node that has an in edge to all original nodes with out degree 0.  
        start_node = Node("Start_", 0, 0,'S')
        end_node = Node("_End", 0, 0,'E')

        for node in self.nodes:
            if len(node.in_edges) == 0:
                node.in_edges.append([start_node, 0, node.weight, 0])
                start_node.out_edges.append([node, 0, node.weight, 0])
                start_node.weight += float(node.weight)

            if len(node.out_edges) == 0:
                node.out_edges.append([end_node, 0, node.weight, 0])
                end_node.in_edges.append([node, 0, node.weight, 0])
                end_node.weight += float(node.weight)

        self.nodes.append(start_node)
        self.nodes.append(end_node)
        self.start = start_node
        self.end = end_node

    def findStartAndEnd3(self):
        '''Adds a dummy node labeled start_node that has an out edge of weight 0 to all original nodes with in degree 0,
        and an out edge of wight -1 to all other origianl nodes.
        Adds a dummy node labeled end_node that has an in edge to all original nodes with out degree 0,
        and an in edge of wight -1 to all other origianl nodes.        
        '''
        start_node = Node("Start_", 0, 0, 'S')
        end_node = Node("_End", 0, 0,'E')

        for node in self.nodes:
            if len(node.in_edges) == 0:
                node.in_edges.append([start_node, 0, node.weight, 0])
                start_node.out_edges.append([node, 0, node.weight, 0])
                start_node.weight += float(node.weight)
            else:
                node.in_edges.append([start_node, 0, node.weight, -1])
                start_node.out_edges.append([node, 0, node.weight, -1])
                start_node.weight += float(node.weight)

            if len(node.out_edges) == 0:
                node.out_edges.append([end_node, 0, node.weight, 0])
                end_node.in_edges.append([node, 0, node.weight, 0])
                end_node.weight += float(node.weight)
            else:
                node.out_edges.append([end_node, 0, node.weight, -1])
                end_node.in_edges.append([node, 0, node.weight, -1])
                end_node.weight += float(node.weight)

        self.nodes.append(start_node)
        self.nodes.append(end_node)
        self.start = start_node
        self.end = end_node



    def printNodes(self):
        ##  Prints out each node with all in-edges and out-edges on the same line.
        print('Nodes:\n' , [[e.name, e.weight, e.L ] for e in self.nodes])
        print('\n')
        for each in self.nodes:
            if len(each.out_edges) == 0:
                list1 = [[each.in_edges[i][0].name,  each.in_edges[i][1], each.in_edges[i][2], each.in_edges[i][3]] for i in range(0, len(each.in_edges))]
                print(each.name,"   out edges: None", "   in edges:",  list1)
            if len(each.in_edges) == 0:
                list1 = [[each.out_edges[i][0].name,  each.out_edges[i][1], each.out_edges[i][2], each.out_edges[i][3]] for i in range(0, len(each.out_edges))]
                print(each.name,"   out edges:", list1, "   in edges: None")
            if len(each.out_edges) != 0 and len(each.in_edges) != 0:
                list_out = [[each.out_edges[i][0].name,  each.out_edges[i][1], each.out_edges[i][2], each.out_edges[i][3]] for i in range(0, len(each.out_edges))]
                list_in = [[each.in_edges[i][0].name,  each.in_edges[i][1], each.in_edges[i][2], each.in_edges[i][3]] for i in range(0, len(each.in_edges))]
                print(each.name,"   out edges:", list_out, "in edges:",list_in)


    def printNodesSmall(self):
        ##  Prints out each node with with it's out edges on the same line.
        for each in self.nodes:
            if len(each.out_edges) != 0:
                list_out = [[each.out_edges[i][0].name] for i in range(0, len(each.out_edges))]
                print(each.name,"   out edges:", list_out)


    def findEdges(self):
        ## finds all edges in the graph and updates all the relevant data structures.
        for node in self.nodes:
            for edge in node.out_edges:
                new_edge = Edge(node, edge[0], edge[1], edge[2], edge[3])
                self.edges.append(new_edge)

                edge_info = (node, edge[0], edge[1], edge[2], edge[3])
                self.edges2.append(edge_info)
                self.edge_weights.append(edge[2])


    def filter_update(self, new_edge_weights):
        # updates all edge weights after filtering
        for node in self.nodes:
            node.out_edges = []
            node.in_edges = [] 
        i = 0
        for edge in self.edges2:
            edge[0].out_edges.append([edge[1], edge[2], new_edge_weights[i], edge[4]])
            edge[1].in_edges.append([edge[0], edge[2], new_edge_weights[i], edge[4]])
            i += 1

    def filter_update_incnodes(self, new_weights,m,n):
        '''Updates all the nodes and edges in the graph with the new copycount information generated from the 
        minimum cost flow filter.  
        new_weights : the new copycount values
        m : the number of edges
        n : the number of nodes
        '''
        i = 0
        for edge in self.edges2:
            if overwrite_normalization:
                edge[0].out_edges.append([edge[1], edge[2], new_weights[i], edge[3]])  #overwrite normalization with true copy-count of the original edge
                edge[1].in_edges.append([edge[0], edge[2], new_weights[i], edge[3]])
            else:
                edge[0].out_edges.append([edge[1], edge[2], new_weights[i], edge[4]])
                edge[1].in_edges.append([edge[0], edge[2], new_weights[i], edge[4]])
            i += 1
        ct = 0
        for node in self.nodes:
            if node is not self.start and node is not self.end:
                node.out_edges = [] 

                node.in_edges = [] 
                node.weight=new_weights[m+ct]
                ct +=1

    def search(self):
        '''Searches for all nodes in the graph that have more than one in edge OR more than one out edge and adds them to
        to the list of nodes to be reduced by path_decompose.  
        If use_Y_paths is on, only sends nodes to path_decompose if the node has one in edge and more than one out edge
        Sorts nodes in topological order.
        '''
        del self.tobereduced[:]
        for node in self.nodes:
            if (len(node.in_edges)!=0)  and (len(node.out_edges)!=0) and (node is not self.start) and (node is not self.end):
                if len(node.in_edges) >1 or len(node.out_edges)>1: #SEARCH FOR any Y nodes
                    if use_Y_paths and len(node.in_edges)<=1:  #Search for left-Y NODES and X-nodes
                        continue
                    self.tobereduced.append(node)
        self.tobereduced.sort(key=lambda x: int(x.name.split("_")[0]), reverse=False)

            
    def algorithm2(self):
        # Runs sparse flow algorithm on graph to simplify graph such that all nodes have an in-degree of 1.
        done = False
        cycle_limit = 1
        in_practice_wau = worry_abt_unique
        algo_iteration= 0
        node_iter = 0 

        for each in self.nodes:
            self.constituent_nodes[each] = [each]
        ## If wau == True, you need to worry about condensing.
        while not done:
            in_practice_wau = worry_abt_unique if algo_iteration<cycle_limit else 0   
            new_node_list = copy.copy(self.nodes) # Algorithm will construct new list of nodes for each iteration through all nodes
            
            for node in self.nodes:
                if node == self.start or node == self.end:
                    continue
                if use_Y_paths or algo_iteration == 0:
                    if len(node.in_edges) <= 1:
                        continue
                else:
                    if len(node.in_edges)<=1 and len(node.out_edges)<=1:
                        continue
                node_iter += 1; 
                sys.stdout.write('\r')
                sys.stdout.write(str(time.asctime())+': comp: ' + str(comp) + ', algo_iter: ' + str(algo_iteration)  + ', node_iter: ' +  str(node_iter) + ', node_name: ' +str(node.name)+ ', m: ' + str(len(node.in_edges)) + ', n: ' + str(len(node.out_edges)) + ', Paths: ' + str(len(paths_for_node.get(node,[]))))
                sys.stdout.flush(); 
                if 1:
                    if 1:
                        new_nodes = []  ## list of new nodes produced from decomposition of the current node.
                        inedges = []  ## A vector that contains the connected node of each in-edge.
                        outedges = []  ## A vector that contains the connected node of each out-edge.
                        inedge_vector = []  ## A vector that contains the copycounts of each in-edge.
                        outedge_vector = []  ## A vector that contains the copycounts of each in-edge.
                        inedge_cc = []  ##  A vector that should contain the copycounts of each in node, but currently contains the overlap of the sequence.
                        outedge_cc = []  ##  A vector that should contain the copycounts of each out node, but currently contains the overlap of the sequence.


                        incoming_edge_attributes = {}  ## A dictionay that contains the overlap and normalization information for each in edge.
                        outgoing_edge_attributes = {}  ## A dictionay that contains the overlap and normalization information for each out edge.
               		if len(node.in_edges) == 0:
			    print('Hanging Node!'); 
                            node.in_edges.append([self.start, 0, node.weight, 0])
                            self.start.out_edges.append([node, 0, node.weight, 0])
                            self.start.weight += float(node.weight)

                        if len(node.out_edges) == 0:
                            print('Hanging Node!'); 
                            node.out_edges.append([self.end, 0, node.weight, 0])
                            self.end.in_edges.append([node, 0, node.weight, 0])
                            self.end.weight += float(node.weight)
 
                        for (i,in_edge) in enumerate(node.in_edges):
                            inedges.append(in_edge[0])
                            inedge_vector.append(float(in_edge[2]))
                            inedge_cc.append(float(in_edge[3]))
                            incoming_edge_attributes[in_edge[0]] = [in_edge[1], in_edge[3]]

                        for (j,out_edge) in enumerate(node.out_edges):
                            outedges.append(out_edge[0])
                            outedge_vector.append(float(out_edge[2]))
                            outedge_cc.append(float(in_edge[3]))
                            outgoing_edge_attributes[out_edge[0]] = [out_edge[1], out_edge[3]]
			
			P = matrix(0.,(len(node.in_edges), len(node.out_edges)))
            
                        #  This section of code determines which known paths will be considered when decomposing this node
                        path_bridge_dict = {}
                        paths_for_all = []        
                        if node in paths_for_node:
                            for kp1 in paths_for_node.get(node):
                                kp1_nodes = known_paths[kp1]
                                if kp1_nodes[0] in self.constituent_nodes[node] or kp1_nodes[-1] in self.constituent_nodes[node]:
                                    paths_for_all.append(kp1)

                        for (m,in_node) in enumerate(inedges):
                            for (n,out_node) in enumerate(outedges):
                                path_bridge_dict[(m, n)] = paths_for_all    

                        if paths_for_node.get(node) != None:
                            for (m,in_node) in enumerate(inedges):
                                if paths_for_node.get(in_node) == None:
                                    continue
                                for (n,out_node) in enumerate(outedges):
                                    if paths_for_node.get(out_node)==None:
                                        continue

                                    node_paths_temp = [paths_for_node[node1] for node1 in self.constituent_nodes[node]]
                                    node_paths = []
                                    for each in node_paths_temp:
                                        node_paths = node_paths + each
                                    if len(node_paths) == 0:
                                        continue
                                    
                                    cand_paths = intersect5(paths_for_node[self.constituent_nodes[in_node][-1]], node_paths, paths_for_node[self.constituent_nodes[out_node][0]], paths_for_node[in_node], paths_for_node[out_node])
                                    l_node = self.constituent_nodes[in_node]
                                    r_node = self.constituent_nodes[out_node]
                                    c_node = self.constituent_nodes[node]
                                    for cp in cand_paths:
                                        node_list = known_paths[cp] #for this path
                                        if self.constituent_nodes[node][0] in node_list and self.constituent_nodes[node][-1] in node_list:
                                            tmp1 = [node_list.index(n1) for n1 in self.constituent_nodes[node]]
                                            node_good = True
                                            for (i, each) in enumerate(tmp1):
                                                if i != 0:
                                                    if prev+1 != each:
                                                        node_good = False
                                                prev = each    
                                            
                                            if node_good == True:
                                                
                                                l_good = True
                                                r_good = True 
                                                l_check  = min(tmp1[0], len(self.constituent_nodes[in_node]))
                                                r_check  = min(len(node_list)-tmp1[-1]-1, len(self.constituent_nodes[out_node]))
                                                for y in range(0, l_check):
                                                    if node_list[tmp1[0]-1-y].string != l_node[-1-y].string:
                                                        l_good = False
                                                for y in range(0, r_check):
                                                    if node_list[tmp1[-1]+1+y].string != r_node[y].string:
                                                        r_good = False
                                                if r_good and l_good:
                                                    P[m, n] = 1 
                                                    path_bridge_dict[(m, n)].append(cp)
                                                
                                        
                        #  This line decomposes the node                
                        output = path_decompose(inedge_vector, outedge_vector, inedge_cc, outedge_cc, overwrite_normalization, P,use_GLPK, path_sparsity)
                        



                        temp_matrix = output[0]
                        m = len(inedge_vector)
                        n = len(outedge_vector)
                        in_node_flow = numpy.sum(temp_matrix, 1)
                        out_node_flow = numpy.sum(temp_matrix, 0)
                        
                        nodes_to_eliminate = [node]
                        
                        # This section of the code builds the new nodes formed during decomposition, and implicitly condenses the 1x1 nodes
                        for i in range(0, m):
                            for j in range(0, n):
                                curr_edge_cc = temp_matrix[i][j]
                                if curr_edge_cc != 0:
                                    in_node_out_deg = len(inedges[i].out_edges)
                                    out_node_in_deg = len(outedges[j].in_edges)
                                    
                                    in_1new = in_node_flow[i] == temp_matrix[i][j]
                                    out_1new = out_node_flow[j] == temp_matrix[i][j]
                                    if in_node_out_deg == 1 and in_1new and out_node_in_deg == 1 and out_1new:
                                        in_node = inedges[i]
                                        out_node = outedges[j]
                                        out_attr = outgoing_edge_attributes[outedges[j]]  
                                        in_attr = incoming_edge_attributes[inedges[i]]    
                                        new_node = Node(in_node.string[:-in_attr[0]] + node.string + out_node.string[out_attr[0]:], curr_edge_cc, node.L,node.name+"_["+str(i)+","+str(j)+"]")
                                        ## Make all new edges for new node.
                                        for in_edge in in_node.in_edges:
                                            new_node.in_edges.append(in_edge)
                                            in_edge[0].add_out_edge(new_node, in_edge[1], in_edge[2], in_edge[3])
                                        ## Make all new edges for new node.
                                        for out_edge in out_node.out_edges:
                                            new_node.out_edges.append(out_edge)
                                            out_edge[0].add_in_edge(new_node, out_edge[1], out_edge[2], out_edge[3])
                                        
                                        self.nodes.append(new_node)
                                        self.constituent_nodes[new_node] = self.constituent_nodes[in_node]+self.constituent_nodes[node]+self.constituent_nodes[out_node]
                                        new_nodes.append(new_node)
                                        new_node_list.append(new_node)
                                        
                                        nodes_to_eliminate.append(in_node)
                                        nodes_to_eliminate.append(out_node)
                                        
                                        paths_for_node[new_node] = path_bridge_dict[(i, j)][:path_sparsity]  #only write out path_sparsity no of paths per node.
                                        
                                    elif in_node_out_deg == 1 and in_1new and not (out_node_in_deg == 1 and out_1new):
                                        in_node = inedges[i]
                                        out_node = outedges[j]
                                        out_attr = outgoing_edge_attributes[outedges[j]]  
                                        in_attr = incoming_edge_attributes[inedges[i]]  
                                        new_node = Node(in_node.string[:-in_attr[0]] + node.string, curr_edge_cc, node.L,node.name+"_["+str(i)+","+str(j)+"]")
                                        ## Make all new edges for new node.
                                        for in_edge in in_node.in_edges:
                                            new_node.in_edges.append(in_edge)
                                            in_edge[0].add_out_edge(new_node, in_edge[1], in_edge[2], in_edge[3])                                       
                                        new_node.add_out_edge(outedges[j], out_attr[0], temp_matrix[i][j], out_attr[1])
                                        outedges[j].add_in_edge(new_node, out_attr[0], temp_matrix[i][j], out_attr[1])
                                        
                                        self.nodes.append(new_node)
                                        self.constituent_nodes[new_node] = self.constituent_nodes[in_node]+self.constituent_nodes[node]
                                        new_nodes.append(new_node)
                                        new_node_list.append(new_node)
                                    
                                        nodes_to_eliminate.append(in_node)
                                        
                                        paths_for_node[new_node] = path_bridge_dict[(i, j)][:path_sparsity]
                                        

                                    elif not (in_node_out_deg == 1 and in_1new) and out_node_in_deg == 1 and out_1new:
                                        in_node = inedges[i]
                                        out_node = outedges[j]
                                        out_attr = outgoing_edge_attributes[outedges[j]]  
                                        in_attr = incoming_edge_attributes[inedges[i]]  
                                        new_node = Node(node.string + out_node.string[out_attr[0]:], curr_edge_cc, node.L,node.name+"_["+str(i)+","+str(j)+"]")                                        
                                        ## Make all new edges for new node.
                                        for out_edge in out_node.out_edges:
                                            new_node.out_edges.append(out_edge)
                                            out_edge[0].add_in_edge(new_node, out_edge[1], out_edge[2], out_edge[3])
                                        new_node.add_in_edge(inedges[i], in_attr[0], temp_matrix[i][j], in_attr[1])
                                        inedges[i].add_out_edge(new_node, in_attr[0], temp_matrix[i][j], in_attr[1])
                                        
                                        
                                        self.nodes.append(new_node)
                                            
                                        self.constituent_nodes[new_node] = self.constituent_nodes[node]+self.constituent_nodes[out_node]
                                        new_nodes.append(new_node)
                                        new_node_list.append(new_node)
                                    
                                        nodes_to_eliminate.append(out_node) 
                                        
                                        paths_for_node[new_node] = path_bridge_dict[(i, j)][:path_sparsity]

                                    else: 
                                        out_attr = outgoing_edge_attributes[outedges[j]]  
                                        in_attr = incoming_edge_attributes[inedges[i]]    

                                        new_node = Node(node.string, temp_matrix[i][j], node.L,node.name+"_["+str(i)+","+str(j)+"]")
                                        new_node.add_in_edge(inedges[i], in_attr[0], temp_matrix[i][j], in_attr[1])
                                        inedges[i].add_out_edge(new_node,in_attr[0], temp_matrix[i][j], in_attr[1])
                                        new_node.add_out_edge(outedges[j], out_attr[0], temp_matrix[i][j], out_attr[1])
                                        outedges[j].add_in_edge(new_node, out_attr[0], temp_matrix[i][j], out_attr[1])
                                        
                                        if len(new_node.out_edges) ==0:
                                            print('woohoo')
                                        self.nodes.append(new_node)

                                        self.constituent_nodes[new_node] = self.constituent_nodes[node]
                                        new_nodes.append(new_node)
                                        new_node_list.append(new_node)
                                        
                                        paths_for_node[new_node] = path_bridge_dict[(i, j)][:path_sparsity]
                        
                        ## For each node that was condensed into a new node, delete all it's connections.
                        for old_node in nodes_to_eliminate:
                            for edge in old_node.in_edges:
                                in_node_temp = edge[0]
                                for oedge in in_node_temp.out_edges:
                                    if oedge[0] is old_node:
                                        in_node_temp.out_edges.remove(oedge)

                            for edge in old_node.out_edges:
                                out_node_temp = edge[0]
                                for iedge in out_node_temp.in_edges:
                                    if iedge[0] is old_node:
                                        out_node_temp.in_edges.remove(iedge)
                            old_node.in_edges = []
                            old_node.out_edges = []
                            if old_node not in self.nodes:
                                'alert'
                            else:    
                                new_node_list.remove(old_node)
            self.nodes = new_node_list # update node list after each iteration through all nodes
            self.search() # checks to see if any more nodes need to be reduced
            if len(self.tobereduced) == 0:
                done = True
            else:
                # This is to ensure nodes are run through topologically
                self.nodes.remove(self.start)
                self.nodes.remove(self.end)
                self.nodes.sort(key=lambda x: int(x.name.split("_")[0]), reverse=False)
                self.nodes.append(self.end)
                self.nodes.insert(0, self.start)
                
            algo_iteration += 1
        sys.stdout.write('\n')               

        
    def read_paths_recursive(self,node,str_till_now,overlap,prev_weight):
        '''Reads all paths in graph recursively
        node:  Current node
        str_till_now:  The string seen before this node.
        overlap:  THe amount of bases of overlap between the last node in the path and the current node.
        prev_weight:  The wieght of thw last node in the path.  
        '''
        curr_str=str_till_now+node.string[overlap:]

        if len(node.out_edges) == 0: ## This assumes all paths end at the _END node.
            if curr_str[-4:] != '_End':
                #Return without appending this path
                return
            else:
                curr_str = curr_str[:-4]    
            self.paths_Y.append([curr_str,prev_weight])
            return
        prev_weight = node.weight
        for (i,each) in enumerate(node.out_edges):
            new_node=each[0]
            overlap = int(each[1])
            self.read_paths_recursive(new_node,curr_str,overlap,prev_weight)



    def read_Y_paths(self):
        ''' Uses read_paths_recursive to find all paths if the graph only has Y nodes 
        (a Y node is a node with at most 1 in edge AND 0 or more out edges).
        '''
	sname_head,sname_tail = os.path.split(sample_name) 
        try:
            spath,sname = sname_tail.split('__')
        except ValueError:
            sname = sname_tail

        with open(reconstr_Y_file, 'a') as pathfile:
            self.search()
            if len(self.tobereduced) != 0:
                print('CAUTION:There are still some unresolved nodes')
            self.read_paths_recursive(self.start,'',0,0)
            for (i,path_str_wt) in enumerate(self.paths_Y):
                path_str = path_str_wt[0][6:]
                path_wt = path_str_wt[1]
		if len(path_str):
                	pathfile.write('>Shannon_'+sname + ' ' +comp+'_'+str(i)+"\t"+str(path_wt))
                	pathfile.write("\n"+path_str+"\n") #with weights




    def read_paths(self):
        ## reads paths in the case when there are only 1x1 nodes (ndoes with at most one in edge AND at most one out edge).
        with open(reconstr_file, 'a') as pathfile:
            self.search()
            if len(self.tobereduced) != 0:
                print('wth')

            
            print("number of paths", len(self.start.out_edges))
            for (i,each) in enumerate(self.start.out_edges):
                string = ''  
                weight = each[2]
                overlap = int(each[1])
                node = each[0]
                seen_nodes = [node]
                while len(node.out_edges) != 0:
                    
                    if len(node.in_edges) > 1:
                        'hi'
                    if node not in self.nodes:
                        print('bad boy')
                        pdb.set_trace()
                    string += node.string[overlap:]


                    if node.out_edges[0][0] not in seen_nodes:
                        overlap = int(node.out_edges[0][1])
                        node = node.out_edges[0][0]
                        seen_nodes.append(node)
                    else:
                        pdb.set_trace()
                        print('BEWARE:there are still some loops')
                        break
                
                string += node.string[overlap:]
                if string[-4:] == '_End':
                    string = string[:-4]    
                self.paths.append([string, weight])
                pathfile.write('>'+sample_name + 'Reconst_'+comp+'_'+str(i)+"\t"+str(weight))
                pathfile.write("\n"+string+"\n") #with weights






def buildMatrixIncNodes(graph):
    '''This builds a matrix that wil be used to constrain the minimum cost flow filtering so that the each node has the 
    amount of flow going in as the amount of flow going out, and the flow value on teh node is eual to the floow of the
    edges leaving the node.  This function also builds the vector that has the normailzation information for each node and edge.  
    graph: the current graph we're filtering.
    '''
    adjacency_matrix = []
    node_indices = {}
    node_indices_other = {}
    edge_indices = {}
    edge_indices_other = {}
    graph.findEdges()
    edge_count = 0
    for edge in graph.edges:
        edge_indices[edge_count] = edge
        edge_indices_other[edge] = edge_count
        edge_count += 1
        if edge.L>=0:
            if unit_normalization:
                graph.normalization.append(1)
            elif restored_normalization:
                graph.normalization.append(1)
            else:
                graph.normalization.append(edge.L)
            graph.penalization.append(0)
        else:
            graph.normalization.append(0)
            graph.penalization.append(1)    

    node_count = 0
    for node in graph.nodes:
        if node is not graph.start and node is not graph.end:
            node_indices[node_count] = node
            node_indices_other[node] = node_count
            graph.node_weights.append(node.weight)
            
            if unit_normalization:
                graph.normalization.append(1)
            elif restored_normalization and node.L == 0:
                graph.normalization.append(len(node.string))
            else:
                graph.normalization.append(node.L)

            graph.penalization.append(0)
            node_count += 1
   
    A = spmatrix(0.,[],[],(2*node_count,edge_count+node_count))  
    for edge in graph.edges:
        if edge.start is not graph.start:
            A[node_indices_other[edge.start], edge_indices_other[edge]] = -1.
        if edge.end is not graph.end:
            A[node_indices_other[edge.end], edge_indices_other[edge]] = 1.
            A[node_count+node_indices_other[edge.end],edge_indices_other[edge]] = -1.
    temp_nc = 0 #temporary node count
    for node in graph.nodes:
        if node is not graph.start and node is not graph.end:
            A[node_count+temp_nc,edge_count+temp_nc]=1.
            temp_nc += 1

    return A



def filter_copycounts_inc_nodes(graph):
    '''This fucntion runs the minimum cost flow algorithm to filter the copycounts of the graph.  
    it has the option of either minimizing the normalized l1 norm of the error between the original copycounts and the 
    new copycounts, or minimizing the normalized l2 norm of the error between the original copycounts and the new copycounts.
    '''
    pen_constant = 10 #set this to 1 so that something like 1/10th of the flow is likely to flow through non-existent edges
    (A) = buildMatrixIncNodes(graph)
    [ta,tb]=A.size
    n = int(ta/2)
    m = int(tb)-n
    I = spmatrix(1.0, range(m+n), range(m+n))
    x = []
    for each in graph.edge_weights:
        x.append(float(each))
    for each in graph.node_weights:
        x.append(float(each))
    x_mat = matrix(x,(m+n,1))
    c = matrix(x,(m+n,1))
    L = matrix(graph.normalization,(m+n,1))
    penality = matrix(graph.penalization,(m+n,1))
    L_th = sum(L)/len(L)*0.001;

    for ctr in range(m+n):
        c[ctr] = x_mat[ctr]*L[ctr]
        if L[ctr]<L_th:
            x_mat[ctr]=0

    pen_cost = 1e10  #set ridiculously large number to force penalization to zero
    if run_penalized:
        q = -c+pen_cost*penality 
    else:
        q = -c
    G = -I

    h = - 0*x_mat # implies f>=0.1c
    dims = {'l': G.size[0], 'q': [], 's': []}
    b = matrix(0.,(2*n,1))
    P = spdiag(graph.normalization)

    #Run it unpenalized in order to calculate the scale for the pen_cost
    if use_norm == 'l2':
        sol=solvers.coneqp(P, q, G, h, dims, A, b)
        x=sol['x']
    elif use_norm == 'l1':
        ## L1 norm cvx_opt 
        L_root = L**(.5)
        c_l1 = matrix([[x_mat*0, L_root]])
        A_l1 = sparse([[A], [A*0]])
        b_l1 = b
        h_l1 = matrix([[h, x_mat, -x_mat]])
        G_l1 = sparse([[G, I, G], [0*I, G, G]])
        print('Generated the matrices, running the solver:')
        if use_GLPK:
		sol = solvers.lp(c_l1, G_l1, h_l1, A_l1, b_l1,solver='glpk')
        else:
		sol = solvers.lp(c_l1, G_l1, h_l1, A_l1, b_l1)

        print('Solver finished')
        x_l1 = sol['x']
        x = x_l1[:m+n, :]

    
    opt_val = sol['primal objective']
    
    
    #Run it penalized to obtain the final answer
    if run_penalized:
        pen_cost = pen_constant*abs(opt_val)/sum(x)  #this is the real value of penality
        q = -c+pen_cost*penality  #check if this is a row vector
        sol=solvers.coneqp(P, q, G, h, dims, A, b)
        x=sol['x']

    ''' Check for negative elements '''    
    i = 0
    for element in x:
        if cmp(element, 0) < 0:
            x[i] = 0.0
        i += 1

    y = numpy.array(x)
    graph.filter_update_incnodes(y,m,n)
    #print(y)
    return x

def buildMatrix(graph):
    # Same as buildMatrix_inc_nodes except doesn't use edge weights
    # Currently not used
    adjacency_matrix = []
    node_indices = {}
    node_indices_other = {}
    edge_indices = {}
    edge_indices_other = {}
    graph.findEdges()
    for (node_count,node) in enumerate(graph.nodes):
        if node is not graph.start and node is not graph.end:
            node_indices[node_count] = node
            node_indices_other[node] = node_count
    for (edge_count,edge) in enumerate(graph.edges):
        edge_indices[edge_count] = edge
        edge_indices_other[edge] = edge_count
    A = matrix(0.,(node_count,edge_count))
    for edge in graph.edges:
        if edge.start is not graph.start:
            A[node_indices_other[edge.start], edge_indices_other[edge]] = -1.
        if edge.end is not graph.end:
            A[node_indices_other[edge.end], edge_indices_other[edge]] = 1.
    return A

def filter_copycounts(graph):
    # Same as filter_copycounts_inc_nodes except doesn't use edge weights
    # Currently not used
    A = buildMatrix(graph)
    [n,m]=A.size
    I = spmatrix(1.0, range(m), range(m))
    c = matrix(map(float,graph.edge_weights),(m,1))
    q = -c  #check if this is a row vector
    G = -I
    h = matrix(0.,(m,1)) # zero matrix
    dims = {'l': G.size[0], 'q': [], 's': []}
    b = matrix(0.,(n,1))
    x=solvers.coneqp(I, q, G, h, dims, A, b)['x']
    y = numpy.array(x)
    graph.filter_update(y)
    return x



# Script to run smoothing, aprse flow algorithm, and output transcripts
# ------------------------------------------------------

    
if comp == '-1':
    single_nodes_to_fasta()
    sys.exit(0)


if use_file:
    graph2 = Graph()
    ParseNodeFile(nodes_file, graph2)
    ParseEdgeFile(edges_file, graph2)
    ParseKnownPathsFile(KnownPathsFile, graph2)
else:
    graph2 = graph1



if run_penalized:
    graph2.findStartAndEnd3()
else:
    graph2.findStartAndEnd2()

if len(graph2.nodes) <= 3:
    if use_Y_paths:
        graph2.read_Y_paths()
    else:
        graph2.read_paths()
    sys.exit(0)
    








if debug_mode:
    graph2.printNodes()
    pdb.set_trace()
    raw_input()
if use_smoothing:
	print('before smoothing')
	new_edge_weights2 = filter_copycounts_inc_nodes(graph2)
	graph2.filter_update(new_edge_weights2)
	print('after smoothing')

if debug_mode:
    graph2.printNodes()
    raw_input()



#DEBUG
for node in graph2.nodes:
    if (node is not graph2.start) and (node is not graph2.end):
        if len(node.out_edges)==0 or len(node.in_edges)==0:
            print('findStartAndEnd2 not working')
            raw_input()




t_start = time.time()
graph2.algorithm2()
t_elapsed = (time.time() - t_start)
#print('after running algorithm' + ' : ' + str(comp) +  " time taken: " + str(t_elapsed) )
#print('after running algorithm')
if debug_mode:
	graph2.printNodes()



if use_Y_paths:
    graph2.read_Y_paths()
else:
    graph2.read_paths()
#print("finished writing file")
#print("No unique solution: " + str(graph2.no_unique_solution)  + ' : ' + str(comp))











                                    
                                    
                                    
                                    
