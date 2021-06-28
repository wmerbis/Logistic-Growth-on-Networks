import networkx as nx
import networkx.algorithms.isomorphism as iso 
import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np
import sympy as sym
from sympy import Sum, factorial, binomial, oo, IndexedBase, Function,symbols
from sympy.functions import exp
from sympy import sympify
import time
import csv
from sympy.matrices import Matrix


t, s = symbols('t s', real=True)

x, y, z = symbols('x y z', real=True)

def coloring(g) :
    '''
    Creates list of node colors, depending on whether the nodes of the network 'g' are in the states "X", "I" or "B"
    '''
    nodecolor=[]
    for i in g.nodes:
        if g.nodes[i]['state'] == "B" :
            nodecolor.append("deepskyblue")
        elif g.nodes[i]['state'] == 'I' :
            nodecolor.append("tomato")
        elif g.nodes[i]['state'] == 'X':
            nodecolor.append("darkred")
        else :
            nodecolor.append("yellow")
    
    return nodecolor


def draw_diagrams(G):
    '''Draws a diagram of the graph G displaying initial infected (sources) as dark red nodes, infected as red nodes and sinks 
    as blue nodes. '''
    plt.figure(figsize=(2, 2))
    nx.draw(G, node_color=coloring(G), with_labels=False)
    #plt.title("sign = "+str(G.graph['coeff']))
    plt.show()


def are_the_same(g,h):
    '''Checks whether two graphs g and h correspond to equivalent diagrams, taking into account the nodes states'''
    ndm=iso.categorical_node_match('state',"S") 
    return nx.is_isomorphic(g,h,node_match=ndm)


def plot_from_ls(ls):
    '''Plots all diagrams in the list ls'''
    counter=0
    for g in ls:
        print("The "+str(counter)+"-th diagram is:")
        draw_diagrams(g)
        print(" ")
        counter+=1

def initialize_G(G,ilist):
    '''Initializes a nx.Graph() G with nodes from ilist as initial infected and all others susceptible.'''

    G.graph['coeff']=1
    G.graph['fun']=1
            
    for i in G.nodes:
        G.nodes[i]['state']="S"
    
  
    for i in ilist:
        G.nodes[i]['state']='X'
        

def add_new_node(g,i):
    '''Adds a new blue node to the directed diagram g with an edge from node i to the new node. If the node i is blue 'B', 
    it changes state to infected 'I'. Returns the new DiGraph gr'''
    #make copy of g
    gr=nx.DiGraph(coeff=g.graph['coeff'],fun=0)
    gr.add_edges_from(list(g.edges.data()))
    gr.add_nodes_from(list(g.nodes.data()))
    
    #add susceptible and edge to node i in g
    x=len(gr)
    gr.add_node(x,state='B')
    gr.add_edge(i,x)
    
    # change the state of the node i if it's blue
    if gr.nodes[i]['state'] == 'B' :
        gr.nodes[i]['state']='I'
    
    return gr

def add_int_edge(g,i,j):
    '''Adds an internal edge to the directed diagram g from node i to node j. If node i in blue 'B', the node changes state to infected 'I'. 
    The diagrams sign (coeff) is flipped. Returns the new DiGraph gr if node j is blue, returns False otherwise. '''
    #make copy of g
    gr=nx.DiGraph(coeff=g.graph['coeff'],fun=0)
    gr.add_edges_from(g.edges.data())
    gr.add_nodes_from(g.nodes.data())
    
    si = gr.nodes[i]['state']
    sj = gr.nodes[j]['state']
    if i==j: return False
    
    #add internal edge between i and j if it doesn't exist already
    if gr.has_edge(i,j): return False
    else:
        if (si,sj)==('I','B') or (si,sj)==('X','B'):
            gr.add_edge(i,j)
            gr.graph['coeff']=gr.graph['coeff']*(-1)
            #print('rule2')
            return gr
        
        elif (si,sj)==('B','B'):
            gr.add_edge(i,j)
            gr.nodes[i]['state']='I'
            gr.graph['coeff']=gr.graph['coeff']*(-1)
            #print('rule4')
            return gr
        else:
            return False



def grow_fund_GG(ls_lastord,m):
    '''Returns the fundamental graph of diagrams (GG) from the order inputted by the list ls_lastord, to the order (lastord + m). 
    The nodes of GG are directed diagrams and the parental relationships between diagrams are the directed edges of GG. The edge weights denotes the 
    multiplicities of the parent diagram in the recursive relation. The final GG will only contain the relations between the last order that 
    was memorized (lastord) and the last order wanted (lastord+m).'''
    pert_fund_sgr=[]
    pert_fund_sgr.append(ls_lastord)
    GG=nx.DiGraph()
    GG.add_nodes_from(ls_lastord)
    L=len(pert_fund_sgr[-1][-1].to_undirected().edges)
    
    #For loop for adding new diagrams
    start_time=time.time()
    for n in range(m):
        order_time=time.time()
        aux_fund_sgr=[]
        for graph in pert_fund_sgr[n]:
            #Create list of blue nodes
            b_ls=[i for i in graph.nodes if graph.nodes[i]['state'] == 'B']
            
            #For each node, add an edge to a new blue node:
            for i in graph.nodes:
                new_graph=diagram_to_directed(add_new_node(graph,i))     
                
                #Check if the new diagram is already in the list of diagrams. If it is, add an edge in GG to the existing diagram
                check=0
                for gr in reversed(aux_fund_sgr):
                    if are_the_same(gr,new_graph):
                        check+=1
                        GG.add_edge(gr,graph,weight=0)
                        break
                
                #If the diagram is new, add it to the list and create an edge in GG
                if check==0:
                    GG.add_edge(new_graph,graph,weight=0)
                    aux_fund_sgr.append(new_graph)
                    
            #For each blue node, add an internal edge to it if it is not already there
            for i in b_ls:
                for j in graph.nodes:
                    sec_graph=False
                    if i!=j and not graph.has_edge(j,i):
                        sec_graph = diagram_to_directed(add_int_edge(graph,j,i))
                    
                    if sec_graph :
                        #Check if the new diagram is already in the list of diagrams. If it is, add an edge in GG to the existing diagram
                        check=0
                        for gr in reversed(aux_fund_sgr):
                            if are_the_same(gr,sec_graph):
                                check+=1
                                GG.add_edge(gr,graph,weight=0)
                                break
                        #If the diagram is new, add it to the list and create an edge in GG
                        if check==0:
                            GG.add_edge(sec_graph,graph,weight=0)
                            aux_fund_sgr.append(sec_graph)
                            
        #Append new diagrams to the list of all diagrams
        pert_fund_sgr.append(aux_fund_sgr)
        order_fin=time.time()
        print('All {} graphs at order {} analyzed in {} sec'.format(len(pert_fund_sgr[n]),L+n,order_fin-order_time))
    print('{} nodes added in {} sec'.format(len(GG),order_fin-start_time))
    mid_time=time.time()
    #For loop for adding parent multiplicities and signs as edge weights
    for n in range(1,m+1):
        order_time=time.time()
        for son in pert_fund_sgr[n]:
            parent_ls=get_parents(son)
            for pp in parent_ls:
                for parent in GG[son]:
                    if are_the_same(pp,parent):
                        sign = parent.graph['coeff']*son.graph['coeff']
                        GG.edges[son,parent]['weight']+=sign
                        break
        order_fin=time.time()
        print("Parents multiplicities at order {} added in {} sec".format(L+n,order_fin-order_time))
    print("Finished in {} sec".format(order_fin-mid_time))
    return GG

def get_parents(graph):
    '''Returns the list of all parents of the directed diagram (graph). The parents consist out of all physical diagrams with 
    one less edge connecting to a blue node.'''
    par=[]
    b_ls = [i for i in graph.nodes if graph.nodes[i]['state'] == 'B']

    #For each blue node we remove each edge to a blue node and append the possible physical diagrams to the list par
    for l in b_ls:
        ls_neig=list(graph.predecessors(l))
        
        #If the blue node has only 1 neighbor, we remove the blue node altogether
        if len(ls_neig)==1:
            gr=nx.DiGraph(coeff=graph.graph['coeff'],fun=0)
            gr.add_nodes_from(graph.nodes.data())
            gr.add_edges_from(graph.edges.data())
            gr.remove_node(l)
            ls_neig_out=list(gr.successors(ls_neig[0]))
            ls_neig_in=list(gr.predecessors(ls_neig[0]))
            ls_doub_arr=[x for x in ls_neig_out if x in ls_neig_in]
            if nx.is_weakly_connected(gr):
                #If the neighbor of the blue node has no outgoing edges, it should change state to blue (unless it is a source)
                if len(ls_neig_out)==0:
                    if gr.nodes[ls_neig[0]]['state']!='X':
                        gr.nodes[ls_neig[0]]['state']='B'
                    if len(gr.nodes)==1:
                        par.append(gr)
                    
                    else:    
                        gr_aux=diagram_to_directed(gr)
                        par.append(gr_aux)
                else:
                    #If the outgoing edges of the neighbor are not double directed it should stay infected.
                    if not ls_doub_arr:
                        gr_aux=diagram_to_directed(gr)
                        par.append(gr_aux)
                    else:
                        #If all outgoing edges are double directed, the neighbor of the blue node can be either 'B' or 'I'
                        if len(ls_neig_out)==len(ls_doub_arr):
                            gr2=nx.DiGraph(coeff=gr.graph['coeff'],fun=0)
                            gr2.add_nodes_from(gr.nodes.data())
                            gr2.add_edges_from(gr.edges.data())
                            gr_aux=diagram_to_directed(gr)
                            par.append(gr_aux)
                            
                            
                            if gr2.nodes[ls_neig[0]]['state']!='X':
                                gr2.nodes[ls_neig[0]]['state']='B'
                                
                            gr2_aux=diagram_to_directed(gr2)
                            par.append(gr2_aux)
                        else:
                            gr_aux=diagram_to_directed(gr)
                            par.append(gr_aux)
        #If the blue node hase multiple neighbors, we instead remove only edge to the blue node. We repeat the above steps for changing the neighbors state
        else:
            for neig in ls_neig:
                gr=nx.DiGraph(coeff=graph.graph['coeff'],fun=0)
                gr.add_nodes_from(graph.nodes.data())
                gr.add_edges_from(graph.edges.data())
                gr.remove_edge(neig,l) 
                ls_neig_out=list(gr.successors(neig))
                ls_neig_in=list(gr.predecessors(neig))
                ls_doub_arr=[x for x in ls_neig_out if x in ls_neig_in]
                if nx.is_weakly_connected(gr):
                    #Change the sign of the diagram:
                    gr.graph['coeff']=gr.graph['coeff']*(-1)
                    
                    #If the neighbor of the blue node has no outgoing edges, it should change state to blue (unless it is a source)
                    if len(ls_neig_out)==0:
                        if gr.nodes[neig]['state']!='X':
                            gr.nodes[neig]['state']='B'
                        gr_aux=diagram_to_directed(gr)
                        par.append(gr_aux)
                    else:
                        #If the outgoing edges of the neighbor are not double directed it should stay infected.
                        if not ls_doub_arr:
                            gr_aux=diagram_to_directed(gr)
                            par.append(gr_aux)
                        else:
                            #If all outgoing edges are double directed, the neighbor of the blue node can be either 'B' or 'I'
                            if len(ls_neig_out)==len(ls_doub_arr):
                                gr2=nx.DiGraph(coeff=gr.graph['coeff'],fun=0)
                                gr2.add_nodes_from(gr.nodes.data())
                                gr2.add_edges_from(gr.edges.data())
                                
                                gr_aux=diagram_to_directed(gr)
                                par.append(gr_aux)
                                
                                if gr2.nodes[neig]['state']!='X':
                                    gr2.nodes[neig]['state']='B'
                                    
                                gr2_aux=diagram_to_directed(gr2)
                                par.append(gr2_aux)
                            else:
                                gr_aux=diagram_to_directed(gr)
                                par.append(gr_aux)
                    
    return par    

    
def paths_to_edges(paths):
    '''Creates a set of edges from a list of paths '''
    edges=set()
    for node_list in paths:
        for i in range(len(node_list)-1):
            edge=(node_list[i],node_list[i+1])
            edges.add(edge)
    
    return edges

def diagram_to_directed(diagram):
    '''Takes an undirected, but colored diagram and returns a directed diagram with edges pointing from source to sink. 
    All edges without a path from source to sink are purged from the diagram. Only works for diagrams with a single source'''
    diagram=diagram.to_undirected()
    pat_0ls=[i for i in diagram.nodes() if diagram.nodes[i]['state']=='X']
    blue_ls=[i for i in diagram.nodes() if diagram.nodes[i]['state']=='B']
    edge_list=set()
    for blue_node in blue_ls:
        ls_nodes=[x for x in diagram.nodes() if diagram.nodes[x]['state']=='I' or diagram.nodes[x]['state']=='S' ]
        ls_nodes.append(blue_node)
        ls_nodes.append(pat_0ls[0])
        sgr=diagram.subgraph(ls_nodes)
        paths=list(nx.all_simple_paths(sgr,pat_0ls[0],blue_node))
        
        edge_list=edge_list.union(paths_to_edges(paths))
    
    gr=nx.DiGraph()
    gr.add_edges_from(edge_list)
    gr.graph=diagram.graph
    for i in gr.nodes:
        gr.nodes[i]['state']=diagram.nodes[i]['state']
        
    return gr



def fix_coeff(diag):
    '''Computes the sign of the diagram diag and stores is a the graph attribute 'coeff'.'''
    gr=diag.to_undirected()
    m = len(gr.edges)
    n= len(gr.nodes)
    nX = len([1 for i in gr.nodes if gr.nodes[i]['state']=='X'])
    sign = (-1)**(m+n-nX)
    diag.graph['coeff']=sign


def parent_cont(GG,diag):
    '''Returns the sum of the contributions from the parent diagrams of a diag, given the graph of diagrams GG. 
    Takes into account the multiplicities and signs stored in the edge weights of GG'''
    #ls_GG=list(GG.nodes)
    parent_res=0
    if len(diag.edges)==0:
        return 1
    else:
        for parent in GG[diag]:
            coef_parent=GG.edges[diag,parent]['weight']
            parent_res+=coef_parent*parent.graph['fun']
        
        return sym.simplify(parent_res)


def fundamental_res(GG):
    '''Integrates the result as a function of t for all diagrams in GG. Stores the answer for each diagram as the graph attribute 'fun'.'''
    
    #ls_GG_nod=list(GG.nodes)
    
    for diag in GG.nodes:
        if diag.graph['fun'] == 0:
        
            c=sum([diag.degree(j) for j in diag.nodes if diag.nodes[j]['state']=='B'])
            d=parent_cont(GG,diag)

            diag.graph['fun']=sym.simplify(exp(-c*t)*sym.integrate(exp(c*s)*d.subs(t,s), (s,0,t)))



def diagram_to_tuple(diag):
    '''
    Converts a single diagram diag into a tuple of tuples which can be seen as a 
    matrix with the node states on the diagonal and the adjacency matrix of the 
    diagram on the off-diagonal elements
    '''
    g=diag.to_undirected()
    n= len(g)
    mat=[]
    for i in g.nodes:
        row=tuple()
        for j in g.nodes:
            if i==j: 
                row+=tuple(g.nodes[i]['state'])
            else:
                if g.has_edge(i,j):
                    row+=tuple('1')
                else: 
                    row+=tuple('0')
        mat.append(row)
    
    return tuple(mat)

def diagram_to_Matrix(g):
    '''
    Converts a single diagram g into a SymPy Matrix with the node states on the diagonal and the adjacency matrix of the 
    diagram as the off-diagonal elements. States are encoded as SymPy symbols x for 'X', y for 'I' and z for 'B'.
    '''
    
    n= len(g)
    mat=[]
    for i in g.nodes:
        row=[]
        for j in g.nodes:
            if i==j: 
                if g.nodes[i]['state']=='X':
                    row.append(x)
                elif g.nodes[i]['state']=='I':
                    row.append(y)
                else:
                    row.append(z)
            else:
                if g.has_edge(i,j):
                    row.append(1)
                else: 
                    row.append(0)
        mat.append(row)
    
    return Matrix(mat)

def diagram_to_det(g):
    '''Returns the determinant of the diagrams matrix. States are encoded as SymPy symbols x for 'X', y for 'I' and z for 'B'. '''
    a=diagram_to_Matrix(g.to_undirected())
    deta=a.det()
    return deta

def det_to_tuple(det):
    '''Returns the diagrams identifier as a tuple of three numbers, obtained by evaluating the determinant det 
    at three different points (x,y,z).'''
    numb1=int(det.evalf(subs={x:3,y:37,z:17}))
    numb2=int(det.evalf(subs={x:13,y:5,z:53}))
    numb3=int(det.evalf(subs={x:43,y:11,z:23}))
    tup= (numb1,numb2,numb3)
    return tup


def tuple_to_diagram(tup):
    '''
    Converts a tuple of tuples (as a matrix) into a nx.Graph of the corresponding undirected diagram
    '''
    y=np.array(tup)
    gr=nx.Graph()
    gr.add_nodes_from([i for i in range(len(y))])
    
    for i in range(len(y)):
        for j in range(i,len(y)):
            if i==j :
                gr.nodes[i]['state'] = y[i,i]            
            elif y[i,j]=="1" : gr.add_edge(i,j,weight=1)
    
    return gr
    
def write_diagrams(filename,GG):
    '''Writes the diagrams in the graph of graphs into the file `filename` '''
    
    w = csv.writer(open(filename, "w"))
    
    for diagram in GG.nodes:
        key=diagram_to_tuple(diagram)
        det_diag=diagram_to_det(diagram)
        identifier=det_to_tuple(det_diag)
        order=len(diagram.to_undirected().edges)
        w.writerow([order,key,diagram.graph['coeff'],diagram.graph['fun'],det_diag,identifier])


def write_diagrams_from_ls(filename,ls_GG):
    '''Writes the diagrams from a list of diagrams organized by order	 into the file `filename` '''
    
    w = csv.writer(open(filename, "w"))
    
    for i in range(len(ls_GG)):
        for sgr in ls_GG[i]:
            diagram=sgr.to_undirected()
            key=diagram_to_tuple(diagram)
            det_diag=diagram_to_det(diagram)
            identifier=det_to_tuple(det_diag)
            order=len(diagram.edges)
            w.writerow([order,key,diagram.graph['coeff'],diagram.graph['fun'],det_diag,identifier])


def lookup_diagram(filename,sgr):
    '''Reads the file `filename` and returns a list of nx.Graphs of the diagrams
    with the sign and function of t as graph attributes'''
    t = symbols('t', real=True)
    x, y, z = symbols('x y z', real=True)
    
    diagram=sgr.to_undirected()
    reader = csv.reader(open(filename))
    order_diag=len(diagram.edges())
    diag_det=diagram_to_det(diagram)
    diag_id=det_to_tuple(diag_det)
    count=0
    
    for row in reader:
        order_entry=eval(row[0])
        
        if order_diag == order_entry:
            entry_id=eval(row[-1])
            if entry_id==diag_id:
                entry_diag=tuple_to_diagram(eval(row[1]))
                
                if are_the_same(entry_diag,diagram):
                    sgr.graph['coeff'] =  eval(row[2])
                    sgr.graph['fun'] = sympify(row[3]).subs('t',t)
                    count=1
                    break
                else :
                    continue

    if count==0: 
        print("OOPS! Diagram not found!")
        draw_diagrams(diagram)

        
def read_diagrams(filename,upto_order):
    '''Reads the file `filename` and returns a list of nx.Graphs of the diagrams
    with the sign and function of t as graph attributes'''
    n, m = symbols('n m', integer=True)
    t = symbols('t', real=True)
    
    reader = csv.reader(open(filename))
    
    ls_GG=[]
    for i in range(upto_order+1):
        ls_GG.append([])
    
    for row in reader:
        order_diag=eval(row[0])
               
        if order_diag <= upto_order:
            diag_tup = eval(row[1])
            diag_dict = {'coeff' :  eval(row[2])}
            diag_dict['fun'] = sympify(row[3]).subs('t',t)
            gr=diagram_to_directed(tuple_to_diagram(diag_tup))
            gr.graph=diag_dict
            ls_GG[order_diag].append(gr)
        else:
            break    
    
    return ls_GG
    
def read_diagrams_at_order(filename,order):
    '''Reads the file `filename` and returns a list of nx.Graphs of the diagrams
    with the sign and function of t as graph attributes'''
    n, m = symbols('n m', integer=True)
    t = symbols('t', real=True)
    
    reader = csv.reader(open(filename))
    
    ls_GG=[]
    
    for row in reader:
        order_diag=eval(row[0])
               
        
        if order_diag == order:
            diag_tup = eval(row[1])
            diag_dict = {'coeff' :  eval(row[2])}
            diag_dict['fun'] = sympify(row[3]).subs('t',t)
            gr=diagram_to_directed(tuple_to_diagram(diag_tup))
            gr.graph=diag_dict
            ls_GG.append(gr)
        elif order_diag > order:
            break 
    
    return ls_GG


def append_diagram_to_file(filename,sgr):
    '''Appends `diagram` as a new row to the file `filename` '''
    diagram=sgr.to_undirected()
    key=diagram_to_tuple(diagram)
    order=len(diagram.edges)
    det_diag=diagram_to_det(diagram)
    identifier=det_to_tuple(det_diag)   
    newline = [order,key,diagram.graph['coeff'],diagram.graph['fun'],det_diag,identifier]
    
    # Open file in append mode
    with open(filename, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(newline)

    
def op_contributions(G,pat_0,blue_node):
    '''Returns a list of ccontributions to the expectation value for the node blue_node given that pat_0 is the initial infected node in a network G '''
    #Initialize the states of the nodes in the network
    G.nodes[pat_0]['state']='X'
    G.nodes[blue_node]['state']='B'
    for i in G.nodes:
        if i != pat_0 and i !=blue_node:
            G.nodes[i]['state']='I'
    
    #Add the greatest possible contribution to the list
    diag=diagram_to_directed(G)
    fix_coeff(diag)
    contributions=[diag]
    diag.graph['fun']=0
    number=0
    
    #For each diagram in the list of contributions, remove one edge and add it to the list if it is a valid and unique contribution
    for diagram in contributions:
        edge_list=list(diagram.edges())
        number+=1
        for edge in edge_list:
            #Create a copy of the diagram and remove one edge from the edge_list
            gr=nx.Graph(coeff='?', fun=0)
            gr.add_edges_from(edge_list)
            gr.add_nodes_from(diagram.nodes.data())
            gr.remove_edge(edge[0],edge[1])
            if nx.is_connected(gr):

                gr=diagram_to_directed(gr)

                #Check whether the labelled diagram is already in the list of contributions
                count=0
                for sgr in contributions:
                    if sgr.edges==gr.edges:
                        count=1
                        break

                #If the diagram is not in the list, add it to the list of contributions
                if count==0 and nx.is_weakly_connected(gr):
                    fix_coeff(gr)
                    contributions.append(gr)
        #print("Analyzing diagram {} out of list of {} contributions".format(number,len(contributions)))
    
    return contributions

def sort_contributions(contributions, max_order):
    '''Sorts a list of diagrams (contributions) into a list of lists of diagrams organized by the order of the diagram. 
    max_order should be the highest order diagram in the list.'''
    ls_cont=[]
    for i in range(max_order+1):
        ls_cont.append([])
    
    for diagram in contributions:
        order_diag= len(diagram.to_undirected().edges)
        ls_cont[order_diag].append(diagram)
    
    
    return ls_cont

def check_ineq_sgr(ls_GG_nod):
    '''Returns a list of the inequivalent diagrams in the list of graphs ls_GG_nod.'''
    ls_inequiv=[]
    for i in range(len(ls_GG_nod)):
        check=0
        for j in range(i):
            if are_the_same(ls_GG_nod[i],ls_GG_nod[j]):
                ls_GG_nod[i].graph=ls_GG_nod[j].graph
                check=1
                break
        if check==0:
            ls_inequiv.append(ls_GG_nod[i])
    
    return ls_inequiv

    
def compute_op(G,pat_0,blue):
    '''Computes the infected expectation value of the node 'blue' as a function of t for logistic (SI) growth on 
    a graph G with pat_0 initially infected.
    Returns the expectation value as a SymPy object and the list of all contributing diagrams and their parents.'''
    origin_time=time.time()
    order_list=7
    
    #Collect all contributing diagrams, and sort them by order
    print("Computing contributing diagrams")
    contributing_digraphs=op_contributions(G,pat_0,blue)
    ls_cont=sort_contributions(contributing_digraphs,len(contributing_digraphs[0].to_undirected().edges))
    cont_time=time.time()
    print("There are {} contributing diagrams, computed in {}s".format(len(contributing_digraphs),cont_time-origin_time))
    
    #Add only the diagrams of inequivalent topology to the graph of diagrams
    print("Making list of inequivalent diagrams")
    pert_fund_sgr=[]
    
    for ls in ls_cont:
        ls_ineq=check_ineq_sgr(ls)
        pert_fund_sgr.append(ls_ineq)
        
    no_ineq_diag=sum([len(i) for i in pert_fund_sgr])
    print("There are {} inequivalent diagrams".format(no_ineq_diag))
    
    print("Computing graph of graphs for diagrams at order > {}".format(order_list))
    begin=time.time()
    GG=nx.DiGraph()

    #For loop for contructing the graph of diagrams: find all the parents of the contributing diagrams (and their parents)
    max_order=len(pert_fund_sgr)-1
    for i in reversed(range(order_list+1,max_order+1)):
        start=time.time()
        for diagram in pert_fund_sgr[i]:
            par=get_parents(diagram)
            unique_par=check_ineq_sgr(par)

            for p in unique_par:
                sign = diagram.graph['coeff']*p.graph['coeff']

                check=0
                for sgr in pert_fund_sgr[i-1]:
                    if are_the_same(sgr,p):
                        check=1
                        GG.add_edge(diagram,sgr,weight=0)
                        for p2 in par:
                            if are_the_same(p,p2):
                                GG.edges[diagram,sgr]['weight']+=sign
                        break
                    else: continue

                if check==0:
                    pert_fund_sgr[i-1].append(p)

                    GG.add_edge(diagram,p,weight=0)
                    for p2 in par:
                        if are_the_same(p,p2):
                            GG.edges[diagram,p]['weight']+=sign

        mid_time=time.time()                        
        print("Finished {} diagrams at order {} in t = {} s".format(len(pert_fund_sgr[i]),i,mid_time-start))
    end_time=time.time()
    
    print("GG computed in t={}".format(end_time-begin))
    
    len_pfs=sum([len(i) for i in pert_fund_sgr])
    no_integrals=sum([len(pert_fund_sgr[i]) for i in range(order_list+1,len(pert_fund_sgr))])
    print("There are {} diagrams of which {} unknown".format(len_pfs,no_integrals))
    
    print("Looking up known contributions")
    #Look up the known contributions from the file 'diagrams.csv'
    st=time.time()
    for i in range(min(order_list+1,max_order+1)):
        print('Looking up order ',i)
        for diagram in pert_fund_sgr[i]:
            lookup_diagram('diagrams.csv',diagram)
    et=time.time()
    print('Assigned known (order {}) results to diagrams in {} seconds'.format(order_list,et-st))
    
    #Integrate unknown diagrams
    print("Integrating unknown diagrams in the graph of graphs.")
    st=time.time()
    for i in range(order_list+1,max_order+1):
        st_order=time.time()
        for diagram in pert_fund_sgr[i]:
            c=sum([diagram.degree(j) for j in diagram.nodes if diagram.nodes[j]['state']=='B'])
            d=parent_cont(GG,diagram)
            diagram.graph['fun']=sym.simplify(exp(-c*t)*sym.integrate(exp(c*s)*d.subs(t,s), (s,0,t)))
        et_order=time.time()
        print("{} diagrams at order {} integrated in t={}s".format(len(pert_fund_sgr[i]),i,et_order-st_order))
    et=time.time()
    print('Computed the results for unknown diagrams (order > {}) in {} seconds'.format(order_list,et-st))
    print('Computing one-point function by adding contributing diagrams')
    
    #Add all functions from the list of contributing diagrams 
    fin_res=0
    batch=50
    count=0
    for elem in contributing_digraphs:
        count+=1
        fin_res+= elem.graph['fun']
        if count==batch:
            fin_res=sym.simplify(fin_res)
            batch+=50
    fin_res=sym.simplify(fin_res)
    et=time.time()
    print('Done! The final result:',fin_res)
    print('It took {} seconds.'.format(et-origin_time))
    
    return fin_res, pert_fund_sgr



def tot_op_contributions(G,pat_0):
    '''Returns a list of ccontributions to the expectation value for the total infected probability 
    given that pat_0 is the initial infected node in a network G'''
    initialize_G(G,[pat_0])
    contributions=[]
    
    #Initialize the list of contributions s.t. pat_0 is initially infected and any other node is a sink (blue node)
    for i in G.nodes:
        if i != pat_0:
            gr=nx.Graph(coeff=1,fun=0)
            gr.add_edges_from(G.edges.data())
            gr.add_nodes_from(G.nodes.data())
            gr.nodes[i]['state']='B'
    
            diag=diagram_to_directed(gr)
            for j in diag.nodes:
                if j != pat_0 and j !=i:
                    diag.nodes[j]['state']='I'
    
            fix_coeff(diag)
            contributions.append(diag)
    
    number=0
    
    #Systematically remove edges to find all contributing diagrams to the total expectation value
    for diagram in contributions:
        edge_list=list(diagram.edges())
        number+=1
        for edge in edge_list:
            gr=nx.Graph(coeff=1, fun=0)
            gr.add_edges_from(edge_list,weight=1)
            gr.add_nodes_from(diagram.nodes.data())
            gr.remove_edge(edge[0],edge[1])
            if nx.is_connected(gr):

                diag=diagram_to_directed(gr)
                if len(diag.edges)==0:
                    break
                    
                count=0
                for sgr in contributions:
                    if sgr.edges==diag.edges:
                        count=1
                        break

                if count==0 and nx.is_weakly_connected(diag):
                    fix_coeff(diag)
                    contributions.append(diag)
    
    return contributions

def compute_tot_op(G,pat_0):
    '''Computes the total infected expectation value as a function of t for logistic (SI) growth on 
    a graph G with pat_0 initially infected.
    
    Returns the sum of all individual expectation values as a SymPy object, as well as the list of 
    contributing diagrams and their parents'''
    origin_time=time.time()
    order_list=10
    max_order=len(G.edges)
    N=int(len(G.nodes))
    
    print("Computing contributing diagrams")
    contributing_digraphs=tot_op_contributions(G,pat_0)
    ls_cont=sort_contributions(contributing_digraphs,max_order)
    cont_time=time.time()
    print("There are {} contributing diagrams, computed in {}s".format(len(contributing_digraphs),cont_time-origin_time))
    
    print("Making list of inequivalent diagrams")
    pert_fund_sgr=[]
    
    for ls in ls_cont:
        ls_ineq=check_ineq_sgr(ls)
        pert_fund_sgr.append(ls_ineq)
        
    no_ineq_diag=sum([len(i) for i in pert_fund_sgr])
    print("There are {} inequivalent diagrams".format(no_ineq_diag))
    
    print("Computing graph of graphs for diagrams at order > {}".format(order_list))
    begin=time.time()
    GG=nx.DiGraph()

#     max_order=len(pert_fund_sgr)-1
    for i in reversed(range(order_list+1,max_order+1)):
        start=time.time()
        for diagram in pert_fund_sgr[i]:
            par=get_parents(diagram)
            unique_par=check_ineq_sgr(par)

            for p in unique_par:
                #print("Parent found")
                sign = diagram.graph['coeff']*p.graph['coeff']

                check=0
                for sgr in pert_fund_sgr[i-1]:
                    if are_the_same(sgr,p):
                        check=1
                        #print("Already got this parent")
                        GG.add_edge(diagram,sgr,weight=0)
                        for p2 in par:
                            if are_the_same(p,p2):
                                GG.edges[diagram,sgr]['weight']+=sign
                        break
                    else: continue

                if check==0:
                    #print("Found new parent")
                    pert_fund_sgr[i-1].append(p)

                    GG.add_edge(diagram,p,weight=0)
                    for p2 in par:
                        if are_the_same(p,p2):
                            GG.edges[diagram,p]['weight']+=sign

        mid_time=time.time()                        
        print("Finished {} diagrams at order {} in t = {} s".format(len(pert_fund_sgr[i]),i,mid_time-start))
    end_time=time.time()
    
    print("GG computed in t={}".format(end_time-begin))
    
    len_pfs=sum([len(i) for i in pert_fund_sgr])
    no_integrals=sum([len(pert_fund_sgr[i]) for i in range(order_list+1,len(pert_fund_sgr))])
    print("There are {} diagrams of which {} unknown".format(len_pfs,no_integrals))
    print("Looking up known contributions")
    st=time.time()
    for i in range(min(order_list+1,max_order+1)):
        print('Looking up order ',i)
        for diagram in pert_fund_sgr[i]:
            lookup_diagram('diagrams.csv',diagram)
    et=time.time()
    print('Assigned known (order {}) results to diagrams in {} seconds'.format(order_list,et-st))
       
    print("Integrating unknown diagrams in the graph of graphs.")
    st=time.time()
    for i in range(order_list+1,max_order+1):
        st_order=time.time()
        for diagram in pert_fund_sgr[i]:
            #if diagram.graph['fun'] == '?':
            c=sum([diagram.degree(j) for j in diagram.nodes if diagram.nodes[j]['state']=='B'])
            d=parent_cont(GG,diagram)
            diagram.graph['fun']=sym.simplify(exp(-c*t)*sym.integrate(exp(c*s)*d.subs(t,s), (s,0,t)))
        et_order=time.time()
        print("{} diagrams at order {} integrated in t={} s".format(len(pert_fund_sgr[i]),i,et_order-st_order))
    et=time.time()
    print('Computed the results for unknown diagrams (order > {}) in {} seconds'.format(order_list,et-st))
    print('Computing one-point function by adding contributing diagrams')
        
    fin_res=1
    count=0
    batch=50
    
    for elem in contributing_digraphs:
        count+=1
        fin_res+= elem.graph['fun']
        if count==batch:
            fin_res=sym.simplify(fin_res)
            batch+=50
       
    
    fin_res=sym.simplify(fin_res)
    et=time.time()
    print('Done! The final result:',fin_res)
    print('It took {} seconds'.format(et-origin_time))
    
    return fin_res, pert_fund_sgr


def lookup_ff_op(G,pat_0,blue):
    '''Looks up the diagrams contributing to the infected expectation value of the node 'blue' 
    as a function of t for logistic (SI) growth on the graph of Florentine families with pat_0 initially infected'''
    origin_time=time.time()
    
    
    print("Computing contributing diagrams")
    contributing_digraphs=op_contributions(G,pat_0,blue)
    ls_cont=sort_contributions(contributing_digraphs,len(contributing_digraphs[0].to_undirected().edges))
    cont_time=time.time()
    print("There are {} contributing diagrams, computed in {}s".format(len(contributing_digraphs),cont_time-origin_time))
    
    print("Making list of inequivalent diagrams")
    pert_fund_sgr=[]
    
    for ls in ls_cont:
        ls_ineq=check_ineq_sgr(ls)
        pert_fund_sgr.append(ls_ineq)
        
    no_ineq_diag=sum([len(i) for i in pert_fund_sgr])
    print("There are {} inequivalent diagrams".format(no_ineq_diag))
    max_order=len(pert_fund_sgr)-1
    print("Looking up contributions")
    st=time.time()
    for i in range(max_order+1):
        #print('Looking up order ',i)
        for diagram in pert_fund_sgr[i]:
            lookup_diagram('diagrams_ff.csv',diagram)
    et=time.time()
    print('Assigned known results to diagrams in {} seconds'.format(et-st))
       
    fin_res=0
    batch=50
    count=0
    for elem in contributing_digraphs:
        count+=1
        fin_res+= elem.graph['fun']
        if elem.graph['fun'] == 0 or elem.graph['fun']=='?':
            print("Graph not found!")
            draw_diagrams(elem)
            print("Diagram tuple:",diagram_to_tuple(elem))
            print("Diagram fun:",elem.graph['fun'])
            
        if count==batch:
            fin_res=sym.simplify(fin_res)
            batch+=50
    
    fin_res=sym.simplify(fin_res)
    et=time.time()
    print('Done! The final result:',fin_res)
    print('It took {} seconds.'.format(et-origin_time))
    
    return fin_res


def flip_BX(ex):
    '''Flips all blue nodes with all initial infected nodes in the (undirected!) diagram ex.'''
    ls_blue=[i for i in ex.nodes if ex.nodes[i]['state']=='B']
    if len(ls_blue)==0: return ex
    ls_init=[i for i in ex.nodes if ex.nodes[i]['state']=='X']
    sign=int((-1)**(len(ls_blue)-len(ls_init)))
    
    gr=nx.Graph(coeff=sign*ex.graph['coeff'],fun=sign*ex.graph['fun'])
    gr.add_nodes_from(ex.nodes.data())
    gr.add_edges_from(ex.edges.data())
    for i in gr.nodes:
        if gr.nodes[i]['state']=='B':
            gr.nodes[i]['state']='X'
        elif gr.nodes[i]['state']=='X':
            gr.nodes[i]['state']='B'
        else: continue
    return gr

def branch_at_pat0(diagram):
    '''Checks whether a diagram branches at an initially infected node. Returns True if diagram branches, False otherwise'''
    patient0=[i for i in diagram.nodes if diagram.nodes[i]['state']=='X']
    
    gr= nx.Graph()
    gr.add_nodes_from(diagram.nodes.data())
    gr.add_edges_from(diagram.edges.data())
    gr.graph=diagram.graph
    gr.remove_nodes_from(patient0)
    if len(gr.nodes)==0 :
        return False

    elif not nx.is_connected(gr) : 
        return True
    else :
        return False
        
def give_init_branch(diagram):
    '''Gives a list of diagrams of branches from patient 0 of the input `diagram`. Does not compute
    the contributions from the branch diagrams'''
    
    if not branch_at_pat0(diagram): return [diagram]
    
    patient0=[i for i in diagram.nodes if diagram.nodes[i]['state']=='X']
    ls_comp=[]
    gr= nx.Graph()
    gr.add_nodes_from(diagram.nodes.data())
    gr.add_edges_from(diagram.edges.data())
    gr.graph=diagram.graph
    gr.remove_nodes_from(patient0)
    
    if not nx.is_connected(gr):
        #pat0_nbrs=set(diagram[i])
        components=nx.connected_components(gr)

        for comp_nodes in components:
            sgr=gr.subgraph(comp_nodes)
            comp_cont=nx.Graph(coeff='?',fun=0)
            comp_cont.add_edges_from(sgr.edges.data())
            comp_cont.add_nodes_from(sgr.nodes.data())
            #comp_cont.add_node(i,state='X')
            
            for i in patient0:
                for j in set(diagram[i]) & comp_nodes:
                    comp_cont.add_edge(i,j,weight=1)
                    comp_cont.nodes[i]['state']='X'
            
            fix_coeff(comp_cont)
            ls_comp.append(comp_cont)

    return ls_comp




def has_blue_loop(diagram):
    '''Checks whether diagram has a blue loop. Returns True if it does and False otherwise.'''
    ls_bl_deg=[diagram.degree(i) for i in diagram.nodes if diagram.nodes[i]['state']=='B']
    for i in ls_bl_deg:
        if i>1: return True
    
    return False 

def cut_open_loops(diagram):
    '''Opens up blue loops in the diagram. Returns the diagram with all blue loops opened up.'''
    if not has_blue_loop(diagram) :return diagram
    
    ls_bl=[i for i in diagram.nodes if diagram.nodes[i]['state']=='B']
    gr= nx.Graph(coeff=diagram.graph['coeff'],fun=diagram.graph['fun'])
    gr.add_nodes_from(diagram.nodes.data())
    gr.add_edges_from(diagram.edges.data())
    
    for i in ls_bl:
        if diagram.degree(i)>1:
            
            ls_nbr=[j for j in diagram[i]]
            
            
            for k in range(1,len(ls_nbr)):
                gr.remove_edge(i,ls_nbr[k])
                gr.add_node('nn'+str(i)+str(k),state='B')
                gr.add_edge(ls_nbr[k],'nn'+str(i)+str(k),weight=1)
                gr.graph['coeff'] = - gr.graph['coeff']
                gr.graph['fun'] = - gr.graph['fun']
    
    return gr          

def has_init_loop(diagram):
    '''Checks whether diagram has a loop containing a source node. Returns True if it does, False otherwise '''
    patient0=[i for i in diagram.nodes if diagram.nodes[i]['state']=='X']
    
    if all([diagram.degree(i)<=1 for i in patient0]): return False
    elif branch_at_pat0(diagram): return False
    else:      
        return True 

def cut_init_loops(diagram):
    '''Opens up loops at the source node. Returns a diagrams with the source loop opened up'''
    if not has_init_loop(diagram) :return diagram
    
    ls_init=[i for i in diagram.nodes if diagram.nodes[i]['state']=='X']
    
    gr= nx.Graph(coeff=diagram.graph['coeff'],fun=diagram.graph['fun'])
    gr.add_nodes_from(diagram.nodes.data())
    gr.add_edges_from(diagram.edges.data())
    
       
    for i in ls_init:
        if diagram.degree(i)>1:
            
            ls_nbr=[j for j in diagram[i]]
            
            
            for k in range(1,len(ls_nbr)):
                gr.remove_edge(i,ls_nbr[k])
                gr.add_node('nn'+str(i)+str(k),state='X')
                gr.add_edge(ls_nbr[k],'nn'+str(i)+str(k),weight=1)

    return gr          


def is_cuttable(diagram):
    '''Checks whether diagram has a source node connected with a single edge to a subgraph with 2 or more blue nodes. Returns boolean.'''
    #gr=cut_open_loops(diagram)
    ls_blue_deg=[diagram.degree(i) for i in diagram.nodes if diagram.nodes[i]['state']=='B']
    pat0=[i for i in diagram.nodes if diagram.nodes[i]['state']=='X']
    if len(pat0)>1: return False
    elif diagram.degree(pat0[0])> 1 or sum(ls_blue_deg) < 2 :
        return False
    else: 
        return True

def decomposes(diagram):
    '''Checks whether diagram can be decomposed using any of the symmetry relations. Returns boolean.'''
    if has_blue_loop(diagram) :
        return True
    elif nx.is_tree(diagram) : 
        return True
    
    elif branch_at_pat0(diagram):
        return True
    
    elif branch_at_pat0(flip_BX(diagram)):
        return True
    
    elif is_cuttable(diagram):
        return True
    
    elif is_cuttable(flip_BX(diagram)):
        return True
    else:
        return False




print("SI_script loaded!")


