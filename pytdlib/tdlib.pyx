r"""
Interface with TdLib (algorithms for tree decompositions)

This module defines functions based on TdLib, a
library that implements algorithms for tree
decompositions written by Lukas Larisch. 

Definition:

A tree decomposition of a graph G is a pair (T, b) consisting of a tree T and 
a function b: V(T) -> 2^{V(G)} associating with each node t in V(T) a set of 
vertices b(t), that are subsets of V(G) such that 

(T1) for every edge e in E(G) there is a node t in V(T) with e is a subset of 
b(t), and

(T2) for all v in V(G) the set b^{-1} := {t in V(T): v in b(t)} is non-empty 
and connected in T.

The width of (T, b) is defined as max{|b(t)|-1: t in V(T) }.
The treewidth of G is defined as the minimum width over all tree decompositions 
of G.

Some known results:

    - Trees have treewidth 1
    - Cycles have treewidth 2
    - Series-parallel graphs have treewidth at most 2
    - Cliques must be contained in some bag of a tree decomposition

Computing the treewidth or a tree decomposition of a given graph is NP-hard in
general.

This module containes the following functions** :

    - preprocessing   
        Applies save reduction rules to a given graph

    - PP_MD
        Applies preprocessing followed by the minDegree-heuristic to a given 
        graph

    - PP_FI_TM
        Applies preprocessing followed by the minFill-heuristic followed by 
        triangulation minimization to a given graph

    - preprocessing_glue_bags
        Glues the given bags with a given tree decomposition according to 
        subset-relation

    - lower_bound
        Computes a lower bound with respect to treewidth of a given graph

    - treedecomposition_exact_cutset
        Computes a tree decomposition of exact width of a given graph 
        (faster than the dynamic version in practice)

    - treedecomposition_exact_dynamic
        Computes a tree decomposition of exact width of a given graph 
        (asymptotically faster than the greedy version)

    - seperator_algorithm
        Computes a 4-approximate tree decomposition of a given graph

    - minDegree_ordering
        Computes an elimination ordering according to the minDegree-heuristic
    
    - fillIn_ordering
        Computes an elimination ordering according to the minFill-heuristic

    - MSVS
        Possibly reduces the width of a given tree decomposition with help 
        of minimal seperating vertex sets
    
    - minimalChordal
        Possibly reduces the width of a given tree decomposition by 
        triangulation minimization

    - treedec_to_ordering
        Computes an elimination ordering out of a tree decomposition
    
    - ordering_to_treedec
        Computes a tree decomposition out of an elimination ordering

    - trivial_decomposition
        Returns a trivial decomposition of a given graph
    
    - get_width
        Returns the width of a given tree decomposition
    
    - is_valid_treedecomposition
        Checks, if a tree decomposition is valid with respect to a given graph

AUTHOR: Lukas Larisch (now): Initial version
-------
"""

from libcpp.vector cimport vector

from tdlib cimport gc_preprocessing
from tdlib cimport gc_PP_MD
from tdlib cimport gc_PP_FI_TM
from tdlib cimport gc_preprocessing_glue_bags

from tdlib cimport gc_deltaC_min_d
from tdlib cimport gc_deltaC_max_d
from tdlib cimport gc_deltaC_least_c

from tdlib cimport gc_LBN_deltaC
from tdlib cimport gc_LBNC_deltaC
from tdlib cimport gc_LBP_deltaC
from tdlib cimport gc_LBPC_deltaC

from tdlib cimport gc_exact_decomposition_cutset
from tdlib cimport gc_exact_decomposition_dynamic

from tdlib cimport gc_seperator_algorithm

from tdlib cimport gc_minDegree_ordering
from tdlib cimport gc_fillIn_ordering

from tdlib cimport gc_is_valid_treedecomposition
from tdlib cimport gc_trivial_decomposition

##############################################################
############ GRAPH/DECOMPOSITION ENCODING/DECODING ###########
#the following will be used implicitly do the translation
#between the python graph encoding and TdLib graph encoding,
#which is based on BGL

cdef cython_make_tdlib_graph(pyV, pyE, vector[unsigned int] &V, vector[unsigned int] &E):
    for v in pyV:
        V.push_back(v)

    if len(pyE) > 0:
        #tuple representation
        if str(type(pyE[0])) == "<type 'tuple'>":
            for u,v in pyE:
                E.push_back(u)
                E.push_back(v)

        #list representation
        elif str(type(pyE[0])) == "<type 'list'>":
            for e in pyE:
                E.push_back(e[0])
                E.push_back(e[1])

        #internal representation (unfolded tuple/list representation)
        elif str(type(pyE[0])) == "<type 'int'>":
            for e in pyE:
                E.push_back(e)

cdef cython_make_tdlib_decomp(pyV, pyE, vector[vector[int]] &V, vector[unsigned int] &E):
    for t in pyV:
        V.push_back(t)

    #tuple representation
    if len(pyE) > 0 and str(type(pyE[0])) == "<type 'tuple'>":
        for u,v in pyE:
            E.push_back(u)
            E.push_back(v)

    #list representation
    elif len(pyE) > 0 and str(type(pyE[0])) == "<type 'list'>":
        for e in pyE:
            E.push_back(e[0])
            E.push_back(e[1])

    #internal representation (unfolded tuple/list representation)
    elif len(pyE) > 0 and str(type(pyE[0])) == "<type 'int'>":
        for e in pyE:
            E.push_back(e)

##############################################################
############ PREPROCESSING ###################################

def preprocessing(V, E):
    """
    Returns a possibly smaller instance of a given graph G and an encoding 
    of the parts of a tree decomposition, that could be computed so far.  
    If the treewidth of G is at most 3, preprocessing will reduce G
    to the empty graph and deliver the full list of bags, from which a 
    tree decomposition of exact width can be made of. 
    Otherwise, the algorithm may return a smaller instance of the original
    graph and the results of the reductions, that have been made so far
    as a list of bags.

    INPUT:

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - A tuple (V(G'), E(G'), encoded_bags, lb) where (V(G'), E(G') is the 
      reduced instance of G, encoded_bags is an encoding of the bags, 
      that belong to a tree decomposition of the whole graph G and lb
      is a lower bound on tw(G) iff V(G') is not empty, otherwise tw(G). 

    EXAMPLES:

        V_T, E_T, bags, lb = tdlib.preprocessing(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G
    cdef vector[vector[int]] c_bags;

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = -1

    py_lb = gc_preprocessing(V_G, E_G, c_bags, c_lb)

    return V_G, E_G, c_bags, py_lb

def PP_MD(V, E):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced 
    instance G' of G will be processed by the minDegree heuristic, which
    successivly eliminates a vertex of minimal degree. The returned tree
    decomposition then may be of non-optimal width. 

    INPUT:

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - A tuple (V, E, width)

    EXAMPLES:

        V_T, E_T, width = tdlib.PP_MD(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = -1

    py_lb = gc_PP_MD(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, py_lb

def PP_FI_TM(V, E):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced 
    instance G' of G will be processed by the fillIn heuristic, which
    successivly eliminates a vertex, that will cause least new edges
    within the elimination process. The resulting tree decomposition 
    will be postprocessed by the minimalChordal algorithm, that may
    reduce the width of the tree decomposition.

    INPUT:

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - A tuple (V, E, width)

    EXAMPLES:

        V_T, E_T, width = tdlib.PP_FI_TM(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = -1

    py_lb = gc_PP_FI_TM(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, py_lb


def preprocessing_glue_bags(V, E, bags):
    """
    Glues a list of bags returned by preprocessing with T according to the 
    subset-relation. This can be used to make a tree decomposition out of the 
    results of tdlib.preprocessing() for a given graph G.

    INPUTS:

    - V : a list of bags

    - E : a list of edges

    OUTPUT:

    - A tree decomposition (V, E), containing the given bags

    EXAMPLES::

        V_T, E_T, bags, lb = tdlib.preprocessing(V_G, E_G)
        ...
        V_T, E_T = tdlib.preprocessing_glue_bags(V_T2, E_T2, bags)
    """

    cdef vector[unsigned int] E_T
    cdef vector[vector[int]] V_T, c_bags

    for i in range(0, len(bags)):
        c_bags.push_back(bags[i])
    
    gc_preprocessing_glue_bags(V_T, E_T, c_bags)

    return V_T, E_T


##############################################################
############ LOWER BOUNDS ####################################

def lower_bound(V, E, algorithm = "deltaC_least_c"):
    """
    Calls one of the following algorithms to compute a lower bound on the 
    treewidth of a given graph:

        - deltaC_min_d
        - deltaC_max_d
        - deltaC_least_c
        - LBN_deltaC
        - LBNC_deltaC
        - LBP_deltaC
        - LBPC_deltaC

    INPUTS:

    - V : a list of vertices

    - E : a list of edges

    - algorithm -- (default: 'deltaC_least_c') specifies the algorithm to use 
                   for computing a lower bound on the treewidth of G. The 
                   algorithms have to be choosen among the above mentioned. 

    OUTPUT:

    - A lower bound on the treewidth of G

    EXAMPLES:
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        lb = tdlib.lower_bound(g, "deltaC_max_d")
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
    """

    cdef vector[unsigned int] V_G, E_G
    cython_make_tdlib_graph(V, E, V_G, E_G)
    cdef int c_lb = 0

    if(algorithm == "deltaC_min_d"):
        c_lb = gc_deltaC_min_d(V_G, E_G)
    elif(algorithm == "deltaC_max_d"):
        c_lb = gc_deltaC_max_d(V_G, E_G)
    elif(algorithm == "deltaC_least_c"):
        c_lb = gc_deltaC_least_c(V_G, E_G)
    elif(algorithm == "LBN_deltaC"):
        c_lb = gc_LBN_deltaC(V_G, E_G)
    elif(algorithm == "LBNC_deltaC"):
        c_lb = gc_LBNC_deltaC(V_G, E_G)
    elif(algorithm == "LBP_deltaC"):
        c_lb = gc_LBP_deltaC(V_G, E_G)
    elif(algorithm == "LBPC_deltaC"):
        c_lb = gc_LBPC_deltaC(V_G, E_G)
    else:
        print("Invalid lower bound algorithm")
        return -2

    py_lb = c_lb

    return py_lb


##############################################################
############ EXACT ALGORITHMS ################################

def exact_decomposition_cutset(V, E, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound 
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed. This algorithm is faster than 
    exact_decomposition_dynamic in practical, but asymptotical slower.

    INPUTS:

    - V : a list of vertices

    - E : a list of edges

    - lb : a lower bound to the treewidth of the given graph, 
           e.g. computed by lower_bound (default: '-1')

    OUTPUT:

    - A tuple (V, E, width), where (V, E) is a treedecomposition G of tw(G), 
      if the given lower bound was not greater than tw(G), otherwise a 
      treedecomposition of width 'lb'. 

    EXAMPLES:

        V, E, width = sage.graphs.tdlib.exact_decomposition_cutset(V, E)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = lb 

    gc_exact_decomposition_cutset(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, get_width(V_T, E_T)

def exact_decomposition_dynamic(V, E, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound 
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed.

    INPUTS:

    - V : a list of vertices

    - E : a list of edges

    - lb : a lower bound to the treewidth of the given graph, 
           e.g. computed by lower_bound (default: '-1')

    OUTPUT:

    - A tuple (V, E, width), where (V, E) is a treedecomposition G of tw(G), 
      if the given lower bound was not greater than tw(G), otherwise a 
      treedecomposition of width 'lb'. 

    EXAMPLES:

        V, E, width = sage.graphs.tdlib.exact_decomposition_dynamic(V, E)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = lb 

    gc_exact_decomposition_dynamic(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, get_width(V_T, E_T)

##############################################################
############ APPROXIMATIVE ALGORITHMS ########################

def seperator_algorithm(V, E):
    """
    Computes a tree decomposition of a given graph using nearly balanced 
    seperators. The returned width is at most 4*tw(G)+1.

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - A tuple (V, E, width)

    EXAMPLES:

        V, E, width = tdlib.seperator_algorithm(V, E)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    gc_seperator_algorithm(V_G, E_G, V_T, E_T);

    return V_T, E_T, get_width(V_T, E_T)

def minDegree_ordering(V, E):
    """
    Computes an elimination ordering of a given graph based on the minDegree 
    heuristic.

    INPUT:

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - An elimination ordering on (V, E) according to the minDegree-heuristic.

    EXAMPLES:

        O = sage.graphs.tdlib.minDegree_ordering(V, E)
    """

    cdef vector[unsigned int] V_G, E_G, elim_ordering
    cython_make_tdlib_graph(V, E, V_G, E_G)

    gc_minDegree_ordering(V_G, E_G, elim_ordering);
    
    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i]);

    return py_elim_ordering


def fillIn_ordering(V, E):
    """
    Computes an elimination ordering of a given graph based on the fillIn 
    heuristic.

    INPUT:

    - V : a list of vertices

    - E : a list of edges

    OUTPUT:

    - An elimination ordering on (V, E) according to the fillIn-heuristic.

    EXAMPLES:

        O = sage.graphs.tdlib.fillIn_ordering(V, E)
    """

    cdef vector[unsigned int] V_G, E_G, elim_ordering
    cython_make_tdlib_graph(V, E, V_G, E_G)

    gc_fillIn_ordering(V_G, E_G, elim_ordering)
    
    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i])

    return py_elim_ordering


##############################################################
############ MISC ############################################

def ordering_to_treedec(V, E, O):
    """
    Applies an elimination ordering on a graph and returns 
    the resulting tree decomposition.

    INPUTS:

    - V : a list of vertices

    - E : a list of edges

    - O : an elimination ordering

    OUTPUT:

    - A tree decomposition, that has been made by applying O on G.

    EXAMPLES:

        O = sage.graphs.tdlib.fillIn_ordering(V, E)
        V, E = tdlib.ordering_to_treedec(V, E, O)
    """

    cdef vector[unsigned int] V_G, E_G, E_T, elim_ordering
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    for i in range(0, len(O)):
        elim_ordering.push_back(O[i])

    gc_ordering_to_treedec(V_G, E_G, V_T, E_T, elim_ordering)

    return V_T, E_T, get_width(V_T, E_T)

def treedec_to_ordering(V, E):
    cdef vector[unsigned int] E_T, elim_ordering
    cdef vector[vector[int]] V_T

    cython_make_tdlib_decomp(V, E, V_T, E_T)

    gc_treedec_to_ordering(V_T, E_T, elim_ordering)

    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i])

    return py_elim_ordering

##############################################################
############ POSTPROCESSING ##################################

def MSVS(pyV_G, pyE_G, pyV_T, pyE_T):
    """
    This may reduce the maximal bag size of a tree decomposition by refinement
    of the bags with help of minimal seperating vertex sets. 

    INPUTS:

        - G : a graph

        - T : a tree decomposition of G

    OUTPUT:

        - A tree decomposition of G with possibly smaller width than T

    EXAMPLES:

        T = tdlib.trivial_decomposition(G)
        T_ = sage.graphs.tdlib.MSVS(G, T)
        MSVS reduced the width by 12, new width: 2
    """
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(pyV_G, pyE_G, V_G, E_G)
    cython_make_tdlib_decomp(pyV_T, pyE_T, V_T, E_T)

    old_width = get_width(pyV_T, pyE_T)

    gc_MSVS(V_G, E_G, V_T, E_T)

    new_width = get_width(V_T, E_T)

    return V_T, E_T, new_width


def minimalChordal(V, E, O):
    """
    Returns an alternative elimination E' ordering than the given elimination 
    ordering E, which may cause a lower width than E, when applied to the
    input graph for computing a tree decomposition.

    INPUTS:

    - ``G`` -- a generic graph

    - ``O`` -- an elimination ordering on ``G``

    OUTPUT:

    - An elimination ordering on ``G`` that may cause a lower width of the tree decomposition, that can be made out of
      it, than the width, that ``O`` will cause.

EXAMPLES::

        sage: g = graphs.HigmanSimsGraph()
        sage: g.show()
        sage: o1 = g.vertices()
        sage: t1 = sage.graphs.tdlib.ordering_to_treedec(g, o1)
        sage: t1.show(vertex_size=100)
        sage: o2 = sage.graphs.tdlib.minimalChordal(g, o1)
        sage: t2 = sage.graphs.tdlib.ordering_to_treedec(g, o2)
        sage: t2.show(vertex_size=100)
    """
    cdef vector[unsigned int] V_G, E_G, old_elim_ordering, new_elim_ordering
    cython_make_tdlib_graph(V, E, V_G, E_G)

    for l in range(0, len(O)):
        old_elim_ordering.push_back(O[l])

    gc_minimalChordal(V_G, E_G, old_elim_ordering, new_elim_ordering)
    
    py_new_elim_ordering = []
    cdef i;
    for i in range(0, len(new_elim_ordering)):
        py_new_elim_ordering.append(new_elim_ordering[i])

    return py_new_elim_ordering


##############################################################
############ MISC ############################################

def is_valid_treedecomposition(pyV_G, pyE_G, pyV_T, pyE_T, message=True):
    """
    Checks, if the definition of a tree decomposition holds for
    a tree decomposition and a graph.

    INPUTS:

    - G : a generic graph

    - T : a tree decomposition

    EXAMPLES:

        V_T, E_T, lb = sage.graphs.tdlib.seperator_algorithm(V_G, E_G)
        status = tdlib.is_valid_treedecomposition(V_G, E_G, V_T, E_T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(pyV_G, pyE_G, V_G, E_G)
    cython_make_tdlib_decomp(pyV_T, pyE_T, V_T, E_T)

    cdef c_status;
    c_status = gc_is_valid_treedecomposition(V_G, E_G, V_T, E_T);
    
    py_status = c_status

    if(message):
        if(py_status == 0):
            pass
        elif(py_status == -1):
            print("Invalid tree decomposition: tree decomposition is not a tree")
        elif(py_status == -2):
            print("Invalid tree decomposition: not all vertices covered")
        elif(py_status == -3):
            print("Invalid tree decomposition: not all edges covered")
        elif(py_status == -4):
            print("Invalid tree decomposition: some encoded vertices are not connected in the tree decomposition")
        else:
            pass

    return py_status

def trivial_decomposition(V, E):
    """
    Returns a trivial tree decomposition of the given graph.

    INPUT:

    - G : a graph

    OUTPUT:

    - A trivial tree decomposition of ``G``

    EXAMPLES:

        sage: g = graphs.RandomGNP(10, 0.05)
        sage: t = sage.graphs.tdlib.trivial_decomposition(g)
        sage: t
        Treedecomposition of width 9 on 1 vertices
        sage: t.show(vertex_size=5000)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)
    
    gc_trivial_decomposition(V_G, E_G, V_T, E_T)

    return V_T, E_T

def get_width(V, E):
    """
    Returns the width of a given tree decomposition.

    INPUT:

    - T -- a tree decomposition

    OUTPUT:

    - The width of T

EXAMPLES::

        sage: g = graphs.RandomGNP(10, 0.05)
        sage: t = sage.graphs.tdlib.trivial_decomposition(g)
        sage: sage.graphs.tdlib.get_width(t)
    """
    return gc_get_width(V)

