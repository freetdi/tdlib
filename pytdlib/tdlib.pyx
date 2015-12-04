r"""
Interface with TdLib (algorithms for tree decompositions)

This module defines functions based on TdLib, a
library that implements algorithms for tree
decompositions written by Lukas Larisch. 

Definition:

A tree decomposition of a graph G is a pair (T, b) consisting of a graph T and 
a function b: V(T) -> 2^{V(G)} associating with each vertex t in V(T) a set of 
vertices b(t), called bags of T, that are subsets of V(G) such that 

(T1) T is a tree,

(T2) every vertex v in V(G) is contained in some bag of T,

(T3) for every edge e in E(G) there is a node t in V(T) with e is a subset of 
b(t), and

(T4) for all v in V(G) the set b^{-1} := {t in V(T): v in b(t)} is non-empty 
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

    - treedecomposition_exact_cutset_decision
        Computes a tree decomposition of exact width of a given graph
        and returns true, if tw(G) <= k. Otherwise false will be returned

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
from tdlib cimport gc_exact_decomposition_cutset_decision
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
        if isinstance(pyE[0], tuple):
            for u,v in pyE:
                E.push_back(u)
                E.push_back(v)

        #list representation
        elif isinstance(pyE[0], list):
            for e in pyE:
                E.push_back(e[0])
                E.push_back(e[1])

        #internal representation (unfolded tuple/list representation)
        elif isinstance(pyE[0], int):
            for e in pyE:
                E.push_back(e)

cdef cython_make_tdlib_decomp(pyV, pyE, vector[vector[int]] &V, vector[unsigned int] &E):
    for t in pyV:
        V.push_back(t)

    #tuple representation
    if isinstance(pyE[0], tuple):
        for u,v in pyE:
            E.push_back(u)
            E.push_back(v)

    #list representation
    elif isinstance(pyE[0], list):
        for e in pyE:
            E.push_back(e[0])
            E.push_back(e[1])

    #internal representation (unfolded tuple/list representation)
    elif isinstance(pyE[0], int):
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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - V_G : a list of vertices of the reduced input graph

    - E_G : a list of edges of the reduced input graph

    - bags

    - lb : a lower bounds on treewidth of (V_G, E_G)

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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, width = tdlib.PP_MD(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = -1

    gc_PP_MD(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, get_width(V_T, E_T)

def PP_FI_TM(V, E):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced 
    instance G' of G will be processed by the fillIn heuristic, which
    successivly eliminates a vertex, that will cause least new edges
    within the elimination process. The resulting tree decomposition 
    will be postprocessed by the minimalChordal algorithm, that may
    reduce the width of the tree decomposition.

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, width = tdlib.PP_FI_TM(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = -1

    gc_PP_FI_TM(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, get_width(V_T, E_T)


def preprocessing_glue_bags(V, E, bags):
    """
    Glues a list of bags returned by preprocessing with T according to the 
    subset-relation. This can be used to make a tree decomposition out of the 
    results of tdlib.preprocessing() for a given graph G.

    INPUTS:

    - V_T : a list of vertices of the input treedecomposition

    - E_T : a list of edges of the input treedecomposition

    - bags : a list of bags, e.g. as returned by preprocessing

    OUTPUTS:

    - V_T' : a list of vertices of the input treedecomposition

    - E_T' : a list of edges of the input treedecomposition

    EXAMPLES:

        V_G2, E_G2, bags, lb = tdlib.preprocessing(V_G1, E_G1)
        ...
        V_T1, E_T1 = ...
        ...
        V_T2, E_T2 = tdlib.preprocessing_glue_bags(V_T1, E_T1, bags)
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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - algorithm -- (default: 'deltaC_least_c') specifies the algorithm to use 
                   for computing a lower bound on the treewidth of G. The 
                   algorithms have to be choosen among the above mentioned. 

    OUTPUT:

    - lb : a lower bound on the treewidth of (V_G, E_G)

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

    return c_lb


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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - lb : a lower bound to the treewidth of (V_G, E_G),
           e.g. computed by lower_bound (default: '-1')

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, width = tdlib.exact_decomposition_cutset(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_lb = lb 

    gc_exact_decomposition_cutset(V_G, E_G, V_T, E_T, c_lb)

    return V_T, E_T, get_width(V_T, E_T)

def exact_decomposition_cutset_decision(V, E, k):
    """
    Computes a tree decomposition of exact width, if tw(G)  k. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed. This algorithm is faster than 
    exact_decomposition_dynamic in practical, but asymptotical slower.

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - lb : a lower bound to the treewidth of (V_G, E_G),
           e.g. computed by lower_bound (default: '-1')

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - status : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, status = tdlib.exact_decomposition_cutset(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)

    cdef int c_k = k 

    is_leq = gc_exact_decomposition_cutset_decision(V_G, E_G, V_T, E_T, c_k)
    
    if is_leq is 0:
        rtn = True
    else:
        rtn = False

    return V_T, E_T, rtn

def exact_decomposition_dynamic(V, E, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound 
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed.

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - lb : a lower bound to the treewidth of (V_G, E_G),
           e.g. computed by lower_bound (default: '-1')

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, width = tdlib.exact_decomposition_dynamic(V_G, E_G)
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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T, width = tdlib.seperator_algorithm(V_G, E_G)
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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - O : an elimination ordering on (V_G, E_G)

    EXAMPLES:

        O = tdlib.minDegree_ordering(V_G, E_G)
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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - O : an elimination ordering on (V_G, E_G)

    EXAMPLES:

        O = tdlib.fillIn_ordering(V_G, E_G)
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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - O : an elimination ordering on (V_G, E_G)

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    - width : the width of (V_T, E_T)

    EXAMPLES:

        O = tdlib.fillIn_ordering(V_G, E_G)
        V_T, E_T, width = tdlib.ordering_to_treedec(V_G, E_G, O)
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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - V_T : a list of vertices of the input treedecomposition

    - E_T : a list of edges of the input treedecomposition

    OUTPUTS:

    - V_T' : a list of vertices of a treedecomposition

    - E_T' : a list of edges of a treedecomposition

    - width : the width of (V_T', E_T')

    EXAMPLES:

        V_T1, E_T1 = tdlib.trivial_decomposition(V_G, E_G)
        V_T2, E_T2, width = tdlib.MSVS(V_G, E_G, V_T1, E_T1)
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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - O : an elimination ordering on (V_G, E_G)

    OUTPUT:

    - An elimination ordering on G that may cause a lower width of the 
      treedecomposition, that can be made out of it, than the width, 
      that O will cause.

    EXAMPLES:

        V_T1, E_T1, w1 = tdlib.ordering_to_treedec(V_G, E_G, O1)
        O2 = tdlib.minimalChordal(V_G, E_G, O1)
        V_T2, E_T2, w2 = tdlib.ordering_to_treedec(V_G, E_G, O2)
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

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    - V_T : a list of vertices of the input treedecomposition

    - E_T : a list of edges of the input treedecomposition

    - message : outputs error message iff (V_T, E_T) is invalid with
                respect to (V_G, E_G) (optional)

    OUTPUT:

    - error_code :     0, if (V_T, E_T) is valid with respect to (V_G, E_G)
                   <= -1, if (V_T, E_T) is not a tree
                   <= -2, if not all vertices of (V_G, E_G) are covered
                   <= -3, if not all edges of (V_G, E_G) are covered
                   == -4, if condition (T4) of a treedecomposition is 
                             not satisfied

    EXAMPLES:

        V_T, E_T, lb = tdlib.seperator_algorithm(V_G, E_G)
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

    INPUTS:

    - V_G : a list of vertices of the input graph

    - E_G : a list of edges of the input graph

    OUTPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    EXAMPLES:

        V_T, E_T = tdlib.trivial_decomposition(V_G, E_G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(V, E, V_G, E_G)
    
    gc_trivial_decomposition(V_G, E_G, V_T, E_T)

    return V_T, E_T

def get_width(V, E):
    """
    Returns the width of a given tree decomposition.

    INPUTS:

    - V_T : a list of vertices of a treedecomposition

    - E_T : a list of edges of a treedecomposition

    OUTPUT:

    - width : the width of (V_T, E_T)

    EXAMPLES:

        V_T, E_T = tdlib.trivial_decomposition(V_G, E_G)
        width = tdlib.get_width(V_T, E_T)
    """

    return gc_get_width(V)

