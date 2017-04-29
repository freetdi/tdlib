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

    - PP_FI
        Applies preprocessing followed by the minFill-heuristic to a given
        graph

    - PP_FI_TM
        Applies preprocessing followed by the minFill-heuristic followed by
        triangulation minimization to a given graph

    - lower_bound
        Computes a lower bound with respect to treewidth of a given graph

    - exact_decomposition_cutset
        Computes a tree decomposition of exact width of a given graph
        (faster than the dynamic version in practice)

    - exact_decomposition_cutset_decision
        Computes a tree decomposition of exact width of a given graph
        and returns true, if tw(G) <= k. Otherwise false will be returned

    - exact_decomposition_dynamic
        Computes a tree decomposition of exact width of a given graph
        (asymptotically faster than the greedy version)

    #- treedecomposition_exact_branch_and_bound
        Computes a tree decomposition of exact width of a given graph

    - separator_algorithm
        Computes a 4-approximate tree decomposition of a given graph

    - minDegree_decomp
        Computes a tree decomposition according to the minDegree-heuristic

    - fillIn_decomp
        Computes a tree decomposition according to the minFill-heuristic

    - minDegree_ordering
        Computes an elimination ordering according to the minDegree-heuristic

    - fillIn_ordering
        Computes an elimination ordering according to the minFill-heuristic

    - MSVS
        Possibly reduces the width of a given tree decomposition with help
        of minimal seperating vertex sets

    - minimalChordal_ordering
        Possibly reduces the width of a given tree decomposition by
        triangulation minimization. Input: elimination ordering

    - minimalChordal_decomp
        Possibly reduces the width of a given tree decomposition by
        triangulation minimization. Input: tree decomposition

    #- random_branch_and_bound

    #- random_elimination_orderings

    - max_clique_with_treedecomposition
        Computes a maximal clique with help of a tree decomposition

    - max_independent_set_with_treedecomposition
        Computes a maximal independent set with help of a tree decomposition

    - max_vertex_cover_with_treedecomposition
        Computes a minimal vertex cover with help of a tree decomposition

    - min_dominating_set_with_treedecomposition
        Computes a minimal dominating set with help of a tree decomposition

    - min_coloring_with_treedecomposition
        Computes a minimal coloring with help of a tree decomposition

    - treedec_to_ordering
        Computes an elimination ordering out of a tree decomposition

    - ordering_to_treedec
        Computes a tree decomposition out of an elimination ordering

    - trivial_decomposition
        Returns a trivial decomposition of a given graph

    - is_valid_treedecomposition
        Checks, if a tree decomposition is valid with respect to a given graph

    - get_width
        Returns the width of a given tree decomposition

    - generic_elimination_search1 TODO: cleanup
        ....

    - generic_elimination_search2
        ....

    - generic_elimination_search3
        ....

    - generic_elimination_search4
        ....

AUTHOR: Lukas Larisch (now): Initial version
-------
"""

from libcpp.vector cimport vector

from Graph import Graph
from Decomp import Decomp


##############################################################
############ GRAPH/DECOMPOSITION ENCODING/DECODING ###########
#the following will be used implicitly do the translation
#between the python graph encoding and TdLib graph encoding,
#which is based on the BGL.

cdef cython_make_tdlib_graph(pyV, pyE, vector[unsigned int] &V, vector[unsigned int] &E):
    labels_map = list()
    labels_dict_inv = dict()

    for i in range(0, len(pyV)):
        V.push_back(i)
        labels_map.append(pyV[i])
        labels_dict_inv[pyV[i]] = i

    if len(pyE) > 0:
        #tuple representation
        if isinstance(pyE[0], tuple):
            for u,v in pyE:
                E.push_back(labels_dict_inv[u])
                E.push_back(labels_dict_inv[v])

        #list representation
        elif isinstance(pyE[0], list):
            for e in pyE:
                E.push_back(labels_dict_inv[e[0]])
                E.push_back(labels_dict_inv[e[1]])

        #internal representation (unfolded tuple/list representation)
        elif isinstance(pyE[0], int):
            for e in pyE:
                E.push_back(labels_dict_inv[e])

    return labels_map


cdef cython_make_tdlib_decomp(pyV, pyE, vector[vector[int]] &V, vector[unsigned int] &E, inv_labels_dict=dict()):
    labels_dict = dict()
    labels_map = list()

    if(len(inv_labels_dict) == 0):
        i = int(0)
        for bag in pyV:
            for v in bag:
                if v not in labels_dict:
                    labels_dict[v] = i
                    labels_map.append(v)
                    i += 1
    else:
        labels_dict = inv_labels_dict

    try:
        for bag in pyV:
            bag_ = []
            for v in bag:
                bag_.append(labels_dict[v])
            V.push_back(bag_)
    except KeyError:
        return False

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

    return labels_map


def apply_labeling(X, labels_map):
    X_ = list()
    if len(X) > 0:
        #X is a list of vertices/edges
        if isinstance(X[0], int):
            for x in X:
                X_.append(labels_map[x])
        #X is a list of bags
        if isinstance(X[0], list):
            for x in X:
                 Y = list()
                 for x_ in x:
                     Y.append(labels_map[x_])
                 X_.append(Y)

    return X_


def inverse_labels_dict(labels_map):
    inv_dict = dict()
    for i in range(0, len(labels_map)):
        inv_dict[labels_map[i]] = i
    return inv_dict

#do not change without modifying python_tdlib.cpp
def graphtype_to_uint(string):
    if string == "boost_graph_undirected":
        return 0
    elif string == "boost_graph_directed":
        return 1
    else:
        raise Exception


##############################################################
############ PREPROCESSING ###################################

def preprocessing(G):
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

    - G : input graph

    OUTPUTS:

    - G' : reduced input graph

    - bags

    - lb : a lower bounds on treewidth of G'

    EXAMPLES:

        G_, bags, lb = tdlib.preprocessing(G)
    """

    cdef vector[unsigned int] V_G, E_G
    cdef vector[vector[int]] c_bags;

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = -1

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    py_lb = gc_preprocessing(V_G, E_G, c_bags, c_lb, graphtype)

    V_G_ = apply_labeling(V_G, labels_map)
    c_bags_ = apply_labeling(c_bags, labels_map)

    return Graph(V_G_, E_G), c_bags_, py_lb


def PP_MD(G):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced
    instance G' of G will be processed by the minDegree heuristic, which
    successivly eliminates a vertex of minimal degree. The returned tree
    decomposition then may be of non-optimal width.

    INPUTS:

    - G : input graph (must provide the methods vertices() and edges())

    OUTPUTS:

    - T : treedecomposition

    - width : the width of T

    EXAMPLES:

        G = Graph([1,2,3], [[1,2],[2,3]])
        T, width = tdlib.PP_MD(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = -1

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_PP_MD(V_G, E_G, V_T, E_T, c_lb, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def PP_FI(G):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced
    instance G' of G will be processed by the fillIn heuristic, which
    successivly eliminates a vertex, that will cause least new edges
    within the elimination process.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.PP_FI(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = -1

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_PP_FI(V_G, E_G, V_T, E_T, c_lb, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def PP_FI_TM(G):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced
    instance G' of G will be processed by the fillIn heuristic, which
    successivly eliminates a vertex, that will cause least new edges
    within the elimination process. The resulting tree decomposition
    will be postprocessed by the minimalChordal algorithm, that may
    reduce the width of the tree decomposition.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.PP_FI_TM(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = -1

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_PP_FI_TM(V_G, E_G, V_T, E_T, c_lb, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


##############################################################
############ LOWER BOUNDS ####################################

def lower_bound(G, algorithm = "deltaC_least_c"):
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

    - G : input graph

    - algorithm -- (default: 'deltaC_least_c') specifies the algorithm to use
                   for computing a lower bound on the treewidth of G. The
                   algorithms have to be choosen among the above mentioned.

    OUTPUT:

    - lb : a lower bound on the treewidth of G

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
    cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    cdef int c_lb = 0

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    if(algorithm == "deltaC_min_d"):
        c_lb = gc_deltaC_min_d(V_G, E_G, graphtype)
    elif(algorithm == "deltaC_max_d"):
        c_lb = gc_deltaC_max_d(V_G, E_G, graphtype)
    elif(algorithm == "deltaC_least_c"):
        c_lb = gc_deltaC_least_c(V_G, E_G, graphtype)
    elif(algorithm == "LBN_deltaC"):
        c_lb = gc_LBN_deltaC(V_G, E_G, graphtype)
    elif(algorithm == "LBNC_deltaC"):
        c_lb = gc_LBNC_deltaC(V_G, E_G, graphtype)
    elif(algorithm == "LBP_deltaC"):
        c_lb = gc_LBP_deltaC(V_G, E_G, graphtype)
    elif(algorithm == "LBPC_deltaC"):
        c_lb = gc_LBPC_deltaC(V_G, E_G, graphtype)
    else:
        print("Invalid lower bound algorithm")
        return -2

    return c_lb


##############################################################
############ EXACT ALGORITHMS ################################

def exact_decomposition_cutset(G, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed. This algorithm is faster than
    exact_decomposition_dynamic in practical, but asymptotical slower.

    INPUTS:

    - G : input graph

    - lb : a lower bound to the treewidth of G,
           e.g. computed by lower_bound (default: '-1')

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.exact_decomposition_cutset(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = lb

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_exact_decomposition_cutset(V_G, E_G, V_T, E_T, c_lb, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def exact_decomposition_cutset_decision(G, k):
    """
    Computes a tree decomposition of exact width, if tw(G)  k. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed. This algorithm is faster than
    exact_decomposition_dynamic in practical, but asymptotical slower.

    INPUTS:

    - G : input graph

    - k : parameter for the question 'tw(G) < k'

    OUTPUTS:

    - status : True if tw(G) < k, else False.

    EXAMPLES:

        status = tdlib.exact_decomposition_cutset_decision(G, 3)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_k = k

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    is_leq = gc_exact_decomposition_cutset_decision(V_G, E_G, V_T, E_T, c_k, graphtype)

    if is_leq is 0:
        rtn = True
    else:
        rtn = False

    return rtn

def exact_decomposition_dynamic(G, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed.

    INPUTS:

    - G : input graph

    - lb : a lower bound to the treewidth of G,
           e.g. computed by lower_bound (default: '-1')

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.exact_decomposition_dynamic(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef int c_lb = lb

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_exact_decomposition_dynamic(V_G, E_G, V_T, E_T, c_lb, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


##############################################################
############ APPROXIMATIVE ALGORITHMS ########################

def seperator_algorithm(G):
    """
    Computes a tree decomposition of a given graph using nearly balanced
    seperators. The returned width is at most 4*tw(G)+1.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.seperator_algorithm(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_seperator_algorithm(V_G, E_G, V_T, E_T, graphtype);

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def minDegree_decomp(G):
    """
    Computes a tree decomposition of a given graph based on the minDegree heuristic.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.minDegree_decomp(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_minDegree_decomp(V_G, E_G, V_T, E_T, graphtype);

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)

def boost_minDegree_decomp(G):
    """
    Computes a tree decomposition of a given graph based on the (boost-)minDegree heuristic.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.boost_minDegree_decomp(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_boost_minDegree_decomp(V_G, E_G, V_T, E_T, graphtype);

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def fillIn_decomp(G):
    """
    Computes a tree decomposition of a given graph based on the fillIn heuristic.
     INPUTS:

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, width = tdlib.fillIn_decomp(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_fillIn_decomp(V_G, E_G, V_T, E_T, graphtype);

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def minDegree_ordering(G):
    """
    Computes an elimination ordering of a given graph based on the minDegree
    heuristic.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - O : an elimination ordering on (V_G, E_G)

    EXAMPLES:

        O = tdlib.minDegree_ordering(G)
    """

    cdef vector[unsigned int] V_G, E_G, elim_ordering
    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_minDegree_ordering(V_G, E_G, elim_ordering, graphtype);

    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i]);

    elim_ordering_ = apply_labeling(py_elim_ordering, labels_map)

    return elim_ordering_


def fillIn_ordering(G):
    """
    Computes an elimination ordering of a given graph based on the fillIn
    heuristic.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - O : an elimination ordering on (G)

    EXAMPLES:

        O = tdlib.fillIn_ordering(G)
    """

    cdef vector[unsigned int] V_G, E_G, elim_ordering
    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_fillIn_ordering(V_G, E_G, elim_ordering, graphtype)

    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i])

    elim_ordering_ = apply_labeling(py_elim_ordering, labels_map)

    return elim_ordering_


##############################################################
############ POSTPROCESSING ##################################

def MSVS(G, T):
    """
    This may reduce the maximal bag size of a tree decomposition by refinement
    of the bags with help of minimal seperating vertex sets.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUTS:

    - T' : a treedecomposition of G

    - width : the width of T'

    EXAMPLES:

        T1 = tdlib.trivial_decomposition(G)
        T2, width = tdlib.MSVS(G, T1)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_MSVS(V_G, E_G, V_T, E_T, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T_ = Decomp(V_T_, E_T)

    return T_, get_width(T_)


def minimalChordal_ordering(G, O):
    """
    Returns an alternative elimination ordering E' than the given elimination
    ordering E, which may cause a lower width than E, when applied to the
    input graph for computing a tree decomposition.

    INPUTS:

    - G : input graph

    - O : an elimination ordering on G

    OUTPUT:

    - An elimination ordering on G that may cause a lower width of the
      treedecomposition, that can be made out of it, than the width,
      that O will cause.

    EXAMPLES:

        O1 = range(0, len(V))
        O2 = tdlib.minimalChordal_ordering(G, O1)
        T, w = tdlib.ordering_to_treedec(G, O2)
    """

    cdef vector[unsigned int] V_G, E_G, old_elim_ordering, new_elim_ordering
    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)

    for l in range(0, len(O)):
        old_elim_ordering.push_back(inv_labels_dict[O[l]])

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_minimalChordal(V_G, E_G, old_elim_ordering, new_elim_ordering, graphtype)

    py_new_elim_ordering_ = []
    cdef i;
    for i in range(0, len(new_elim_ordering)):
        py_new_elim_ordering_.append(new_elim_ordering[i])

    py_new_elim_ordering = []
    for i in range(0, len(py_new_elim_ordering_)):
        py_new_elim_ordering.append(labels_map[py_new_elim_ordering_[i]])

    return py_new_elim_ordering


def minimalChordal_decomp(G, T):
    """
    Returns an alternative elimination ordering E' than the given elimination
    ordering E, which may cause a lower width than E, when applied to the
    input graph for computing a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUT:

    - T' : a tree decomposition of G of possibly lower width than the width of T

    EXAMPLES:

        T1, w1 = tdlib.minDegree_decomp(G)
        T2, w2 = tdlib.minimalChordal_decomp(G, T1)
    """

    O1 = treedec_to_ordering(T)
    O2 = minimalChordal_ordering(G, O1)

    return ordering_to_treedec(G, O2)


##############################################################
############ APPLICATIONS ####################################

def max_clique_with_treedecomposition(G, T):
    """
    Computes a maximum clique with help of a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUT:

    - C:    a maximum clique in G

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        C = tdlib.max_clique_with_treedecomposition(G, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T, C_
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_max_clique_with_treedecomposition(V_G, E_G, V_T, E_T, C_, graphtype);

    py_C = []
    cdef i;
    for i in range(0, len(C_)):
        pyCi = C_[i]
        py_C.append(pyCi)

    return py_C


def max_independent_set_with_treedecomposition(G, T):
    """
    Computes a maximum sized independent set with help of a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUT:

    - IS:    a maximum sized independent set in G

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        IS = tdlib.max_independent_set_with_treedecomposition(V, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T, IS_
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_max_independent_set_with_treedecomposition(V_G, E_G, V_T, E_T, IS_, graphtype);

    py_IS = []
    cdef i;
    for i in range(0, len(IS_)):
        pyISi = IS_[i]
        py_IS.append(pyISi)

    return py_IS


def min_vertex_cover_with_treedecomposition(G, T):
    """
    Computes a minimum vertex cover with help of a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUT:

    - VC:    a minimal vertex cover in G

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        VC = tdlib.min_vertex_cover_with_treedecomposition(G, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T, VC_
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_min_vertex_cover_with_treedecomposition(V_G, E_G, V_T, E_T, VC_, graphtype);

    py_VC = []
    cdef i;
    for i in range(0, len(VC_)):
        pyVCi = VC_[i]
        py_VC.append(pyVCi)

    return py_VC


def min_dominating_set_with_treedecomposition(G, T):
    """
    Computes a minimal dominating set based on a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUTS:

    - DS : a list of vertices of a minimal dominating set in G

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        DS = tdlib.min_dominating_set_with_treedecomposition(G, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T, DS
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_min_dominating_set_with_treedecomposition(V_G, E_G, V_T, E_T, DS, graphtype)

    py_DS = []
    cdef i;
    for i in range(0, len(DS)):
        py_DS.append(DS[i])

    return py_DS


def min_coloring_with_treedecomposition(G, T):
    """
    Computes a minimum coloring with help of a tree decomposition.

    INPUTS:

    - G : input graph

    - T : a treedecomposition of G

    OUTPUT:

    - VC:    a minimal coloring of G

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        VC = tdlib.min_coloring_with_treedecomposition(G, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T, C_

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False):
        return

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_min_coloring_with_treedecomposition(V_G, E_G, V_T, E_T, C_, graphtype);

    py_C = []
    cdef i;
    for i in range(0, len(C_)):
        pyCi = []
        for j in range(0, len(C_[i])):
            pyCij = C_[i][j]
            pyCi.append(pyCij)
        py_C.append(pyCi)

    return py_C


##############################################################
############ MISC ############################################

def ordering_to_treedec(G, O):
    """
    Applies an elimination ordering on a graph and returns
    the resulting tree decomposition.

    INPUTS:

    - G : input graph

    - O : an elimination ordering on G

    OUTPUTS:

    - T : treedecomposition of G

    - width : the width of T

    EXAMPLES:

        O = tdlib.fillIn_ordering(G)
        T, width = tdlib.ordering_to_treedec(G)
    """

    if(sorted(G.vertices()) != sorted(O)):
        print("error: an elimination ordering must be a permutation of the vertices!")
        return

    cdef vector[unsigned int] V_G, E_G, E_T, elim_ordering
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    inv_labels_dict = inverse_labels_dict(labels_map)

    for i in range(0, len(O)):
        elim_ordering.push_back(inv_labels_dict[O[i]])

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_ordering_to_treedec(V_G, E_G, V_T, E_T, elim_ordering, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def treedec_to_ordering(T):
    """
    Converts a treedecomposition to an elimination ordering.

    INPUTS:

    - T : a treedecomposition

    OUTPUTS:

    - O : an elimination ordering

    EXAMPLES:

        T = tdlib.fillIn_decomp(G)
        O = tdlib.treedec_to_ordering(T)
    """

    cdef vector[unsigned int] E_T, elim_ordering
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T)

    gc_treedec_to_ordering(V_T, E_T, elim_ordering)

    py_elim_ordering_ = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering_.append(elim_ordering[i])

    py_elim_ordering = []
    for i in range(0, len(py_elim_ordering_)):
        py_elim_ordering.append(labels_map[py_elim_ordering_[i]])

    return py_elim_ordering


def trivial_decomposition(G):
    """
    Returns a trivial tree decomposition of the given graph.

    INPUTS:

    - G : input graph

    OUTPUTS:

    - T : a treedecomposition of G

    - width : the width of T

    EXAMPLES:

        T, w = tdlib.trivial_decomposition(G)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    gc_trivial_decomposition(V_G, E_G, V_T, E_T, graphtype)

    V_T_ = apply_labeling(V_T, labels_map)

    T = Decomp(V_T_, E_T)

    return T, get_width(T)


def is_valid_treedecomposition(G, T, message=True):
    """
    Checks, if the definition of a tree decomposition holds for
    a tree decomposition and a graph.

    INPUTS:

    - G : input graph

    - T : treedecomposition

    - message : outputs error message iff T is invalid with
                respect to G (optional)

    OUTPUT:

    - status : True if T is a valid treedecomposition of G else False.

    EXAMPLES:

        T, w = tdlib.seperator_algorithm(G)
        status = tdlib.is_valid_treedecomposition(G, T)
    """

    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)
    inv_labels_dict = inverse_labels_dict(labels_map)
    rtn = cython_make_tdlib_decomp(T.vertices(), T.edges(), V_T, E_T, inv_labels_dict)

    if(rtn is False and T.vertices() != [[]] and T.edges() != []):
        if message:
            print("error: labels_dict is corrupted (possible reason: there is no bijective mapping 'bags -> vertices'")
        return False

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    cdef c_status;
    c_status = gc_validate_treedecomposition(V_G, E_G, V_T, E_T, graphtype);

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

    return py_status == 0


def get_width(T):
    """
    Returns the width of a given tree decomposition.

    INPUTS:

    - T : a treedecomposition

    OUTPUT:

    - width : the width of T

    EXAMPLES:

        T = tdlib.trivial_decomposition(G)
        width = tdlib.get_width(T)
    """

    return gc_get_width(T.vertices())


def generic_elimination_search1(G, max_nodes, max_orderings):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    cdef unsigned max_nodes_c = max_nodes
    cdef unsigned max_orderings_c = max_orderings

    gc_generic_elimination_search1(V_G, E_G, graphtype, max_nodes_c, max_orderings_c)


def generic_elimination_search2(G, max_nodes, max_orderings):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    cdef unsigned max_nodes_c = max_nodes
    cdef unsigned max_orderings_c = max_orderings

    gc_generic_elimination_search2(V_G, E_G, graphtype, max_nodes_c, max_orderings_c)


def generic_elimination_search3(G, max_nodes, max_orderings):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    cdef unsigned max_nodes_c = max_nodes
    cdef unsigned max_orderings_c = max_orderings

    gc_generic_elimination_search3(V_G, E_G, graphtype, max_nodes_c, max_orderings_c)


def generic_elimination_search_p17(G, max_nodes, max_orderings):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    labels_map = cython_make_tdlib_graph(G.vertices(), G.edges(), V_G, E_G)

    cdef unsigned graphtype = graphtype_to_uint(G.graphtype())

    cdef unsigned max_nodes_c = max_nodes
    cdef unsigned max_orderings_c = max_orderings

    gc_generic_elimination_search_p17(V_G, E_G, graphtype, max_nodes_c, max_orderings_c)

