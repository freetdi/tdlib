from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "python_tdlib.hpp":

##############################################################
############ PREPROCESSING ###################################

    int gc_preprocessing(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                         vector[vector[int]] &bags, int lb, unsigned graphtype)
    int gc_PP_MD(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                 vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)
    int gc_PP_FI(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                 vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)
    int gc_PP_FI_TM(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                    vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)

##############################################################
############ LOWER BOUNDS ####################################

    int gc_deltaC_min_d(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)
    int gc_deltaC_max_d(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)
    int gc_deltaC_least_c(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)

    int gc_LBN_deltaC(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)
    int gc_LBNC_deltaC(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)
    int gc_LBP_deltaC(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)
    int gc_LBPC_deltaC(vector[unsigned int] &V, vector[unsigned int] &E, unsigned graphtype)

##############################################################
############ EXACT ALGORITHMS ################################

    int gc_exact_decomposition_cutset(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                      vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)
    int gc_exact_decomposition_cutset_decision(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                      vector[vector[int]] &V_T, vector[unsigned int] &E_T, int k, unsigned graphtype)
    int gc_exact_decomposition_dynamic(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                      vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)
    #int gc_exact_decomposition_branch_and_bound(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
    #                                            vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)

##############################################################
############ APPROXIMATIVE ALGORITHMS ########################

    int gc_seperator_algorithm(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    void gc_minDegree_decomp(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                               vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    void gc_boost_minDegree_decomp(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                               vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    void gc_fillIn_decomp(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                            vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    void gc_minDegree_ordering(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                               vector[unsigned int] &elim_ordering, unsigned graphtype)
    void gc_fillIn_ordering(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                            vector[unsigned int] &elim_ordering, unsigned graphtype)
    #int gc_random_branch_and_bound(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
    #                        vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)
    #int gc_random_elimination_orderings(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
    #                                    vector[vector[int]] &V_T, vector[unsigned int] &E_T, int lb, unsigned graphtype)

##############################################################
############ POSTPROCESSING ##################################

    int gc_MSVS(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)

    void gc_minimalChordal(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                           vector[unsigned int] &old_elimination_ordering,
                           vector[unsigned int] &new_elimination_ordering, unsigned graphtype)

##############################################################
############ APPLICATIONS ####################################

    void gc_max_clique_with_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                              vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                                              vector[unsigned int] &C, unsigned graphtype)
    void gc_max_independent_set_with_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                                       vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                                                       vector[unsigned int] &IS, unsigned graphtype)
    void gc_min_vertex_cover_with_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                                    vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                                                    vector[unsigned int] &VC, unsigned graphtype)
    void gc_min_coloring_with_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                                    vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                                                    vector[vector[int]] &col, unsigned graphtype)
    void gc_min_dominating_set_with_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                                       vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                                                       vector[unsigned int] &DS, unsigned graphtype)

##############################################################
############ MISC ############################################

    int gc_ordering_to_treedec(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                               vector[vector[int]] &V_T, vector[unsigned int] &E_T,
                               vector[unsigned int] &elim_ordering, unsigned graphtype)
    void gc_treedec_to_ordering(vector[vector[int]] &V, vector[unsigned int] &E,
                                vector[unsigned int] &elim_ordering);
    int gc_validate_treedecomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                      vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    void gc_trivial_decomposition(vector[unsigned int] &V_G, vector[unsigned int] &E_G,
                                  vector[vector[int]] &V_T, vector[unsigned int] &E_T, unsigned graphtype)
    int gc_get_width(vector[vector[int]] &V_T)


##############################################################
############ GENERIC ELIMINATION SEARCH ######################

    void gc_generic_elimination_search1(vector[unsigned int] &V_G, vector[unsigned int] &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings)
    void gc_generic_elimination_search2(vector[unsigned int] &V_G, vector[unsigned int] &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings)
    void gc_generic_elimination_search3(vector[unsigned int] &V_G, vector[unsigned int] &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings)
