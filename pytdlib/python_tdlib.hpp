/* TdLib interface */


/* PREPROCESSING */

int gc_preprocessing(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                     std::vector<std::vector<int> > &bags, int lb, unsigned graphtype);

int gc_PP_MD(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb, unsigned graphtype);

int gc_PP_FI(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb, unsigned graphtype);

int gc_PP_FI_TM(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb, unsigned graphtype);

/* LOWER BOUNDS */

int gc_deltaC_min_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_deltaC_max_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_deltaC_least_c(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_LBN_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_LBNC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_LBP_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);

int gc_LBPC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype);


/* EXACT TREE DECOMPOSITIONS */

int gc_exact_decomposition_cutset(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb, unsigned graphtype);

int gc_exact_decomposition_cutset_decision(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int k, unsigned graphtype);

int gc_exact_decomposition_dynamic(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                   std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb, unsigned graphtype);


/* APPROXIMATIVE TREE DECOMPOSITIONS */

int gc_seperator_algorithm(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

int gc_minDegree_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                         std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

int gc_boost_minDegree_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                         std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);


int gc_fillIn_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                      std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

void gc_minDegree_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                           std::vector<unsigned int> &elim_ordering, unsigned graphtype);

void gc_fillIn_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                        std::vector<unsigned int> &elim_ordering, unsigned graphtype);


/* POSTPROCESSING */

int gc_MSVS(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
            std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

void gc_minimalChordal(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                       std::vector<unsigned int> &old_elimination_ordering,
                       std::vector<unsigned int> &new_elimination_ordering, unsigned graphtype);


/* APPLICATIONS */

void gc_max_clique_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                          std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                          std::vector<unsigned int> &C, unsigned graphtype);

void gc_max_independent_set_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                   std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                   std::vector<unsigned int> &IS, unsigned graphtype);

void gc_min_vertex_cover_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                std::vector<unsigned int> &VC, unsigned graphtype);

void gc_min_dominating_set_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                  std::vector<unsigned int> &DS, unsigned graphtype);

void gc_min_coloring_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                            std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                            std::vector<std::vector<int> > &col, unsigned graphtype);


/* MISC */

int gc_ordering_to_treedec(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                           std::vector<unsigned int> &elim_ordering, unsigned graphtype);

void gc_treedec_to_ordering(std::vector<std::vector<int> > &V, std::vector<unsigned int> &E,
                            std::vector<unsigned int> &elim_ordering);

int gc_trivial_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

int gc_validate_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, unsigned graphtype);

int gc_get_width(std::vector<std::vector<int> > &V_T);


/* Generic elimination search */

void gc_generic_elimination_search1(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings);
void gc_generic_elimination_search2(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings);
void gc_generic_elimination_search3(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings);
void gc_generic_elimination_search4(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, unsigned graphtype, unsigned max_nodes, unsigned max_orderings);

