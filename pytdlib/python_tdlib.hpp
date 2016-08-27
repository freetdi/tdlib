/* TdLib interface */


/* PREPROCESSING */

int gc_preprocessing(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                     std::vector<std::vector<int> > &bags, int lb);

int gc_PP_MD(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);

int gc_PP_FI(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);

int gc_PP_FI_TM(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);

/* LOWER BOUNDS */

int gc_deltaC_min_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_deltaC_max_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_deltaC_least_c(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_LBN_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_LBNC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_LBP_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);

int gc_LBPC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G);


/* EXACT TREE DECOMPOSITIONS */

int gc_exact_decomposition_cutset(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);

int gc_exact_decomposition_cutset_decision(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int k);

int gc_exact_decomposition_dynamic(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                   std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb);


/* APPROXIMATIVE TREE DECOMPOSITIONS */

int gc_seperator_algorithm(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

int gc_minDegree_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                         std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

int gc_boost_minDegree_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                         std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);


int gc_fillIn_decomp(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                      std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

void gc_minDegree_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                           std::vector<unsigned int> &elim_ordering);

void gc_fillIn_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                        std::vector<unsigned int> &elim_ordering);


/* POSTPROCESSING */

int gc_MSVS(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
            std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

void gc_minimalChordal(std::vector<unsigned int> &V, std::vector<unsigned int> &E,
                       std::vector<unsigned int> &old_elimination_ordering,
                       std::vector<unsigned int> &new_elimination_ordering);


/* APPLICATIONS */

void gc_max_clique_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                          std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                          std::vector<unsigned int> &C);

void gc_max_independent_set_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                   std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                   std::vector<unsigned int> &IS);

void gc_min_vertex_cover_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                std::vector<unsigned int> &VC);

void gc_min_dominating_set_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                                  std::vector<unsigned int> &DS);

void gc_min_coloring_with_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                            std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                                            std::vector<std::vector<int> > &col);


/* MISC */

int gc_ordering_to_treedec(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                           std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T,
                           std::vector<unsigned int> &elim_ordering);

void gc_treedec_to_ordering(std::vector<std::vector<int> > &V, std::vector<unsigned int> &E,
                            std::vector<unsigned int> &elim_ordering);

int gc_trivial_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                             std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

int gc_is_valid_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G,
                                  std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T);

int gc_get_width(std::vector<std::vector<int> > &V_T);

