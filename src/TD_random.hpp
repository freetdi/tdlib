
#ifndef TD_RANDOM

namespace treedec{

#include "TD_elimination_orderings.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename G_t>
int randomly_try_some_elimination_orderings(G_t &G, unsigned int count = 5){
    std::vector<unsigned int> elim_ordering(boost::num_vertices(G));
    for(unsigned int i=0; i<elim_ordering.size(); ++i){
        elim_ordering[i] = i;
    }

    std::vector<std::vector<unsigned int> > elimination_orderings(count);
    for(unsigned int i=0; i<elimination_orderings.size(); ++i){
        std::random_shuffle(elim_ordering.begin(), elim_ordering.end()); // using built-in random generator
        elimination_orderings[i] = elim_ordering;
    }

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    int min_width = INT_MAX;

    #pragma omp parallel for
    for(unsigned int i = 0; i < count; i++){
        G_t H;
        boost::copy_graph(G, H); // ..(H, G)..?! "unavoidable"?
        int width_i = get_width_of_elimination_ordering(H, elimination_orderings[i], idxMap);
        //std::cout << "width_" << i << ": " << width_i << std::endl;
        //compute minimum over all widths (shared min_width?!)
    }

    return min_width; //also return the elimination ordering causing minimal width?
}

} //namespace treedec

#endif //TD_RANDOM
