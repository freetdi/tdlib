
#ifndef TD_RANDOM

#include "TD_elimination_orderings.hpp"

namespace treedec{

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename G_t>
int randomly_try_some_elimination_orderings(G_t &G,
       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &best_elim_ordering,
       unsigned int count = 5)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> elim_ordering(boost::num_vertices(G));
    for(unsigned int i=0; i<elim_ordering.size(); ++i){
        elim_ordering[i] = i;
    }

    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > elimination_orderings(count);
    for(unsigned int i=0; i<elimination_orderings.size(); ++i){
         //Use the built-in random generator.
        std::random_shuffle(elim_ordering.begin(), elim_ordering.end());
        elimination_orderings[i] = elim_ordering;
    }

    int min_width = INT_MAX;

    #pragma omp parallel for
    for(unsigned int i = 0; i < count; i++){
        G_t H;
        boost::copy_graph(G, H); // ..(H, G)..?! "unavoidable"?
        int width_i = treedec::get_width_of_elimination_ordering(H, elimination_orderings[i]);

        if(width_i < min_width){
            best_elim_ordering = MOVE(elimination_orderings[i]);
        }
        //std::cout << "width_" << i << ": " << width_i << std::endl;
        //compute minimum over all widths (shared min_width?!)
    }

    return min_width; //also return the elimination ordering causing minimal width?
}

} //namespace treedec

#endif //TD_RANDOM
