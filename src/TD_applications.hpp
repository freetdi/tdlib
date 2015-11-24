// Lukas Larisch, 2014 - 2015
//
// (c) 2014-2015 Goethe-Universität Frankfurt
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//
// Offers functionality to solve hard problems with help of treedecompositions
//

#define HAVE_CLIQUER

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <vector>
#include <boost/graph/adjacency_list.hpp>

#ifdef HAVE_CLIQUER
#include "cliquer.h"
#endif


namespace treedec{

namespace app{

template <typename G_t>
class max_clique_base{
    public:
        virtual void max_clique(G_t &G, std::vector<unsigned int> &result) = 0;
};

#ifdef HAVE_CLIQUER

template <typename G_t>
class max_clique_cliquer : public max_clique_base<G_t>{
    public:
        void max_clique(G_t &G, std::vector<unsigned int> &result){
            graph_t *h;
            set_t s;
            clique_options *opts;

            std::vector<unsigned int> id_map;
            treedec::reorder_ids_graph(G, id_map);

            h = graph_new(boost::num_vertices(G));

            typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
            for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                    if(G[*nIt].id > G[*vIt].id)
                        GRAPH_ADD_EDGE(h, G[*vIt].id, G[*nIt].id);
                }
            }

        /* Initialize out clique_options */
        opts=(clique_options*)malloc(sizeof(clique_options));
        opts->time_function=NULL;
        opts->output=stderr;
        opts->reorder_function=NULL;
        opts->reorder_map=NULL;
        opts->user_function=NULL;
        opts->user_data=NULL;
        opts->clique_list=NULL;
        opts->clique_list_length=0;

            /* max-clique call */
        s=clique_unweighted_find_single(h, 0, 0, TRUE, opts);

/*
            std::cout << "induced subgraph: " << std::endl;
            for(unsigned int i = 0; i < boost::num_vertices(G); i++)
                std::cout << id_map[G[i].id] << " ";

            std::cout << std::endl;
*/

            int size=set_size(s);
            for(unsigned int i = 0; i < SET_MAX_SIZE(s); i++){
                if(SET_CONTAINS(s,i))
                    result.push_back(id_map[i]);
            }

/*
            std::cout << "maxclique: " << std::endl;
            for(unsigned int i = 0; i < result.size(); i++)
                std::cout << result[i] << " ";
            std::cout << std::endl;
*/

            /* free all stuff */
            free(opts);
            graph_free(h);
        }

};

#endif


template <typename G_t, typename T_t>
void max_clique_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, max_clique_base<G_t> &mcb){
    for(unsigned int i = 0; i < boost::num_vertices(T); i++){
        G_t H;
        induced_subgraph(H, G, T[i].bag);
        std::vector<unsigned int> result;
        mcb.max_clique(H, result);
        if(result.size() > global_result.size())
            global_result = result;
    }
}

template <typename G_t, typename T_t>
void max_vertex_cover(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void max_independent_set(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void max_coloring(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void hamiltonian_cycle(G_t &G, T_t &T){
}

}

}

#endif


