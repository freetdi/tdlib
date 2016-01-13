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
// Implements some exact algorithms for hard problems. The algorithms are
// used in TD_applications.
//


#ifndef TD_HARD_PROBLEMS
#define TD_HARD_PROBLEMS

#ifdef HAVE_CLIQUER
#define new new_foo
//you have to use -fpermissive in compilation 
extern "C"{
#include <cliquer/cliquer.h>
#undef new
}
#endif

#include <tdlib/TD_misc.hpp>

namespace treedec{

namespace np{

template <typename G_t>
class max_clique_base{
    public:
        virtual void max_clique(G_t &G, std::vector<unsigned int> &result) = 0;
};

template <typename G_t>
class max_independent_set_base{
    public:
        virtual void max_independent_set(G_t &G, std::vector<unsigned int> &result) = 0;
};

template <typename G_t>
class min_vertex_cover_base{
    public:
        virtual void min_vertex_cover(G_t &G, std::vector<unsigned int> &result) = 0;
};

template <typename G_t>
class min_coloring_base{
    public:
        virtual void min_coloring(G_t &G, std::vector<std::vector<unsigned int> > &result) = 0;
};

template <typename G_t>
class min_dominating_set_base{
    public:
        virtual void min_dominating_set(G_t &G, std::vector<unsigned int> &result) = 0;
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

            /* Initialize clique_options */
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

            int size=set_size(s);
            for(unsigned int i = 0; i < SET_MAX_SIZE(s); i++){
                if(SET_CONTAINS(s,i))
                    result.push_back(id_map[i]);
            }

            /* free all stuff */
            free(opts);
            graph_free(h);
        }

};

template <typename G_t>
class max_independent_set_cliquer : public max_independent_set_base<G_t>{
    public:
        void max_independent_set(G_t &G, std::vector<unsigned int> &result){
            G_t cG;
            complement_graph(cG, G);

            treedec::np::max_clique_cliquer<G_t> A;
            A.max_clique(cG, result);
        }

};

#endif //HAVE_CLIQUER

} //namespace np

} //namespace treedec


#endif
