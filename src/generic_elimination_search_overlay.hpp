#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY
#define GENERIC_ELIMINATION_SEARCH_OVERLAY

#include <boost/graph/adjacency_list.hpp>

namespace treedec{

namespace gen_search{

template <typename UnderlyingG_t, typename OverlayG_t>
class overlay{
public:
    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vdU;

    overlay(UnderlyingG_t &G_input)
      : G(G_input)
    {
        _active = std::vector<bool>(boost::num_vertices(G_input), true);
    }

    overlay(UnderlyingG_t &G_input, std::vector<bool> &active_input) //e.g. after PP
      : G(G_input), _active(active_input)
    {
    }

    /*const*/ UnderlyingG_t &underlying(){
        return G;
    }

    std::vector<bool> &active(){
        return _active;
    }

    unsigned eliminate(vdU elim_vertex)
    {
        _active[elim_vertex]= false;

        unsigned actual_degree = 0;

        typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator nIt1, nIt2, nEnd;
        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(elim_vertex, G); nIt1 != nEnd; ++nIt1){
            if(!_active[*nIt1]){
                continue;
            }

            ++actual_degree;

            nIt2 = nIt1;
            ++nIt2;
            for(; nIt2 != nEnd; ++nIt2){
                if(!_active[*nIt2]){
                    continue;
                }
                if(!boost::edge(*nIt1, *nIt2, G).second){
                    boost::add_edge(*nIt1, *nIt2, G);
                    boost::add_edge(*nIt2, *nIt1, G);
                    _changes_container.push_back(*nIt1);
                    _changes_container.push_back(*nIt2);
                }
            }
        }
        return actual_degree;
    }

    void undo_eliminate(vdU elim_vertex)
    {
        _active[elim_vertex]= true;
        unsigned c = _changes_container.size() >> 1;
        for(unsigned i = 0; i < c; ++i){
            typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor v1 = _changes_container.back();
            _changes_container.pop_back();
            typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor v2 = _changes_container.back();
            _changes_container.pop_back();

            boost::remove_edge(v1, v2, G);
            boost::remove_edge(v2, v1, G);
        }
    }

private:
    /*const*/ UnderlyingG_t &G;
    OverlayG_t O;
    std::vector<bool> _active;

    std::vector<vdU> _changes_container;
};

} //namespace gen_search

} //namespace treedec

#endif //guard
