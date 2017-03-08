#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY
#define GENERIC_ELIMINATION_SEARCH_OVERLAY

#include <boost/graph/adjacency_list.hpp>
#include <stack>

namespace treedec{

namespace gen_search{

template <typename UnderlyingG_t, typename OverlayG_t> //UnderlyingG_t should be gala_vec_sorted, Overlay should be gala_vec_unsorted
class overlay{
public:
    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vdU;

    overlay(UnderlyingG_t &G_input)
      : G(G_input)
    {
        _active = std::vector<bool>(boost::num_vertices(G_input), true);
        for(unsigned i = 0; i < boost::num_vertices(G_input); ++i)
        {
            boost::add_vertex(O);
        }
    }

    overlay(UnderlyingG_t &G_input, std::vector<bool> &active_input) //e.g. after PP
      : G(G_input), _active(active_input)
    {
        for(unsigned i = 0; i < boost::num_vertices(G_input); ++i)
        {
            boost::add_vertex(O);
        }
    }

    const UnderlyingG_t &underlying(){
        return G;
    }

    const std::vector<bool> &active(){
        return _active;
    }


    /* TODO:
        -actual degree as in DEGREE..
        -Underlying should be const, vec and sorted
        -N(elim_vertex) = N_U(elim_vertex) + N_O(elim_vertex)
        -make clique:
          -sort N(elim vertex)
          - (binsearch in Underlying, linear in Overlay)
          -> NOT deg-many binsearch on the whole outedgevec! (because N(elim_v) is sorted, the search range reduces!)
        -add edges just in Overlay
        -changes_container should be a stack of pair<uint, vec<uint> > with |vec<uint>| = actual_degree
          -> meaning of pair<uint, vec<uint> >: first: modified vertex in overlay, second: #addition edges
          -> undo is stack.back(), then resize overlay[pair.first] according to vec<uint>[i], then stack.pop()
    */
    unsigned eliminate(vdU elim_vertex)
    {
        _active[elim_vertex]= false;

        _changes_container.push(std::vector<vdU>());

        unsigned actual_degree = 0;

        typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator nIt1, nIt2, nEnd1, nEnd2;
        boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(elim_vertex, G);
        boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(elim_vertex, O);

        for(; nIt1 != nEnd2; ++nIt1){ //this just goes through N_G(elim) and than through N_O(elim), TODO: iterator
            if(nIt1 == nEnd1){
                boost::tie(nIt1, nEnd2) = boost::adjacent_vertices(elim_vertex, O);
                if(nIt1 == nEnd2){
                    break;
                }
            }

            if(!_active[*nIt1]){
                continue;
            }

            ++actual_degree;

            nIt2 = nIt1;
            ++nIt2;

            for(; nIt2 != nEnd2; ++nIt2){
                if(nIt2 == nEnd1){
                    boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(elim_vertex, O);
                    if(nIt2 == nEnd2){
                        break;
                    }
                }

                if(!_active[*nIt2]){
                    continue;
                }

                //TODO: can be further improved..
                if(!boost::edge(*nIt1, *nIt2, G).second && !boost::edge(*nIt1, *nIt2, O).second)
                {
                    boost::add_edge(*nIt1, *nIt2, O);
                    boost::add_edge(*nIt2, *nIt1, O);
                    _changes_container.top().push_back(*nIt1);
                    _changes_container.top().push_back(*nIt2);
                }
            }
        }
        return actual_degree;
    }

    void undo_eliminate(vdU elim_vertex)
    {
        _active[elim_vertex]= true;
        while(!_changes_container.top().empty()){
            typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor v1 = _changes_container.top().back();
            _changes_container.top().pop_back();
            typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor v2 = _changes_container.top().back();
            _changes_container.top().pop_back();

            boost::remove_edge(v1, v2, O);
            boost::remove_edge(v2, v1, O);
        }
        _changes_container.pop();
    }

private:
    const UnderlyingG_t &G;
    OverlayG_t O;
    std::vector<bool> &_active;

    std::stack<std::vector<vdU> > _changes_container;
};

} //namespace gen_search

} //namespace treedec

#endif //guard
