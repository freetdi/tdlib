#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY
#define GENERIC_ELIMINATION_SEARCH_OVERLAY

#include <boost/graph/adjacency_list.hpp>
#include <stack>

namespace treedec{

namespace gen_search{


template<class iter1, class iter2> //TODO: fix this
class concat_iterator{
public:
    concat_iterator(iter1 begin1, iter1 end1, iter2 begin2, iter2 end2)
     : _i1(begin1), _e1(end1), _i2(begin2), _e2(end2){}

    bool operator!=(const iter2& end){
        // warning: only works if end==end_of range2.

        return !(_i1 ==_e1 && _i2 == _e2);
    }

    void operator++(){
        if(_i1!=_e1){
            // still busy with range 1
            ++_i1;
        }
        else{
            ++_i2;
        }
    }

    unsigned operator*(){
        if(_i1!=_e1){
            //still busy with range 1
            return *_i1;
        }

        return *_i2;
    }

    bool is_in_underlying(){
        return _i1!=_e1;
    }

private:
    iter1 _i1, _e1;
    iter2 _i2, _e2;
};


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

        typedef typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator adj1_iterator;
        typedef typename boost::graph_traits<OverlayG_t>::adjacency_iterator adj2_iterator;
        adj1_iterator nIt1, nEnd1;
        adj2_iterator nIt2, nEnd2;
        boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(elim_vertex, G);
        boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(elim_vertex, O);

        concat_iterator<adj1_iterator, adj2_iterator> cIt1(nIt1, nEnd1, nIt2, nEnd2);
        concat_iterator<adj1_iterator, adj2_iterator> cIt2(nIt1, nEnd1, nIt2, nEnd2);

        for(; cIt1 != nEnd2; ++cIt1){
            if(!_active[*cIt1]){
                continue;
            }

            ++actual_degree;

            cIt2 = cIt1;
            ++cIt2;

            for(; cIt2 != nEnd2; ++cIt2){
                if(!_active[*cIt2]){
                    continue;
                }

                //TODO: can be further improved..
                //if cIt1 or cIt2 are not in G, than the first one (! bla) is always true
                if(cIt1.is_in_underlying() && cIt2.is_in_underlying() && boost::edge(*cIt1, *cIt2, G).second){
                    continue;
                }
                if(!boost::edge(*cIt1, *cIt2, O).second)
                {
                    boost::add_edge(*cIt1, *cIt2, O);
                    boost::add_edge(*cIt2, *cIt1, O);
                    _changes_container.top().push_back(*cIt1);
                    _changes_container.top().push_back(*cIt2);
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
