#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY
#define GENERIC_ELIMINATION_SEARCH_OVERLAY

#include <boost/graph/adjacency_list.hpp>
#include <stack>

namespace treedec{

namespace gen_search{

template <typename UnderlyingG_t, typename OverlayG_t>
class overlay{
public:
    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vdU;

    overlay(UnderlyingG_t const &G_input)
      : G(G_input)
    {
        _active = std::vector<BOOL>(boost::num_vertices(G_input), true);
        O = OverlayG_t(boost::num_vertices(G_input));
    }

    overlay(UnderlyingG_t &G_input, std::vector<BOOL> &active_input) //e.g. after PP
      : G(G_input), _active(active_input)
    {
        O = OverlayG_t(boost::num_vertices(G_input));
    }

    const UnderlyingG_t &underlying(){
        return G;
    }

    const std::vector<BOOL> &active(){
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

private:
public: /// bug. accessed from outside.
    const UnderlyingG_t &G;
    OverlayG_t O;
public: // BUG. wrong class
    std::vector<BOOL> &_active;

public: /// bug. accessed from outside.
	 // BUG: inefficient.
    std::stack<std::vector<vdU> > _changes_container;
};


#if 0
template <typename G_t, typename VD_t>
void gala_resize(G_t &G, VD_t v, unsigned num){
    auto &g = G.vertices();
    g[v].resize(g[v].size()-num);
}

template <typename UnderlyingG_t, typename OverlayG_t> //UnderlyingG_t should be gala_vec_sorted, Overlay should be gala_vec_unsorted
class overlay_gala : public overlay<UnderlyingG_t, OverlayG_t>{
private:
	template<class i1, class i2>
	using concat_iterator=draft::concat_iterator<i1, i2>;
public:
    typedef overlay<UnderlyingG_t, OverlayG_t> baseclass;

    typedef typename baseclass::vdU vdU;

    overlay_gala(UnderlyingG_t &G_input)
      : overlay<UnderlyingG_t, OverlayG_t>(G_input)
    {}

    overlay_gala(UnderlyingG_t &G_input, std::vector<BOOL> &active_input) //e.g. after PP
      : overlay<UnderlyingG_t, OverlayG_t>(G_input, active_input)
    {}

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
        baseclass::_active[elim_vertex]= false;

        baseclass::_changes_container.push(std::vector<vdU>());
        _changes_size.push(std::vector<unsigned>(baseclass::_active.size()));

        unsigned actual_degree = 0;

        typedef typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator adj1_iterator;
        typedef typename boost::graph_traits<OverlayG_t>::adjacency_iterator adj2_iterator;
        adj1_iterator nIt1, nEnd1;
        adj2_iterator nIt2, nEnd2;

        boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(elim_vertex, baseclass::G);
        boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(elim_vertex, baseclass::O);

        concat_iterator<adj1_iterator, adj2_iterator> cIt1(nIt1, nEnd1, nIt2, nEnd2);
        concat_iterator<adj1_iterator, adj2_iterator> cIt2(nIt1, nEnd1, nIt2, nEnd2);

        for(; cIt1 != nEnd2; ++cIt1){
            if(!baseclass::_active[*cIt1]){
                continue;
            }

            ++actual_degree;

            baseclass::_changes_container.top().push_back(*cIt1);

            cIt2 = cIt1;
            ++cIt2;

            for(; cIt2 != nEnd2; ++cIt2){
                if(!baseclass::_active[*cIt2]){
                    continue;
                }

                //TODO: can be further improved..
                //if cIt1 or cIt2 are not in G, than the first one (! bla) is always true
                if(!boost::edge(*cIt1, *cIt2, baseclass::G).second && !boost::edge(*cIt1, *cIt2, baseclass::O).second)
                {
                    boost::add_edge(*cIt1, *cIt2, baseclass::O);
                    boost::add_edge(*cIt2, *cIt1, baseclass::O);

                    ++_changes_size.top()[*cIt1];
                    ++_changes_size.top()[*cIt2];
                }
            }
        }
        return actual_degree;
    }

    void undo_eliminate(vdU elim_vertex)
    {
        baseclass::_active[elim_vertex]= true;
        for(unsigned i = 0; i < baseclass::_changes_container.top.size(); ++i){
            vdU v = baseclass::_changes_container.top()[i];
            gala_resize(baseclass::O, v, _changes_size.top()[v]);
        }
        baseclass::_changes_container.pop();
        _changes_size.pop();
    }

private:
    std::stack<std::vector<unsigned> > _changes_size;
};
#endif

} //namespace gen_search

} //namespace treedec

#endif //guard
