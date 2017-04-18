#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY_H
#define GENERIC_ELIMINATION_SEARCH_OVERLAY_H

#include <boost/graph/adjacency_list.hpp>
#include <stack>

#include "iter.hpp"

namespace treedec{

namespace gen_search{

// BUG: inefficient.
template <typename UnderlyingG_t, typename OverlayG_t>
class overlay{
public:
    typedef typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator adj1_iterator;
    typedef typename boost::graph_traits<OverlayG_t>::adjacency_iterator adj2_iterator;

    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vdU;
    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vertex_descriptor;
    typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;

    typedef typename OverlayG_t::directed_category      directed_category;
    typedef typename OverlayG_t::edge_parallel_category edge_parallel_category;
    typedef typename OverlayG_t::traversal_category     traversal_category;

    template<class i1, class i2>
    using concat_iterator=draft::concat_iterator<i1, i2>;

    typedef draft::concat_iterator<adj1_iterator, adj2_iterator> adjacency_iterator;

public: // construct
#if 0 // not yet.
    overlay(UnderlyingG_t const &G_input)
      : _g(G_input),
        _og(boost::num_vertices(G_input)),
        _active(active) // BUG
    {
        _active = std::vector<BOOL>(boost::num_vertices(G_input), true);
        commit();
	assert(_changes_container.size()==1);
    }
#endif

    overlay(UnderlyingG_t const& g) // (, std::vector<BOOL> &active_input) //e.g. after PP
      : _g(g),
        _og(boost::num_vertices(g))
    {
        commit();
	assert(_changes_container.size()==1);
    }

    // private: // BUG
    overlay(const overlay&o)
        : _g(o._g),
          _og(o._og),
          _changes_container(o._changes_container)
    { untested();
            assert(_changes_container.size()==1);
    }

public:
    unsigned num_vertices() const{return boost::num_vertices(_g);}
    std::pair<edge_descriptor, bool> edge(vertex_descriptor a, vertex_descriptor b) const{
        auto e=boost::edge(a, b, _g);
        auto P=std::make_pair(a,b);
        if(e.second){ untested();
            return std::make_pair(P, true);
        }else{
            return std::make_pair(P, boost::edge(a, b, _og).second);
        }
    }
    std::pair<edge_descriptor, bool> add_edge(vertex_descriptor a, vertex_descriptor b){
        auto e=boost::add_edge(a, b, _og);

        _changes_container.top().push_back(a);
        _changes_container.top().push_back(b);
        return std::make_pair(std::make_pair(a, b), e.second);
    }
    void commit(){
        _changes_container.emplace(0);
    }
    void reset(unsigned ref=0);
public:
	std::pair<adjacency_iterator, adjacency_iterator>
	adjacent_vertices(vertex_descriptor v) const{
            using draft::concat_iterator;
		auto p=boost::adjacent_vertices(v, _g);
		auto q=boost::adjacent_vertices(v, _og);
#if 0 //not yet
		auto j=boost::range::join(p,q);
		auto i=j.begin();
		auto e=j.end();
#else
		typedef typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator adj1_iterator;
		typedef typename boost::graph_traits<OverlayG_t>::adjacency_iterator adj2_iterator;
		auto i=concat_iterator<adj1_iterator, adj2_iterator>(p.first, p.second, q.first, q.second);
      auto e=concat_iterator<adj1_iterator, adj2_iterator>(p.second, p.second, q.second, q.second);
#endif
		return std::make_pair(i, e);
	}

    const UnderlyingG_t &underlying() const{
        return _g;
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
    const UnderlyingG_t &_g;
    OverlayG_t _og;
private:
    std::stack<std::vector<vdU> > _changes_container;
}; // overlay

template<class A, class B>
void treedec::gen_search::overlay<A, B>::reset(unsigned i)
{
    assert(i==1); // for now.
    assert(_changes_container.size());
    assert(_changes_container.top().empty()); // for now
    _changes_container.pop();
    assert(_changes_container.size());

    while(!_changes_container.top().empty()){
        auto v1=_changes_container.top().back();
        _changes_container.top().pop_back();
        auto v2=_changes_container.top().back();
        _changes_container.top().pop_back();

        assert(boost::edge(v1, v2, _og).second);
        assert(boost::edge(v2, v1, _og).second);
        boost::remove_edge(v1, v2, _og);
        boost::remove_edge(v2, v1, _og);
    }
}

} // gen_search

} // treedec

namespace boost {

template<class A, class B>
std::pair<typename treedec::gen_search::overlay<A, B>::adjacency_iterator,
          typename treedec::gen_search::overlay<A, B>::adjacency_iterator>
adjacent_vertices(
          //typename treedec::gen_search::overlay<A, B>::vertex_descriptor v,
          unsigned v,
          treedec::gen_search::overlay<A, B> const& o)
{
	return o.adjacent_vertices(v);
}

template<class A, class B>
std::pair<typename treedec::gen_search::overlay<A, B>::edge_descriptor, bool>
edge(unsigned a, unsigned b, treedec::gen_search::overlay<A, B> const& o)
{
	return o.edge(a, b);
}

template<class A, class B>
std::pair<typename treedec::gen_search::overlay<A, B>::edge_descriptor, bool>
add_edge(unsigned a, unsigned b, treedec::gen_search::overlay<A, B>& o)
{
    return o.add_edge(a, b);
}

template<class A, class B>
unsigned num_vertices(treedec::gen_search::overlay<A, B> const& o)
{
	return o.num_vertices();
}


} // boost

namespace treedec {

namespace gen_search {

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

// vim:ts=8:sw=4:et
