#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY_H
#define GENERIC_ELIMINATION_SEARCH_OVERLAY_H

#ifdef USE_GALA
#include <gala/boost.h>
#endif

#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <stack>

#include "iter.hpp"
#include "marker.hpp"
#include "trace.hpp"

namespace treedec{

namespace gen_search{

namespace detail{

struct pcnt{
    pcnt(size_t x)
        : _cnt(1+x)
    {
    }
    template<class T>
    bool operator()(T){
        if(_cnt){
            --_cnt;
        }
        return !_cnt;
    }
    size_t _cnt;
};

//TODO
template<class G>
struct outedge_resize{
    static void do_it(typename boost::graph_traits<G>::vertex_descriptor,
                     size_t, G&)
    {
        //static_assert(sizeof(G)==0);
    }
};

// TODO:: more generic.
#ifdef USE_GALA
template<class G>
struct dvv_config : public gala::graph_cfg_default<G> {
	static constexpr bool is_directed=true;
};
typedef gala::graph<std::vector, std::vector, uint32_t, dvv_config> gdvv;

    BOOST_STATIC_ASSERT(
            !std::is_convertible<typename boost::graph_traits<gdvv>::traversal_category*,
                         boost::bidirectional_graph_tag* >::value);

template<>
struct outedge_resize<gdvv>{
    static void do_it(typename boost::graph_traits<gdvv>::vertex_descriptor x,
                     size_t size, gdvv& g)
    {
        g.out_edges(x).resize(size);
    }
};
#endif

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> baluo;
template<>
struct outedge_resize<baluo>{
static void do_it(typename boost::graph_traits<baluo>::vertex_descriptor x,
		 size_t size, baluo& g)
{
    assert(boost::out_degree(x, g)>=size);
    auto P=pcnt(size);
    boost::remove_out_edge_if(x, P, g);
    assert( size == boost::out_degree(x, g));
}

};

template <typename G_t, typename VD_t>
void delete_top_edges(G_t &G, VD_t v, unsigned howmany){
    auto deg = boost::out_degree(v, G);
    assert(howmany<=deg);
    assert(howmany>=0);
    outedge_resize<G_t>::do_it(v, deg - howmany, G);

    assert( deg-howmany == boost::out_degree(v, G));
}

} // detail


template <typename UnderlyingG_t, typename OverlayG_tt, class ACTMAP>
class overlay{
public:

#ifdef USE_GALA
    typedef detail::gdvv OverlayG_t;
#else
    typedef detail::baluo OverlayG_t;
#endif

    typedef typename boost::graph_traits<UnderlyingG_t>::adjacency_iterator adj1_iterator;
    typedef typename boost::graph_traits<OverlayG_t>::adjacency_iterator adj2_iterator;
    typedef draft::concat_iterator<adj1_iterator, adj2_iterator> all_adj_it;

    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vdU;
    typedef typename boost::graph_traits<UnderlyingG_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<UnderlyingG_t>::vertices_size_type vertices_size_type;
    typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;

    typedef typename boost::graph_traits<UnderlyingG_t>::directed_category      directed_category;
    typedef typename boost::graph_traits<UnderlyingG_t>::edge_parallel_category edge_parallel_category;
    typedef typename boost::graph_traits<UnderlyingG_t>::traversal_category     traversal_category;

    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;

    template<class i1, class i2>
    using concat_iterator=draft::concat_iterator<i1, i2>;

    struct active_filter{
        active_filter(ACTMAP const & v) : _v(v) {}
        bool operator()(vertex_descriptor v) const{
            return _v[v];
        }
        ACTMAP const& _v;
    };
    typedef boost::filter_iterator<active_filter, all_adj_it> adjacency_iterator;

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
    overlay(UnderlyingG_t const& g, ACTMAP m) // (, std::vector<BOOL> &active_input) //e.g. after PP
      : _g(g),
        _og(boost::num_vertices(g)),
        _active(m),
        _degree(boost::num_vertices(g)),
        _marker(boost::num_vertices(g))
    {
        commit();
	assert(_changes_container.size()==1);
        auto vs=boost::vertices(g);
        for(; vs.first!=vs.second; ++vs.first){
            _degree[*vs.first]=boost::out_degree(*vs.first, _g);
        }
    }

private:
    overlay(const overlay&o)
        : _g(o._g),
          _og(o._og),
          _changes_container(o._changes_container),
          _active(o._active),
          _degree(o._degree),
          _marker(o._marker)
    {
        unreachable();
        assert(_changes_container.size()==1); // for now
    }

public:
    unsigned num_vertices() const{return boost::num_vertices(_g);}
    std::pair<edge_descriptor, bool> edge(vertex_descriptor a, vertex_descriptor b) const{
        assert(_active[a]);
        assert(_active[b]);

        auto e=boost::edge(a, b, _g);
        auto P=std::make_pair(a,b);
        if(e.second){
            return std::make_pair(P, true);
        }else{
            return std::make_pair(P, boost::edge(a, b, _og).second);
        }
    }
    std::pair<edge_descriptor, bool> add_edge(vertex_descriptor a, vertex_descriptor b){
        incomplete();
        assert(_active[a]);
        assert(_active[b]);
        auto e=boost::add_edge(a, b, _og);

        _changes_container.top().push_back(a);
        _changes_container.top().push_back(b);
        return std::make_pair(std::make_pair(a, b), e.second);
    }
    void commit(){
        _changes_container.emplace(0);
    }
    void reset(unsigned ref=0);

private:
    void reset_neigh(vertex_descriptor v);

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
        auto i=all_adj_it(p.first, p.second, q.first, q.second);
        auto e=all_adj_it(p.second, p.second, q.second, q.second);
#endif
        typedef boost::filter_iterator<active_filter, all_adj_it> FilterIter;
        active_filter fP(_active);

        FilterIter fb(fP, i, e);
        FilterIter fe(fP, e, e);

        return std::make_pair(fb, fe);
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
    void eliminate(vertex_descriptor g);
     // can only undo the previous elimination.
    vertex_descriptor undo_eliminate();
    vertices_size_type degree(vertex_descriptor v)const{
#ifndef NDEBUG
        auto p=adjacent_vertices(v);
        unsigned cnt=0;
        for(; p.first!=p.second; ++p.first){
            assert(_active[*p.first]);
            ++cnt;
        }

        assert(!_active[v] || cnt==_degree[v]);
#endif
        return _degree[v];
    }

private:
    ACTMAP& active(){
        return _active;
    }

//private:

public: //accessed from outside.
    const UnderlyingG_t &_g;
    OverlayG_t _og;

private:
    std::stack<std::vector<vdU> > _changes_container;
    std::stack<long> _elim_stack; // need union {descr, degree}
    ACTMAP _active; // active and current_ordering -> numbering.
    std::vector<vertices_size_type> _degree; // active and current_ordering -> numbering.
    marker_type _marker;
}; // overlay

template<class A, class B, class C>
void treedec::gen_search::overlay<A, B, C>::reset(unsigned i)
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

template<class A, class B, class C>
void treedec::gen_search::overlay<A, B, C>::reset_neigh(vertex_descriptor v)
{
#if OLDSTUFF
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
#endif
    auto p=adjacent_vertices(v);
    std::vector<long> reverse;

    for(; p.first!=p.second; ++p.first){
        assert(!_elim_stack.empty());
        auto howmany=_elim_stack.top();
        _elim_stack.pop();
        reverse.push_back(howmany);
    }
    assert(_degree[v]==reverse.size());

    auto reverseit=reverse.rbegin();
    auto p2=adjacent_vertices(v); // use p2...
    for(; p2.first!=p2.second; ++p2.first){
        assert(*reverseit+1>=0);
#ifndef NDEBUG
        auto olddegree=boost::out_degree(*p2.first, _og);
#endif
        detail::delete_top_edges(_og, *p2.first, *reverseit + 1);
        assert(boost::out_degree(*p2.first, _og) == olddegree - *reverseit - 1);
        _degree[*p2.first] -= *reverseit;
        assert(_degree[*p2.first]<boost::num_vertices(_g));
        ++reverseit;
    }
} // reset;

} // gen_search

} // treedec

namespace boost {

template<class A, class B, class C>
std::pair<typename treedec::gen_search::overlay<A, B, C>::adjacency_iterator,
          typename treedec::gen_search::overlay<A, B, C>::adjacency_iterator>
  adjacent_vertices(
          unsigned v,
          treedec::gen_search::overlay<A, B, C> const& o)
{
	return o.adjacent_vertices(v);
}

template<class A, class B, class C>
std::pair<typename treedec::gen_search::overlay<A, B, C>::edge_descriptor, bool>
  edge(unsigned a, unsigned b, treedec::gen_search::overlay<A, B, C> const& o)
{
	return o.edge(a, b);
}

template<class A, class B, class C>
std::pair<typename treedec::gen_search::overlay<A, B, C>::edge_descriptor, bool>
  add_edge(unsigned a, unsigned b, treedec::gen_search::overlay<A, B, C>& o)
{
    return o.add_edge(a, b);
}

template<class A, class B, class C>
unsigned num_vertices(treedec::gen_search::overlay<A, B, C> const& o)
{
	return o.num_vertices();
}


} // boost

namespace treedec {

namespace gen_search {

#if 0

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

    vertex_descriptor undo_eliminate()
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

template <class A, class B, class C>
void overlay<A, B, C>::eliminate(
        typename overlay<A, B, C>::vertex_descriptor elim_vertex)
{
    auto nv=boost::num_vertices(_g); (void)nv;
    auto olddegree=_degree; // hmm.
    assert(_active[elim_vertex]);
    active()[elim_vertex] = false;
    (void)degree(elim_vertex);

    unsigned actual_degree = 0; // recompute degree of center vertex?

    auto p=adjacent_vertices(elim_vertex);
    for(; p.first!=p.second; ++p.first){
        assert(actual_degree<_degree[elim_vertex]);
        assert(*p.first<nv);
        assert(*p.first != elim_vertex);
        assert(active()[*p.first]);

        _marker.clear();
        // mark_smaller_neighbours(_marker, *p.first, *this); doesntwork
        auto cnt=mark_smaller_neighbours(_marker, *p.first, *this, active());
        assert(cnt<=_degree[*p.first]); (void)cnt;

        --_degree[*p.first];
        assert(_degree[*p.first]<nv);

        auto q=adjacent_vertices(elim_vertex);
        for(; q.first!=q.second; ++q.first){
            assert(active()[*q.first]);
            if(*q.first>=*p.first){
                // skip. TODO: more efficient skip
            }else if(_marker.is_marked(*q.first)){
                // done.
            }else{
                ++_degree[*p.first];
                ++_degree[*q.first];
                assert(_degree[*p.first]<nv);
                assert(_degree[*q.first]<nv);
                assert(!boost::edge(*p.first, *q.first, _og).second);
                trace2("overlay add", *p.first, *q.first);
                treedec::add_edge(*p.first, *q.first, _og);

                // treedec graph iface enforces this. (.. should)
                assert(boost::edge(*p.first, *q.first, _og).second);
                assert(boost::edge(*q.first, *p.first, _og).second);
            }
        }
        assert(_degree[elim_vertex]<nv);
        ++actual_degree;
    }

    assert(_degree[elim_vertex]==actual_degree);


    auto p2=adjacent_vertices(elim_vertex);
    for(; p2.first!=p2.second; ++p2.first){
        // YUCK. mixing vertex_descriptor and vertices_size_t.
        long delta = _degree[*p2.first] - olddegree[*p2.first];
        _elim_stack.push(delta);
    }

    _elim_stack.push(elim_vertex);

#ifndef NDEBUG
    auto p3=boost::vertices(_g);
    for(; p3.first!=p3.second; ++p3.first){
    }
#endif

//    return actual_degree;
}

// YUCK. this whole eliminationj stuff should happen in generic_base
template <class A, class B, class C>
typename overlay<A, B, C>::vertex_descriptor overlay<A, B, C>::undo_eliminate()
{
    assert(!_elim_stack.empty());
    auto elim_vertex=_elim_stack.top();
    assert(!active()[elim_vertex]);
    active()[elim_vertex] = true;
    _elim_stack.pop();
    reset_neigh(elim_vertex);

#ifndef NDEBUG
    auto p=boost::vertices(_g);
    for(; p.first!=p.second; ++p.first){
        degree(*p.first);
    }
#endif

    return elim_vertex;
}


} //namespace gen_search

} //namespace treedec

#endif //guard

// vim:ts=8:sw=4:et
