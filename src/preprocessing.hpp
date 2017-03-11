// Lukas Larisch, 2014 - 2016
// Felix Salfelder 2016 - 2017
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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

// #define NEGATIVE_TAGS1
// #define NEGATIVE_TAGS2

/*
 * Offers functionality to preprocess a graph, such that after
 * the repeated application of reduction rules, which in case the
 * input graph has tree-width at most 3 allow us to determine it's tree-width exactly
 * and in addition compute the corresponding tree decomposition. If the tree-width
 * is larger, the reduction rules return a possibly smaller instance of the same
 * tree-width as the original graph, a 'partial' tree decomposition and a lower bound
 * with respect to tree-width, such that
 * further algorithms can be applied to the resulting graph.
 *
 *
 * These functions are most likely to be interesting for outside use:
 *
 * -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags)
 * -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &lb)
 *
 */

// TODO:
//  - switch to bs
//  - skeletal?
//  - dormant on?
//  - mass elimination

#ifndef TD_PREPROCESSING
#define TD_PREPROCESSING

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/copy.hpp>

#include "directed_view.hpp"
#include "trace.hpp"
#include "algo.hpp"
#include "degree.hpp"
#include "marker.hpp"
#include "numbering.hpp"
#include "copy.hpp"
#include "config_traits.hpp"
#include "graph.hpp"

namespace treedec{

namespace impl{

namespace draft{

template<class G_t>
struct pp_cfg{
    typedef typename
    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags_type;
    //typedef typename deg_chooser<G_t>::type degs_type;
    typedef typename misc::DEGS<G_t> degs_type;
};

} // draft

// default directed view
// there are better ways
template<class G>
struct directed_view_select{
    typedef treedec::draft::directed_view<G> type;
};

namespace detail{

template<class G_t>
struct PP_degree_config : treedec::degs::default_config<G_t> {
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::type idmap_type;
    typedef boost::iterator_property_map<vertex_descriptor*,
        idmap_type, vertex_descriptor, vertex_descriptor&> degree_type;
    /// static_constexpr use_external_degree?
};

}

template<class V, class N, class G, class M>
V deg_vector_init(V const&, N n, G const& g, M const& m)
{
    V v(n);

    typename boost::graph_traits<G>::vertex_iterator I, vend;

    unsigned i=0;
    for (boost::tie(I, vend)=boost::vertices(g); I!=vend; ++I) {
        assert(m[*I]==i);
        assert(i<v.size());
        trace2("i", i, boost::degree(*I, g));

        v[i++] = boost::degree(*I, g);
    }
    return v;
}

//Recursively apply preprocessing rules and glue corresponding bags with
//current tree decomposition this version stores the resulting bags in a vector
//and does not call further algorithms.
//
//TODO: bags is actually a tree decomposition tree.
// FIXME: use bagsize.
template<class G_t, template<class G_> class CFGT=draft::pp_cfg>
class preprocessing : public treedec::algo::draft::algo1 {
private: // hmm, fetch from CFG.
    constexpr static bool disable_triangle=false;
    constexpr static bool disable_cube=false;
    constexpr static bool disable_buddy=false;
    constexpr static bool disable_simplicial=false;
    constexpr static bool disable_almost_simplicial=false;
    constexpr static bool enable_dormant_nodes=false;
    constexpr static bool enable_mass_elimination=false;
    constexpr static bool enable_edge_cleanup=false;
private:
    // STATUS_COUNTER<"somelabel"> sth;
    // STATUS_RATE("someother"); ... ?
public:
    typedef CFGT<G_t> CFG;
    typedef typename directed_view_select<G_t>::type D_t;
    typedef typename boost::graph_traits<D_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<D_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<D_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<D_t>::edges_size_type edges_size_type;
    typedef typename directed_view_select<G_t>::type directed_view_type;
    typedef typename boost::property_map<D_t, boost::vertex_index_t>::type idmap_type;
    //typedef typename deg_chooser<G_t>::type degs_type;
//    typedef typename treedec::config::get::degree<CFG, std::vector<vertices_size_type> >::type
//        degree_type;
    typedef typename std::vector<vertices_size_type> degree_type;

    typedef boost::iterator_property_map<vertex_descriptor*,
        idmap_type, vertex_descriptor, vertex_descriptor&> degreemap_type;
//    typedef typename CFG<directed_view_type>::degs_type degs_type;
    typedef DEGS<D_t, detail::PP_degree_config> degs_type;
    typedef typename CFG::bags_type BV_t;
    typedef typename boost::graph_traits<D_t>::adjacency_iterator adjacency_iterator;
    typedef typename treedec_chooser<G_t>::type T_t; // FIXME
    typedef typename treedec_traits<T_t>::vd_type vd_type;

    typedef treedec::draft::NUMBERING_1<D_t> numbering_type;
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
private:
    template<class A, class N>
    class adjacency_iterator_filter_ : public adjacency_iterator{
    public:
        adjacency_iterator_filter_(const adjacency_iterator_filter_& o) :
            adjacency_iterator(o), _numbering(o._numbering), _end(o._end)
        { }
    public:
        adjacency_iterator_filter_(A a, N const& n, A e)
            : adjacency_iterator(a), _numbering(n), _end(e){
            skip();
        }
        adjacency_iterator_filter_(A a, N const& n)
            : adjacency_iterator(a), _numbering(n), _end(a){
        }
    public:
        adjacency_iterator_filter_& operator=(adjacency_iterator_filter_ const& o){
            adjacency_iterator::operator=(o);
            _end = o._end;
            assert(&_numbering == &o._numbering);
            return *this;
        }
        vertex_descriptor operator*(){ itested();
            assert(*this!=_end);
            return adjacency_iterator::operator*();
        }
        adjacency_iterator_filter_& operator++(){ itested();
            assert(*this!=_end);
            adjacency_iterator::operator++();
            skip();
            return *this;
        }
    private:
        void skip(){ itested();
            while(*this!=_end){ itested();
                if(_numbering.is_numbered(adjacency_iterator::operator*())){ itested();
                    adjacency_iterator::operator++();
                }else{ itested();
                    return;
                }
            }
        }
    private:
        N const& _numbering;
        adjacency_iterator _end;
    };
    typedef adjacency_iterator_filter_<adjacency_iterator, numbering_type>
            adjacency_iterator_filter;
public:
    preprocessing(G_t &G)
        : algo1("pp"), _g(G),
          _id(boost::get(boost::vertex_index, _g)),
          _degree( /*MOVE*/ (deg_vector_init(_degree, boost::num_vertices(_g), _g, _id) )),
          _degreemap(boost::make_iterator_property_map(&_degree[0], _id, _degree[0])),
          _degs(_g, _degreemap),
          _num_edges(boost::num_edges(_g)),
          _marker(boost::num_vertices(_g)),
          _dormant(boost::num_vertices(_g)),
          _numbering(_g)
    {
        assert(_num_edges ^ 1);
        _num_edges /= 2;
        _lb_bs = 1;
    }
private:
    class ISNUM{
    public:
        ISNUM(numbering_type const& n, D_t const& g) : _n(n), _g(g){}
    public:
        template<class E>
        bool operator()(E e) const{
            vertex_descriptor v = boost::target(e, _g);
            return _n.is_numbered(v);
        }
        //bool operator()(vertex_descriptor v) const{ untested();
        //    return _n.is_numbered(v);
        //}
    private:
        numbering_type const& _n;
        D_t const& _g;
    };
    void clear_out_edges(vertex_descriptor v){
        boost::remove_out_edge_if(v, [](edge_descriptor){ return true; }, _g);
    }
public:
    void set_treewidth(std::pair<int, int> p) { untested();
        _lb_bs = p.first+1;
    }
    void set_treewidth(int l, int){
        _lb_bs = l+1;
    }
    void do_it();
    size_t get_bagsize() const{ untested();
        return _lb_bs;
    }
    void get_bags(BV_t& bags) { // for now

        auto b=_elims.begin();
        for(; b!=_elims.end(); ++b){
            bags.emplace_back();
            auto v=*b;
            auto& B=boost::get<1>(bags.back());
            boost::get<0>(bags.back()) = v;

            trace1("BBBBBBBBBBB", v);

            // yes, need boost::
            auto Is=boost::adjacent_vertices(v, _g);
            for(; Is.first!=Is.second; ++Is.first){
                assert(treedec::is_valid(*Is.first, _g));
                if(_numbering.is_before(v, *Is.first)){
                        // num(v) < num(Is.first)
                        // thats why boost numbers reverse?
                    B.insert(*Is.first);
                    trace1("BBBBBBBBBBB", *Is.first);
                }else{ untested();
                }
            }
            // expensive?
            boost::clear_vertex(v, _g);
        }

#if 0 // fixme. do efficiently within loop above
        // remove edges from unnnumbered nodes.
        auto p=boost::vertices(_g);
        ISNUM p_isnum(_numbering, _g);
        for(; p.first!=p.second; ++p.first){
            if(_numbering.is_numbered(*p.first)){
                // don't touch, it's a bag actually
            }else{
                clear_out_edges(*p.first);
            }
        }
#endif
        assert(bags.size()==_elims.size());
    }
    template<class GG>
    void get_graph(GG& gg) {
        assert(boost::is_directed(_g));
        GG gr;
#if 1
        copy_trace(_g, gr);
#else
        boost::copy_graph(_g, gr);
#endif
        assert(boost::is_undirected(gr));
        trace2("", boost::num_edges(gr), boost::num_edges(_g));
        trace2("", boost::num_vertices(gr), boost::num_vertices(_g));
        assert(boost::num_vertices(gr) == boost::num_vertices(_g));
        assert(boost::num_edges(gr)*2 == boost::num_edges(_g));
        gg=gr; //MOVE(gr);
        assert(boost::num_edges(gg)*2 == boost::num_edges(_g));
    }
    int get_treewidth() {
        return int(_lb_bs)-1;
    }
    void isolate(vertex_descriptor v){ untested();
        unsigned deg = _degree[v];
        _num_edges -= deg;

        auto p=adjacent_vertices(v);
        for(; p.first!=p.second; ++p.first){ untested();
            remove_edge(*p.first, v);
        }
        assert(deg == _degree[v]);
    }
    void redegree(vertex_descriptor v, unsigned mark_needs_update=0)
    { // call degree.hpp redegree?
        auto p=adjacent_vertices(v);
        for(; p.first!=p.second; ++p.first){ untested();
            auto n = *p.first;
            assert(_numbering.is_not_numbered(n));

            if(!enable_dormant_nodes){
                assert(!_dormant.is_marked(n));
            }else if(_dormant.is_marked(n)){
                _dormant.unmark(n);
            }else{
            }
            // hmm, which ones do really need update?
            // (later...)
            _degs.update(n);
        }
    }
    bool is_dormant(vertex_descriptor n) const{
        if(!enable_dormant_nodes){
            assert(!_dormant.is_marked(n));
            return false;
        }
        return _dormant.is_marked(n);
    }
    void set_dormant(vertex_descriptor n){
        if(!enable_dormant_nodes){
            return;
        }else{untested();
            assert(!_dormant.is_marked(n));
            _dormant.mark(n);
            _degs.unlink(n);
        }
    }
    void wake_up_node(vertex_descriptor n){
        trace3("wakeup", n, _degree[n], _dormant.is_marked(n));
#ifndef NDEBUG
        trace2("wakeup", _degs.is_reg(0), _degs.is_reg(1));
#endif
        if(_dormant.is_marked(n)){ untested();
            _dormant.unmark(n);
            // wake up!
            _degs.reg(n);
        }else{
            assert(_degs.is_reg(n));
            // back/front?
            //
            _degs.update(n);
            //_degs.reg(n);
        }
        assert(_degs.is_reg(n));
#ifndef NDEBUG
        trace2("wokeup", _degs.is_reg(0), _degs.is_reg(1));
#endif
    }

    void wake_up_neighs(vertex_descriptor v)
    { // call degree.hpp redegree?
        auto p=adjacent_vertices(v);
        for(; p.first!=p.second; ++p.first){
            auto n = *p.first;
            assert(_numbering.is_not_numbered(n));
            wake_up_node(n);
        }
    }
    bool check_twins_3(vertex_descriptor a, vertex_descriptor b) const;
    void isolate_(vertex_descriptor v)
    {
        addtoelims(v);
        // later
        // treedec::make_clique_and_mark(v, _g, _marker);
        //
        _marker.clear();
        // isolate
        auto p=adjacent_vertices(v);
        for(; p.first!=p.second; ++p.first){
            trace2("in", v, *p.first);
            assert(*p.first!=v);
            _marker.mark(*p.first);
            assert(!_numbering.is_numbered(*p.first));
            remove_edge(*p.first, v);
        }
        _num_edges -= _degree[v];
    }
    void make_neigh_clique(vertex_descriptor v, bool isclique=false)
    {
        isolate_(v);
        if(isclique){ untested();
            return;
        }else{
        }

        auto Is=adjacent_vertices(v);
        auto next=Is.first;
        for(; Is.first!=Is.second; Is.first=next){
            ++next;
            auto Ii=next;
            for(; Ii!=Is.second; ++Ii){
                assert(*Is.first != *Ii);
                if(!_marker.is_marked(*Ii)){ untested();
                }else{
                    // hmm how to avoid this?
                    assert(boost::edge(*Is.first, *Ii, _g).second
                            == boost::edge(*Ii, *Is.first, _g).second);
                    add_edge(*Is.first, *Ii);
                    _num_edges += add_edge(*Ii, *Is.first);
                }
            }
        }
    }
private:
    unsigned add_edge(vertex_descriptor v, vertex_descriptor w){
        assert(v!=w);
	if(!boost::edge(v, w, _g).second){
            boost::add_edge(v, w, _g);
            ++_degree[v]; // outdegree
            return 1;
        }else{ itested();
            return 0;
        }
    }
    bool is_numbered(vertex_descriptor v){ untested();
        return _numbering.is_numbered(v);
    }
    void unlink_1_neighbourhood(vertex_descriptor v){ untested();
        auto pp=adjacent_vertices(v);
        for(; pp.first!=pp.second; ++pp.first){ untested();
            _degs.unlink(*pp.first);
        }
    }
    void remove_edge(vertex_descriptor v, vertex_descriptor w){
	assert(boost::edge(v, w, _g).second);
        assert(_numbering.is_numbered(w));
        // boost::remove_edge(v, w, _g);
        --_degree[v];
        assert(_num_edges);
    }
    struct mark_and_remove_helper{
        mark_and_remove_helper(vertex_descriptor a, vertex_descriptor b,
                marker_type& m, D_t const& g, bool eec)
            : _a(a), _b(b), _marker(m), _count(0),
              _g(g), _enable_cleanup(eec){}

        bool operator()(edge_descriptor e){
            auto x=boost::target(e, _g);
            if(x==_a || x==_b){
                return _enable_cleanup;
            }else{
                _marker.mark(x);
                return false;
            }
        }
        vertex_descriptor _a;
        vertex_descriptor _b;
        marker_type& _marker;
        unsigned _count;
        D_t const& _g;
        bool _enable_cleanup;
    };
    void cube_make_clique(vertex_descriptor u, vertex_descriptor v,
            vertex_descriptor w, vertex_descriptor x,
            vertex_descriptor a, vertex_descriptor b, vertex_descriptor c
            ){ untested();
        assert(boost::is_directed(_g));

        // x has degree three connected to a,b,c. not these
        assert(!boost::edge(u, x, _g).second);
        assert(!boost::edge(v, x, _g).second);
        assert(!boost::edge(w, x, _g).second);

        // these are not connected due to the triangle rule
        assert(!boost::edge(u, v, _g).second);
        assert(!boost::edge(u, w, _g).second);
        assert(!boost::edge(v, w, _g).second);

        _marker.clear();
        auto p=mark_and_remove_helper(a, b, _marker, _g, enable_edge_cleanup);
        boost::remove_out_edge_if(u, p, _g);
        p._b = c;
        boost::remove_out_edge_if(v, p, _g);
        p._a = b;
        boost::remove_out_edge_if(w, p, _g);

        if(enable_edge_cleanup){ untested();
            clear_out_edges(x);
        }else{untested();
        }

        boost::add_edge(u, v, _g);
        boost::add_edge(u, w, _g);
        boost::add_edge(u, x, _g);
        boost::add_edge(v, w, _g);
        boost::add_edge(v, x, _g);
        boost::add_edge(w, x, _g);
        boost::add_edge(v, u, _g);
        boost::add_edge(w, u, _g);
        boost::add_edge(x, u, _g);
        boost::add_edge(w, v, _g);
        boost::add_edge(x, v, _g);
        boost::add_edge(x, w, _g);

        _degree[u] += 3; // hmm actually, +1?
        _degree[v] += 3;
        _degree[w] += 3;
        _degree[x] += 3;
        _num_edges += 6;

    }

    void addtoelims(vertex_descriptor v){
#ifndef NDEBUG
        if(_degs.is_reg(v)){
        }else{ untested();
            assert(!(disable_cube && disable_triangle && disable_buddy));
            //
        }
#endif
        _degs.unlink(v);
        assert(!_degs.is_reg(v));
        _elims.push_back(v);
        _numbering.put(v);
        _numbering.increment();
    }

protected: // rules
//    void Islet();
    void eliminate_vertex_1(vertex_descriptor v);
    void eliminate_vertex_2(vertex_descriptor v);
    bool Triangle(vertex_descriptor v);
    bool Buddy(vertex_descriptor v, vertex_descriptor w);
    bool Cube(vertex_descriptor v);
    bool BothSimplicial(vertex_descriptor v);
private: // graph update stuff.
    void increment_edges(long int n=1){_num_edges+=n;}
    void clear_vertex(vertex_descriptor v){ untested();
        incomplete();
        // does not work like this for vectors.
        boost::clear_vertex(v, _g);
    }
#ifndef nofilter
    std::pair<adjacency_iterator_filter, adjacency_iterator_filter>
    adjacent_vertices(vertex_descriptor v) const{
        auto p=boost::adjacent_vertices(v, _g);
        if(p.first==p.second){untested();
        }else{
        }
        return std::make_pair(adjacency_iterator_filter(p.first, _numbering, p.second),
                adjacency_iterator_filter(p.second, _numbering));
    }
#else
    std::pair<adjacency_iterator, adjacency_iterator>
    adjacent_vertices(vertex_descriptor v) const{ untested();
        auto p=boost::adjacent_vertices(v, _g);
        return p;
    }
#endif
private:
   // G_t _g;
    directed_view_type _g;
    idmap_type _id;
    degree_type _degree; // the vector. for now
    degreemap_type _degreemap;
    degs_type _degs;
    std::deque<vertex_descriptor> _elims;
    unsigned _lb_bs;
    unsigned _num_edges;
    marker_type _marker;
    treedec::draft::sMARKER<vertices_size_type, vertices_size_type> _dormant;
    numbering_type _numbering;
}; // preprocessing

// Check if there exists a degree-0-vertex.
template <typename G_t, typename B_t>
void Islet(G_t &G, B_t &bags, int &low_td)
{ untested();
    typedef typename treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){ untested();
        if(boost::degree(*vIt, G) == 0){ untested();
            typename treedec_traits<T_t>::vd_type vd=get_vd(G, *vIt);
            typename treedec_traits<T_t>::bag_type emptybag;

            auto n=boost::add_vertex(bags);
            boost::get(boost::vertex_underlying, n, bags) = vd;

            if(low_td<0){ untested();
                low_td = 0;
            }else{ untested();
            }
        }else{ untested();
        }
    }
}

// pick degree zero nodes into singleton bags of T.
template <typename G_t, typename T_t>
void Islet(G_t &G, T_t &bags)
{ untested();
    int low = -1;
    Islet(G, bags, low);
}

// check if a and b have the same neighbour set
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::check_twins_3(
        vertex_descriptor a, vertex_descriptor b) const
{
    assert(_degree[a]==3);
    assert(_degree[b]==3);
    auto pa = adjacent_vertices(a);
    auto& Ia=pa.first;
    auto pb = adjacent_vertices(b);
    auto& Ib=pb.first;
    bool ret;
    if(*Ia==*Ib){
        ++Ia; ++Ib;
        if(*Ia==*Ib){
            ++Ia; ++Ib;
            ret = *Ia==*Ib;
        }else{
            vertex_descriptor a=*Ia;
            ++Ia;
            if(*Ia==*Ib){
                ++Ib;
                ret = a==*Ib;
            }else{
                ret = false;
            }
        }
    }else{
        vertex_descriptor a=*Ia;
        ++Ia;
        if(*Ia==*Ib){
            // a=x
            // =xx
            ++Ib;
            if(*Ib == a){ untested();
                ++Ia;
                ++Ib;
                return Ia==Ib;
            }else{
                ++Ia;
                if(*Ia==*Ib){
                    ++Ib;
                    ret = a==*Ib;
                }else{ untested();
                    ret = false;
                }
            }
        }else{
            vertex_descriptor A=*Ia;
            ++Ia;
            // aA?
            // ?xx
            if(*Ia==*Ib){ untested();
                // aA=
                // =xx
                ++Ib;
                if(a==*Ib){ untested();
                    ++Ib;
                    ret = A==*Ib;
                }else if(A==*Ib){ untested();
                    ++Ib;
                    ret = a==*Ib;
                }else{ untested();
                    ret = false;
                }
            }else{
                ret = false;
            }
        }
    }
    // assert(ret==check_twins(a, b, _g));
    return ret;
}

// eliminate vertex: turn neighbours into clique, remove center, update degrees.
template<class G_t, template<class G_> class CFGT>
void preprocessing<G_t, CFGT>::eliminate_vertex_1(
        typename preprocessing<G_t, CFGT>::vertex_descriptor v)
{
    assert(_degree[v]==1);

    // queue for redegree
    assert(_degs.is_reg(v));
    auto f=adjacent_vertices(v).first;

    auto& df=_degree[*f];
    // _degs.unlink(*f, df);

    _num_edges -= 1;

    trace4("========= ex_1", v, *f, _degree[*f], _num_edges);
    addtoelims(v);
    assert(_degs.is_reg(*f));
    assert(_degs.is_reg(*f));
    assert(_degs.is_reg(*f));
    _degs.unlink(*f, df);
    assert(!_degs.is_reg(*f));
    --df;
    assert(!_degs.is_reg(*f));
    _degs.reg(*f, df);
    assert(_degs.is_reg(*f));

    assert(df==_degree[*f]);

    if(_lb_bs < 2){
        _lb_bs = 2;
    }else{
    }
}

template<class G_t, template<class G_> class CFGT>
void preprocessing<G_t, CFGT>::eliminate_vertex_2(
        typename preprocessing<G_t, CFGT>::vertex_descriptor v)
{

    trace1("========= ex_2", v);
    auto f=adjacent_vertices(v).first;
    auto x=*f;

    assert(2==_degree[v]);

    // isolate_(v);
    _num_edges -= 2;
    addtoelims(v);

    _marker.clear();
    _marker.mark(*f);
    bool need_edg=true;

    auto Is=adjacent_vertices(*(++f));
    for(; Is.first!=Is.second; ++Is.first){
        if(_marker.is_marked(*Is.first)){
            need_edg=false;
            break;
        }else{
        }
    }

    if(need_edg) {
        boost::add_edge(x, *f, _g);
        boost::add_edge(*f, x, _g);
        _num_edges+=1;

        // needed?
        wake_up_node(x);
        wake_up_node(*f);
    }else{
        // degree decreases by 1;
        --_degree[*f];
        --_degree[x];

        // needed?
        // actually, just x and *f
        wake_up_node(x);
        wake_up_node(*f);
    }

    if(_lb_bs < 3){
        _lb_bs = 3;
    }else{
    }
}

// Apply the Buddy rule if possible
// (checks if there exists two degree-3-vertices,
//  such that they share their neighbours)
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::Buddy(
        vertex_descriptor v, vertex_descriptor w)
{
    assert(_degree[v]==3);
    assert(_degree[w]==3);
    if(check_twins_3(v, w)){ untested();
        assert(!is_numbered(v));
        assert(!is_numbered(w));
        unlink_1_neighbourhood(v);
        _degs.unlink(w, 3);

        vd_type vd1 = get_vd(_g, v);
        vd_type vd2 = get_vd(_g, w);

        make_neigh_clique(v);
        assert(_degree[v]);
        addtoelims(vd2);
        isolate(w);

#ifdef NEGATIVE_TAGS3
        // TODO inefficent/ too many
        // which edges have been inserted?
        auto Is=adjacent_vertices(w);
        for(; Is.first!=Is.second; ++Is.first){ untested();
            auto It=adjacent_vertices(*Is.first);
            for(; It.first!=It.second; ++It.first){ untested();
                _dormant.unmark(*It.first);
                if(*It.first!=v){
                    _degs.reg(*It.first);
                }
            }
        }
#else

        // redegree and need-update
        redegree(v, 1);
#endif

        if(_lb_bs < 4){
            _lb_bs = 4;
        }else{
        }
        return true;
    }else{
        return false;
    }
}

template<class T, class I>
inline void rearrange_neighs(T* N, T x, I i)
{ untested();
    if(N[0] == x){ untested();
        N[0] = *(++i);
    } else if(N[1] == x){ untested();
        N[1] = *(++i);
    } else{ untested();
    }
}

//Apply the Cube rule if possible.
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::Cube(vertex_descriptor x)
{

    assert(_degree[x]==3);

    if(disable_triangle){
        incomplete();
        return false;
    }else{
    }
    vertex_descriptor a, b, c;

    auto f=adjacent_vertices(x).first;
    a = *f;

    if(_degree[a]!=3){ untested();
        return false;
    }else{
    }
    b = *(++f);
    if(_degree[b]!=3){ untested();
        return false;
    }else{
    }
    c = *(++f);
    if(_degree[c]!=3){ untested();
        return false;
    }else{
    }

    std::set<vertex_descriptor> cbag;
    assign_neighbours(cbag, a, b, c, _g);

    if(cbag.size() != 4){
        return false;
    }else{ untested();
    }

    vertex_descriptor N[6];
    vertex_descriptor* Na=&N[0];
    vertex_descriptor* Nb=&N[2];
    vertex_descriptor* Nc=&N[4];

    f = adjacent_vertices(a).first;
    Na[0] = *f; Na[1] = *(++f);
    rearrange_neighs(Na, x, f);

    f = adjacent_vertices(b).first;
    Nb[0] = *f; Nb[1] = *(++f);
    rearrange_neighs(Nb, x, f);

    f = adjacent_vertices(c).first;
    Nc[0] = *f; Nc[1] = *(++f);
    rearrange_neighs(Nc, x, f);


    typename boost::graph_traits<G_t>::vertex_descriptor u, v, w;

    if(Na[0] == Nb[0]){ untested();
        u = Na[0]; v = Na[1]; w = Nb[1];
    }else if(Na[0] == Nb[1]){ untested();
        u = Na[0]; v = Na[1]; w = Nb[0];
    }else if(Na[1] == Nb[0]){ untested();
        u = Na[1]; v = Na[0]; w = Nb[1];
    }else if(Na[1] == Nb[1]){ untested();
        u = Na[1]; v = Na[0]; w = Nb[0];
    }else{ untested();
        return false;
    }

    if(  (Nc[0] == v && Nc[1] == w)
      || (Nc[1] == v && Nc[0] == w)){ untested();
        assert(boost::edge(a, u, _g).second);
        assert(boost::edge(a, v, _g).second);
        assert(boost::edge(a, x, _g).second);

        addtoelims(a);
        addtoelims(c);
        addtoelims(b);

        // faster.
        isolate(a);
        isolate(b);
        isolate(c);

        cube_make_clique(u, v, w, x,
                         a, b, c);

        wake_up_node(u);
        wake_up_node(v);
        wake_up_node(w);
        wake_up_node(x); // required?? should be wake.
                         // also: simplicial! eliminate.

        wake_up_neighs(u);
        wake_up_neighs(v);
        wake_up_neighs(w);
        // wake_up_neighs(x); u,v,w

        if(_lb_bs < 4){
            _lb_bs = 4;
        }else{
        }

        return true;
    }else{ untested();
        return false;
    }
}

#if 0 // outdated.
//Apply the Simplicial rule, if possible
// (checks if the neighbors of v induce a clique.)
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::Simplicial(vertex_descriptor v)
{
    auto p=adjacent_vertices(v);
    bool isClique=true;
#if 0
    //The neighbourhood of v is a clique, if no "edge miss" occurs.

    adjacency_iterator nIt2;
    auto next=p.first;

    for(; p.first!=p.second; p.first=next){
        ++next;
        nIt2 = next;
        for(; nIt2!=p.second; ++nIt2){
            if(!boost::edge(*p.first, *nIt2, _g).second){
                isClique = false;
                goto DOUBLE_BREAK;
            }
        }
    }
    DOUBLE_BREAK:
#else
    // hmm, maybe slower on sets?
    // (do we want to keep sets?)

    _marker.clear();
    vertices_size_type cnt=0;

    // mark neighbours one aftyer another
    // check if marked neighbours are connected to previous neighs.
    if( p.first!=p.second){
        _marker.mark(*p.first);
        ++p.first;
    }
    for(; p.first!=p.second; ++p.first){
        ++cnt;
        unsigned missing = cnt;

        auto q=adjacent_vertices(*p.first);
        for(; q.first!=q.second; ++q.first){
            if(_marker.is_marked(*p.first)){
                --missing;
            }else{
            }
        }
        if(missing){ untested();
            isClique = false;
            break;
        }else{ untested();
        }
        _marker.mark(*p.first);
    }
#endif


    if(isClique){ untested();
        vd_type vd = get_vd(_g, v);

        unlink_1_neighbourhood(v, _g, _degs);
        _degs.unlink(v);

        auto deg=_degree[v];
        make_neigh_clique(v);
        redegree(v);

        if (unsigned(_lb_tw) < deg){ untested();
            _lb_tw = deg;
        }else{ untested();
        }

        return true;
    }else{ untested();
        return false;
    }
}

//Apply the Almost Simplicial rule if possible.
//same as Simplicial, but specialNeighbour
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::AlmostSimplicial(vertex_descriptor v)
{
    bool isAlmostSimplicial = true;
    bool specialNeighbourFound = false;
    adjacency_iterator next, nIt1, nIt2;
    vertex_descriptor cand1, cand2, specialNeighbour;
    unsigned int missingEdgesCount;

    auto p=adjacent_vertices(v);
    next = p.first;
    for(; p.first!=p.second; nIt1=next){
        nIt2 = ++next;
        missingEdgesCount = 0;
        for(; nIt2!=p.second; ++nIt2){
            if(!specialNeighbourFound){
            }else if(*nIt1 == specialNeighbour){
                continue;
            }else if(*nIt2 == specialNeighbour){
                continue;
            }else{
            }

            if(!boost::edge(*nIt1, *nIt2, _g).second){
                if(specialNeighbourFound){
                    //#special neighbours > 1.
                    isAlmostSimplicial = false;
                    goto DOUBLE_BREAK;
                }
                //*nIt1 or *nIt2 is a special neighbour.
                cand1 = *nIt1;
                cand2 = *nIt2;
                missingEdgesCount++;
            }
        }

        if(missingEdgesCount > 0){
            if(missingEdgesCount == 1){
                //cand2 has to be the special neighbour.
                specialNeighbour = cand2;
            }
            else{
                //cand1 has to be the special neighbour.
                specialNeighbour = cand1;
            }
            specialNeighbourFound = true;
        }
    }

    return isAlmostSimplicial;

    DOUBLE_BREAK:

    if(isAlmostSimplicial){
        vertices_size_type deg_v = _degree[v];
        assert(deg_v);
        assert(_lb_tw>=0);

        if(deg_v <= (vertices_size_type)_lb_tw){
            vd_type vd = get_vd(_g, v);

            unlink_1_neighbourhood(v, _g, _degs);
            _degs.unlink(v);

            make_neigh_clique(v);
            assert(_degree[v]);
            redegree(v);

            return true;
        }else if(vertices_size_type(_lb_tw) < deg_v-1u){ untested();
            _lb_tw = deg_v-1;
            return true;
        }else{
            return false;
        }
    }
    return false;
}
#endif

// Simplicial and AlmostSimplicial in one go.
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::BothSimplicial(vertex_descriptor v)
{
    assert(!boost::edge(v, v, _g).second);
    if(disable_almost_simplicial && disable_simplicial){
        return false;
    }else{
    }
    vertices_size_type deg=_degree[v];
    unsigned cnt=deg-1;
    unsigned missing;
    typename std::make_signed<edges_size_type>::type balance=0;

    _marker.clear();

    // mark neighbours one aftyer another
    // check if marked neighbours are connected to previous neighs.
    auto pp=adjacent_vertices(v);
    for(; pp.first!=pp.second; ++pp.first){
        _marker.mark(*pp.first);
    }

    auto p=adjacent_vertices(v);
    vertex_descriptor special;
    for(; p.first!=p.second; ++p.first){
        trace2("neighbourhood", v, *p.first);
        assert(cnt+3>deg);
        missing = cnt;

        auto q=adjacent_vertices(*p.first);
        for(; q.first!=q.second; ++q.first){
            if(_marker.is_marked(*q.first)){
                --missing;
            }else{
            }
        }

        assert(missing<=boost::num_vertices(_g));

        if(missing>1){
            if(balance>0){
                balance=-1;
                // that would be the second
                // special node.
                break;
            }
            balance += (missing+deg);
            special = *p.first;
        }else if(missing==1){
            if(!balance){
                special = *p.first;
            }else{
            }
            --balance;
        }else{
            // not missing.
        }
        trace2("missing", missing, balance);
    }


    if(!disable_simplicial
     && !balance){
        // did not find any missing.
        // it's a clique!
        vd_type vd = get_vd(_g, v);

        // hmm redegree can be faster.
        // unlink_1_neighbourhood(v);
        auto pp=adjacent_vertices(v);
        for(; pp.first!=pp.second; ++pp.first){
            auto n=*pp.first;
            if(0 && enable_mass_elimination){
                incomplete();
                // degree update missing.
                if(_degree[*pp.first]==deg-1){
                    _numbering.put(n);
                    _numbering.increment();
                    _degs.unlink(n);
                    _num_edges-=(deg-1);
                    // eliminate
                }else{
                    wake_up_node(n);
                }
            }else{
            }
        }

        addtoelims(v);
        // isolate
        auto p=adjacent_vertices(v);
        for(; p.first!=p.second; ++p.first){
            auto n=*p.first;
            assert(n!=v);
            assert(!_numbering.is_numbered(n));
            --_degree[n];
            wake_up_node(n);
        }
        _num_edges -= _degree[v];

        if (_lb_bs < 1+deg){ untested();
            _lb_bs = 1+deg;
        }else{
        }

        return true;
    }else if(disable_almost_simplicial){ untested();
        return false;
    }else if( balance==-2 || edges_size_type(balance)==deg ){ untested();
        trace1(">>>>>>>>>> AlmostSimplicial", v);
        // missing and multimissings are balanced.
        // it's almost simplicial.
        assert(deg);
        assert(_lb_bs>=1);

        if(deg+1 <= _lb_bs){ untested();
            vd_type vd = get_vd(_g, v);

#if 0
            make_neigh_clique(v);
#else
            isolate_(v); // marks all neighbours of v, numbers v
                         // "removes" edges to v
                         //
            assert(_numbering.is_numbered(v));
            unsigned check=0;
            auto Is=adjacent_vertices(special);
            for(; Is.first!=Is.second; ++Is.first){ untested();
                _marker.unmark(*Is.first);
                trace1("unmark sneigh", *Is.first);
                if(*Is.first==v){ untested();
                }else if(is_dormant(*Is.first)){
                    wake_up_node(*Is.first);
                }else{
                    _degs.update(*Is.first);
                }
            }
            _marker.unmark(special);
            // edge from still marked v-neighbours to special
            trace3("need edges", check, balance, _degree[v]);
            auto Iv=adjacent_vertices(v);
            for(; Iv.first!=Iv.second; ++Iv.first){ untested();
                auto n=*Iv.first;
                trace2("adj v", n, _marker.is_marked(n));
                if(_marker.is_marked(n)){ untested();
                    assert(n!=special);
                    assert(!boost::edge(n, special, _g).second);
                    boost::add_edge(n, special, _g);

                    assert(!boost::edge(special, n, _g).second);
                    boost::add_edge(special, n, _g);
                    ++_degree[n];
                    ++_degree[special];
                    increment_edges();
                }else{ untested();
                }
            }

#endif
            assert(_degree[v]);
            redegree(v);

            return true;
        }else if(_lb_bs < deg){ untested();
            _lb_bs = deg;
            return true;
        }else{ untested();
            return false;
        }
    }else{ itested();
        return false;
    }
}

//Apply the Triangle rule if applicable (checks if there exists a
//degree-3-vertex, such that at least one edge exists in its neighbourhood).
//return true, if degs has been modified.
template<class G_t, template<class G_> class CFG>
bool preprocessing<G_t, CFG>::Triangle(vertex_descriptor v)
{
    vertices_size_type deg=_degree[v];
    assert(deg>=3);
    vertex_descriptor N[3];
    auto f=adjacent_vertices(v).first;
    N[0] = *f;
    N[1] = *(++f);
    N[2] = *(++f);

    bool have_edg=false;

    if(boost::edge(N[0], N[1], _g).second){
        have_edg=true;
    }else if(boost::edge(N[0], N[2], _g).second){
        have_edg=true;
        std::swap(N[1], N[2]);
    }else if(boost::edge(N[1], N[2], _g).second){
        have_edg=true;
        std::swap(N[0], N[2]);
    }else{
    }

    if(have_edg) {
        trace4("edg", deg, v, N[0], N[1]);
        assert(_degs.is_reg(v));

        // installs at most 2 edges.
        make_neigh_clique(v);
        assert(!_degs.is_reg(v));

        wake_up_neighs(N[0]);
        wake_up_neighs(N[1]);
        wake_up_neighs(N[2]);

        if(_lb_bs<4){
           _lb_bs = 4;
        }else{
        }
        assert(!_degs.is_reg(v));
        return true;
    }else{
        // there is no edge.
        return false;
    }

    return false;
}


template<class G_t, template<class G_> class CFG>
void preprocessing<G_t, CFG>::do_it()
{
    typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(_g);

    if(num_vert == 0){ untested();
        return;
    }

    const degs_type& cdegs(_degs);

    //Islet rule
    trace1("", cdegs.size());
    assert(cdegs.size());
    if(!cdegs[0].empty()){
        if (_lb_bs==0){
            _lb_bs = 1;
        }else{ untested();
        }
    }else{
    }

    auto const& B=cdegs[0];
    auto I=B.begin();
    auto E=B.end();
    for(; I!=E; ++I){
        _elims.push_back(*I);
        // no need to number...
    }

    unsigned min_ntd = 1;
    while(_num_edges)
    {
        trace2("", boost::num_edges(_g), _num_edges);
        if(min_ntd>1){
            --min_ntd;
        }else{
        }
        assert(min_ntd);

        vertex_descriptor v;
        boost::tie(v, min_ntd) = _degs.pick_min(min_ntd, num_vert);
        assert(!boost::edge(v, v, _g).second);
        trace2("pp main loop ================ ", v, min_ntd);
        assert(treedec::is_valid(v, _g));
        assert(_degree[v] == min_ntd);
        assert(_numbering.is_not_numbered(v));

        if(min_ntd==1){
            auto f=adjacent_vertices(v).first;
            *f;
#ifndef NDEBUG
            trace3("", v, *f, _degs.is_reg(*f));
#endif
            eliminate_vertex_1(v);
            continue;
        }else if(min_ntd==2){
            auto f=adjacent_vertices(v).first;
            ++f;
            *f;
            eliminate_vertex_2(v);
            assert(_degs.is_reg(*f));
            continue;
        }else if(min_ntd==3){
            //degree 3-rules
            auto const& B=cdegs[3];
            unsigned cnt=0;
            auto it1=B.begin();
            auto next=it1;
            for(; it1!=B.end(); it1=next){
                ++next;
                ++cnt; // needed later for Cube.
#ifdef NEGATIVE_TAGS
                if(_dormant.is_marked(*it1)){
                    incomplete();
                    // continue;
                }else{untested();
                }
#endif
                //Triangle
                if(disable_triangle){ untested();
                }else if(Triangle(*it1)){
                    trace1("==============did triangle", *it1);
                    goto NEXT_ITER;
                }else{
                    // graph is unchanged.
                }
                //Buddy
                auto it2=next;
                // check: does buddy require two wake nodes?
                // (probably)
                for(; it2!=B.end(); ++it2){
                    trace2("buddyloop", *it1, *it2);
                    assert(*it1!=*it2);
#ifdef NEGATIVE_TAGS
                    if(_dormant.is_marked(*it2)){
                        incomplete();
                        // continue;
                    }else{untested();
                    }
#endif
                    if(disable_buddy){ untested();
                    }else if(Buddy(*it1, *it2)){ untested();
                        trace1("==============buddy", *it1);
                        goto NEXT_ITER;
                    }else{
                        // graph is unchanged.
                    }
                }
#ifdef NEGATIVE_TAGS
                if(cnt<4){
                    // cannot do anything about this degree 3 node.
                    // (no cube rule will change that)
                    set_dormant(*it1);
                    // puh
                    // _degs.unlink(*it1);
                }
#endif
            }
            it1=B.begin();
            if(cnt>3)
            for(; it1!=B.end(); ++it1){
                if(disable_cube){ untested();
                }else if(Cube(*it1)){ untested();
                    trace1("==============Cube", *it1);
                    goto NEXT_ITER;
                }else{
                    // graph is unchanged.
                }

#ifdef NEGATIVE_TAGS
                // cannot do anything about this degree 3 node.
                _dormant.mark(*it1);
#endif
            }
            goto ARBITRARY_DEGREE;
        }else{
            ARBITRARY_DEGREE:
                ;
        }
        {

            if (_lb_bs < 5){
                _lb_bs = 5;
            }else{
            }

            for(unsigned int i = min_ntd; i < num_vert; ++i){
                auto const& B=cdegs[i];
                auto it=B.begin();
                for(; it != B.end(); ++it){
                    if(_dormant.is_marked(*it)){untested();
                        unreachable();
                    }else{
                    }

                    if(disable_simplicial && disable_almost_simplicial){
                    }else if(BothSimplicial(*it)){
                        goto NEXT_ITER;
                    }else{
                    }
#if 0
                    if(disable_simplicial){
                    }else if(Simplicial(*it)){
                        std::cout<< "s" << *it << "\n";
                        std::cout<< "d" << _degree[*it] << "\n";
                        std::cout<< "d" << boost::degree(*it, _g) << "\n";
                      //  assert(false);
                        goto NEXT_ITER;
                    }else if(disable_almost_simplicial){
                    }else if(AlmostSimplicial(*it)){
                        std::cout<< "a" << *it << "\n";
                        std::cout<< "d" << _degree[*it] << "\n";
                        std::cout<< "d" << boost::degree(*it, _g) << "\n";
                      //  assert(false);
                        goto NEXT_ITER;
                    }else{
                        std::cout << "nope\n";
                    }
#endif
                    set_dormant(*it);
                }
            }
        }
        return;
NEXT_ITER:
        ;
    } // main loop
} // pp::do_it

} //namespace impl

template <typename G_t, typename BV_t>
void preprocessing(G_t &G, BV_t &bags, int &low)
{
    if(boost::num_vertices(G)){
        impl::preprocessing<G_t> A(G);
        A.set_treewidth(low, -1u);
        A.do_it();
        low = A.get_treewidth();
        // obsolete interface. possibly slow
        A.get_bags(bags);
        A.get_graph(G);
    }else{
    }
}

template <typename G_t, typename BV_t>
void preprocessing(G_t &G, BV_t &bags)
{ untested();
    if(boost::num_vertices(G)){ untested();
        impl::preprocessing<G_t> A(G);
        A.do_it();
        // obsolete interface. possibly slow
        A.get_bags(bags);
        A.get_graph(G);
    }else{ untested();
        // it does not work (yet.)
    }
}


} //namespace treedec

#endif //TD_PREPROCESSING

// vim:ts=8:sw=4:et
