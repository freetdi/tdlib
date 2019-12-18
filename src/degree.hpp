// Felix Salfelder 2016
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
// Keeping track of degrees.
//


#ifndef TREEDEC_DEGREE_HPP
#define TREEDEC_DEGREE_HPP

// use local copy
#include "trace.hpp"
#include "config_traits.hpp"
#include "bucket_sorter_bits.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <assert.h>

#include "directed_view.hpp"
//
#include "container.hpp"
#include "degree_config.hpp"
#include "random_generators.hpp" // BUG. see below

// #if __cplusplus >= 201103L
// # include <unordered_set>
// #include <tlx/container/btree_set>
// #endif

#include "platform.hpp"
#include "graph_traits.hpp"

#include <stack>

#ifdef HAVE_TLX_CONTAINER_BTREE_SET_HPP
# include <tlx/container/btree_set.hpp>
#endif


// BUG
namespace misc {

namespace detail {


// FIXME: wrong place
// needs random_generators.hpp
// BUG: used by svbs_random.
template<class G_t>
struct random_deg_config : public treedec::degs::default_config<G_t>{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd_type;

    // FIXME: need to choose bucket backend...
    // (the default does not implement back access)

    template <typename C_t>
    static vd_type pick(C_t const &C){
        bool c = treedec::random::coin();
        if(c){
            return *C.begin();
        }else{ itested();
            return treedec::back(C);
        }
    }
};


} // namespace detail

template<class G_t, template<class G> class CFGT=detail::deg_config >
class DEGS{
    DEGS(const DEGS&)
    { untested();
    }

public: // types
    typedef CFGT<G_t> CFG;
    typedef typename boost::graph_traits<G_t>::vertices_size_type size_type; //??
    typedef typename boost::graph_traits<G_t>::vertices_size_type degree_t; //??
    typedef typename boost::graph_traits<G_t>::vertices_size_type value_type;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    // typedef typename CFG::bag_type bag_type;
    typedef typename std::set<vertex_descriptor> bag_type; // legacy
    // typedef typename CFG::const_bag_type const_bag_type;
    typedef typename bag_type::iterator bag_iterator;
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::const_type idmap_type;
    // hmm, override?
    // typedef typename boost::property_map<G_t, boost::vertex_out_degree_t>::const_type degree_type;
    
    typedef boost::degree_property_map<G_t> fallback_degree_type;

    typedef typename boost::iterator_property_map<degree_t*, idmap_type, value_type, value_type&>
        degreemap_type;

#if 1
    typedef typename treedec::config::get::degree<CFG, fallback_degree_type>::type degree_type;
#else
    typedef degreemap_type degree_type;
#endif

    // TODO/BUG/FIXME?
    // bucket sorter keeps pos_to_vd map, although not always necessary...
    typedef boost::bucket_sorter<value_type,
            vertex_descriptor,
            degreemap_type,
           idmap_type > container_type;
    typedef typename container_type::const_stack internal_bag_type;

    //typedef typename container_type::iterator iterator;
    //typedef typename container_type::const_iterator const_iterator;

public: // construct
    DEGS(const G_t &g): _g(g),
         _vi(boost::get(boost::vertex_index, g)),
         _degree(_g),
         _vals(boost::num_vertices(g)),
         _degs(boost::num_vertices(g), // length
               boost::num_vertices(g) /*-1?*/,  // max_bucket
               boost::make_iterator_property_map(&_vals[0], _vi, size_type()),
               _vi)
    {
        // delegate construction? how?
        init();
    }
    DEGS(const G_t &g, degree_type& d): _g(g),
         _vi(boost::get(boost::vertex_index, g)),
         _degree(d),
         _vals(boost::num_vertices(g)),
         _degs(boost::num_vertices(g), // length
               boost::num_vertices(g) /*-1?*/,  // max_bucket
               boost::make_iterator_property_map(&_vals[0], _vi, size_type()),
               _vi)
    {
        init();
    }

private:
    void init(){

        if(!boost::num_vertices(_g)){
        }
        CFG::alloc_init(boost::num_vertices(_g));
        vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(_g); vIt != vEnd; ++vIt){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *vIt);
            assert(pos<_vals.size());
            _vals[pos] = boost::out_degree(*vIt, _g);
            trace2("deginit", pos, _vals[pos]);
            _degs.push(*vIt);
            assert(is_reg(*vIt));
        }
    }

public:
    DEGS& operator=(const DEGS& d)
    { untested();
        assert(&d._g==&_g);
        _degs = d._degs;
        return *this;
    }
    DEGS& operator=(const DEGS&& d)
    { untested();
        assert(&d._g==&_g);
        _degs = MOVE(d._degs);
        return *this;
    }

public: // queueing
    void unlink(const vertex_descriptor& v, size_t /*x*/)
    {
        assert(treedec::is_valid(v, _g));
#ifndef NDEBUG
        if(!is_reg(v)){
            incomplete();
            // bug in your code?
        }
#endif
        _degs.remove(v);
        // trace1("rmd", v);
        assert(!is_reg(v));
        assert(!_degs.is_known(v));
    }
    void unlink(const vertex_descriptor& v)
    {
        size_t d=boost::out_degree(v,_g);
        unlink(v,d);
    }

    void q_update(const vertex_descriptor& v)
    {
        unlink(v);
        _q.push(v);
    }
    void reg(const vertex_descriptor& v)
    {
        size_t d=boost::out_degree(v,_g);
        reg(v,d);
    }
    void reg(const vertex_descriptor& v, size_t d)
    {
        assert(treedec::is_valid(v,_g));
//        bool n=_fill.insert(std::make_pair(missing_edges,v)).second;
//        assert(n);
//        (void)n;

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos] = d;
        trace2("reg", v, pos);
        _degs.push(v);
    }
    void update(const vertex_descriptor& v)
    {
        _vals[v] = _degree[v]; // BUG. use the same array!!
        trace2("update", v, _vals[v]);
        _degs.update(v);
    }
#ifndef NDEBUG
    bool is_reg(const vertex_descriptor& v) const
    {
        return _degs.is_known(v);
    }
#endif

    void update_queued() {
        while(!_q.empty()){
            reg(_q.top());
            _q.pop();
        }
    }
    void flush() const {
    }

public: // picking
    vertex_descriptor pick(unsigned degree, bool erase=false)
    {
        if(erase){
            incomplete();
        }
        return CFG::pick(_degs[degree]);
    }
    bag_type const detach_bag(unsigned degree)
    {
#if 0 // pre-bucket
        auto& D = _degs[degree];
        bag_type B=MOVE(_degs[degree]);
        D.clear();
#else
        // this is inefficient...
        auto D = _degs[degree];
        bag_type B;
        while(!D.empty()){
            vertex_descriptor t=D.top();
            treedec::push(B, t);
            D.pop();
        }
#endif
        return B;
    }

    // pick a minimum degree vertex within degree range [lower, upper]
    std::pair<vertex_descriptor,degree_t> pick_min(unsigned lower=0, unsigned upper=-1)
    {
        while(_degs[lower].empty()){
            ++lower;
            // min_ntd==num_vert contradicts the outer loop condition
            // (this loop should be safe)
            assert(lower != upper+1);
        }
        vertex_descriptor min_nv;
        min_nv = CFG::pick(_degs[lower]);

        return std::make_pair(min_nv, lower);
    }
    std::pair<vertex_descriptor,degree_t> pick_min(unsigned lower, unsigned upper, bool erase)
    {
        auto p=pick_min(lower,upper);
        if(erase){
            unlink(p.first,p.second);
        }else{untested();
        }
        return p;
    }

#if 0
    size_t num_nodes() const{ untested();
        unsigned N=0;
        for(const_iterator i=_degs.begin(); i!=_degs.end(); ++i) { itested();
            N+=i->size();
        }
        return N;
    }
#endif

    void check()
    { // sometimes required when debugging fancy callbacks :/
#ifdef EXCESSIVE_DEG_DEBUG
            DEGS degs(_g);
            assert(_degs.size()==degs.size());
            assert(size()==boost::num_vertices(_g));
            assert(num_nodes()==degs.num_nodes());

            iterator j=degs._degs.begin();
            unsigned N=0;
            for(iterator i=_degs.begin(); i!=_degs.end();) {
                assert(N<boost::num_vertices(_g));
                unsigned I=i->size();
                unsigned J=j->size(); //actual _g

                if(I>J){
                    std::cerr<<"mismatch " << I << " " << J << "\n";
                    std::cerr<<"extra node " << *i->begin() << " of deg " << N << " in degs\n";
                }else if(I<J){
                    std::cerr<<"mismatch " << I << " " << J << " in " << N << "\n";
                    std::cerr<<"extra node " << *j->begin() << " of deg " << N << " in g\n";
                }
                assert(I==J);
                ++i;
                ++j;
                ++N;
            }
            assert(N==boost::num_vertices(_g));
#endif
    } //void check()

#if 1
    /// hmm better: iterator range?!
    internal_bag_type const operator[](size_t x) const
    {
        return _degs[x];
    }
    size_t size() const
    {
        return _degs.size();
    }
#endif
private:
    bag_type& operator[](size_t x)
    {
        return _degs[x];
    }

private:
    const G_t& _g;
    idmap_type _vi;
    degree_type _degree;
    std::vector<value_type> _vals;
    container_type _degs;
private:
    std::stack<vertex_descriptor> _q;
}; // DEGS

} //namespace misc

//register a 1-neighborhood to DEGS
template<class U, class G_t, class B, class D>
void redegree(U, G_t const &G, B I, B E, D& degree)
{
    for(; I!=E; ++I){
        typename boost::graph_traits<G_t>::vertex_descriptor x=*I;
        assert(treedec::is_valid(x, G));
        size_t deg = boost::out_degree(x, G);
        degree.reg(x, deg);
    }
}

template<class U, class G_t, class B, class D>
inline void redegree(U u, G_t const &G, std::pair<B,B> const& neighborhood, D& degree)
{
    auto I=std::begin(neighborhood);
    auto E=std::end(neighborhood);
    redegree(u, G, I, E, degree);
}

template<class U, class G_t, class B, class D>
inline void redegree(U u, G_t const &G, B const& neighborhood, D& degree)
{
    BOOST_AUTO(I, neighborhood.begin());
    BOOST_AUTO(E, neighborhood.end());
    redegree(u, G, I, E, degree);
}

// obsolete. don't use.
template<class G_t, class D>
void unlink_1_neighbourhood(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G, D &degs){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; ++nIt){
        degs.unlink(*nIt);
    }
}

namespace treedec{
    using misc::DEGS;
}

#endif // guard

// vim:ts=8:sw=4:et
