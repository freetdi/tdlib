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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//
// Keeping track of degrees.
//


#ifndef TD_DEGREE_HPP
#define TD_DEGREE_HPP

// use local copy
#include "trace.hpp"
#include "bucket_sorter_bits.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <assert.h>

#include "container.hpp"
#include "degree_config.hpp"
#include "random_generators.hpp" // BUG. see below

// #if __cplusplus >= 201103L
// # include <unordered_set>
// #include <stx/btree_set>
// #endif

#include "degree_config.hpp"
#include "platform.hpp"
#include "trace.hpp"
#include "graph_traits.hpp"

#include <stack>

#ifdef HAVE_STX_BTREE_SET_H
# include <stx/btree_set>
#endif

namespace misc {

namespace detail {


// FIXME: wrong place
// needs random_generators.hpp
// BUG: used by svbs_random.
template<class G_t>
struct random_deg_config : public misc::detail::deg_config<G_t>{
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

template<class G_t, class CFG=detail::deg_config<G_t> >
class DEGS{
    DEGS(const DEGS&)
    { untested();
    }

public: // types
    typedef typename boost::graph_traits<G_t>::vertices_size_type size_type;
    typedef typename boost::graph_traits<G_t>::vertices_size_type degree_t;
    typedef typename boost::graph_traits<G_t>::vertices_size_type value_type;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename CFG::bag_type bag_type;
    // typedef typename CFG::const_bag_type const_bag_type;
    typedef typename bag_type::iterator bag_iterator;
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::const_type idmap_type;
    typedef typename boost::iterator_property_map<degree_t*, idmap_type, value_type, value_type&>
        degreemap_type;

    // TODO/BUG/FIXME?
    // bucket sorter keeps pos_to_vd map, although not always necessary...
    typedef boost::bucket_sorter<value_type,
            vertex_descriptor,
            degreemap_type,
           idmap_type > container_type;
    typedef typename container_type::const_stack internal_bag_type;

    //typedef std::vector<bag_type> container_type;
    //typedef typename container_type::iterator iterator;
    //typedef typename container_type::const_iterator const_iterator;

public: // construct
    DEGS(const G_t &g): _g(g),
         _vi(boost::get(boost::vertex_index, g)),
         _vals(boost::num_vertices(g)),
         _degs(boost::num_vertices(g), // length
               boost::num_vertices(g) /*-1?*/,  // max_bucket
               boost::make_iterator_property_map(&_vals[0], _vi, size_type()),
               _vi)
    {
        if(!boost::num_vertices(g)){
        }
        CFG::alloc_init(boost::num_vertices(g));
        vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){ itested();
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *vIt);
            assert(pos<_vals.size());
            _vals[pos] = boost::degree(*vIt, g);
            trace2("deginit", pos, _vals[pos]);
            _degs.push(*vIt);
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
    void unlink(const vertex_descriptor& v, size_t)
    {
        assert(treedec::is_valid(v, _g));
        _degs.remove(v);
    }
    void unlink(const vertex_descriptor& v)
    {
        size_t d=boost::degree(v,_g);
        unlink(v,d);
    }

    void q_update(const vertex_descriptor& v)
    {
        unlink(v);
        _q.push(v);
    }
    void reg(const vertex_descriptor& v)
    {
        size_t d=boost::degree(v,_g);
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
    { untested();
        _degs.update(v);
    }

    void update_queued()
    { untested();

        while(!_q.empty()){
            reg(_q.top());
            _q.pop();
        }
    }
    void flush() const
    {
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
        // trace4("", lower, _degs[lower].top(), *_degs[lower].begin(), CFG::pick(_degs[lower]));
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
//    std::vector<status_t> _vals;
    const G_t& _g;
    idmap_type _vi;
    std::vector<value_type> _vals;
    container_type _degs;

private:
    std::stack<vertex_descriptor> _q;
}; // DEGS

} //namespace misc

//register a 1-neighborhood to DEGS
template<class U, class G_t, class B, class D>
void redegree(U, G_t &G, B& neighborhood, D& degree)
{
    BOOST_AUTO(I, neighborhood.begin());
    BOOST_AUTO(E, neighborhood.end());

    for(; I != E ; ++I){

        typename boost::graph_traits<G_t>::vertex_descriptor x=*I;
        assert(treedec::is_valid(x, G));
        size_t deg = boost::degree(x, G);
        degree.reg(x, deg);
    }
}

// obsolete. don't use.
template<class G_t, class D>
void unlink_1_neighbourhood(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G, D &degs){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; ++nIt){
        degs.unlink(*nIt);
    }
}


#endif // guard

// vim:ts=8:sw=4:et
