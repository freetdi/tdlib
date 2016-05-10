// Felix Salfelder 2016
// Lukas Larisch 2016
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

#ifndef TD_FILL_HPP
#define TD_FILL_HPP

#include <boost/graph/graph_traits.hpp>
#include <assert.h>

#include "graph.hpp"

#if __cplusplus >= 201103L
# include <unordered_set>
#endif

// temporary.
#ifndef untested
#define untested()
#endif
#ifndef itested
#define itested()
#endif
#ifndef incomplete
#define incomplete()
#endif


namespace misc {

namespace detail {
// FIXME: not here
template<class G>
struct fill_config{
    typedef typename boost::graph_traits<G>::vertex_descriptor vd_type;
#if __cplusplus < 201103L
    typedef std::set<vd_type> bag_type;
#else
    typedef std::unordered_set<vd_type> bag_type;
#endif
    // typedef stx::btree_set<vd_type> bag_type;
    static void alloc_init(size_t){
    }
    static unsigned num_threads(){return 1;}
};
} // detail
} // misc


namespace misc{ // treedec::misc? hmmm.

template<class G_t, class CFG=detail::fill_config<G_t> >
class FILL{
    FILL(const FILL&)
    { untested();
    }
public:
#if 0
    class firstless{
        public:
            bool operator()( pair_t a, pair_t b){
                if(a.first<b.first){
                    return true;
                }else{
                    return false;
                }
            }
    };
#endif

public: // types
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename CFG::bag_type bag_type;
    typedef typename bag_type::iterator bag_iterator;
    //typedef std::vector<bag_type> container_type;
    typedef std::set<std::pair<size_t, vertex_descriptor> > container_type;
    // typedef typename container_type::iterator iterator;
    // typedef typename container_type::const_iterator const_iterator;
    typedef typename boost::graph_traits<G_t>::vertices_size_type fill_t;

private:
    size_t max_missing_edges() const
    {
        // isolated nodes?
        size_t nv=boost::num_vertices(_g);
        return (nv-1)*(nv-2);
    }

public: // construct
    FILL(const G_t& g): _g(g)
    {
        _vals.resize(boost::num_vertices(g));
     //   CFG::alloc_init(boost::num_vertices(g));
        vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
            if(boost::degree(*vIt, g) > 0){ //skip isolated vertices
                reg(*vIt);
            }
        }
    }

public: // queueing
    void unlink(const vertex_descriptor& v, size_t f)
    {
        int n=_fill.erase(std::make_pair(f,v));
        (void)n;
        assert(n==1);
    }
    void unlink(const vertex_descriptor& v)
    {
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        unlink(v, _vals[pos]);
    }

public:
    void reg(const vertex_descriptor& v, size_t missing_edges)
    {
        bool n=_fill.insert(std::make_pair(missing_edges,v)).second;
        assert(n);
        (void)n;

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos] = missing_edges;
    }
public:
    void reg(const vertex_descriptor v)
    {
        size_t missing_edges=treedec::count_missing_edges(v,_g);
        reg(v, missing_edges);
    }
    void q_decrement(const vertex_descriptor v)
    { incomplete();
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        unlink(v, _vals[pos]);
        reg(v, --_vals[pos]);
    }

public: // picking
    // vertex_descriptor pick(unsigned fill)
    // {
    //     return *_fill[fill].begin();
    // }
    // pick a minimum fill vertex within fill range [lower, upper]
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower=0, unsigned upper=-1) const
    {
        assert(lower==0); // for now.

        BOOST_AUTO(b, _fill.begin());
        return std::make_pair(b->second, b->first);
    }
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower, unsigned upper, bool erase)
    {
        vertex_descriptor p = pick_min(lower,upper);
        if(erase){ untested();
            unlink(p.first,p.second);
        }
        return p;
    }

#if 0
    size_t num_nodes() const
    { untested();
        unsigned N=0;
        for(const_iterator i=_fill.begin(); i!=_fill.end(); ++i) { itested();
            N+=i->size();
        }
        return N;
    }
#endif

#if 0
    bag_type const& operator[](size_t x) const
    {
        return _fill[x];
    }
    size_t size() const
    {
        return _fill.size();
    }
#endif

#if 0
public:
    bag_type& operator[](size_t x)
    {
        return _fill[x];
    }
#endif

private:
    const G_t& _g;
//private: // later.
    container_type _fill;
    std::vector<int> _vals;
}; // FILL

} //namespace misc

namespace noboost{

template<class G_t>
struct fill_chooser{
    typedef typename misc::FILL<G_t> type;
    typedef type fill_type; // transition? don't use.
};

}

#endif //guard

// vim:ts=8:sw=4:et
