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

template <typename G_t>
inline size_t get_missing_edges_count(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const &G)
{
    size_t missing_edges = 0;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(v, G); nIt1 != nEnd; nIt1++){
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            if(!boost::edge(*nIt1, *nIt2, G).second){
                ++missing_edges;
            }
        }
    }
    return missing_edges;
}

template<class G_t, class CFG=detail::fill_config<G_t> >
class FILL{
    FILL(const FILL&)
    { untested();
    }

public: // types
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename CFG::bag_type bag_type;
    typedef typename bag_type::iterator bag_iterator;
    typedef std::vector<bag_type> container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef typename boost::graph_traits<G_t>::vertices_size_type fill_t;

public: // construct
    FILL(const G_t& g): _fill(boost::num_vertices(g)*boost::num_vertices(g)), _g(g)
    {
        _vals.resize(boost::num_vertices(g));
        CFG::alloc_init(boost::num_vertices(g));
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
        int n=_fill[f].erase(v);
        (void)n;
        //assert(n==1); // oops?
    }
    void unlink(const vertex_descriptor& v)
    {
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        unlink(v, _vals[pos]);
    }

private:
    void reg(const vertex_descriptor& v, size_t missing_edges)
    {
        bool n=_fill[missing_edges].insert(v).second;
        // assert(n); // BUG? multiple regs!
        (void)n;

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos] = missing_edges;
    }
public:
    void reg(const vertex_descriptor& v)
    {
        size_t missing_edges = get_missing_edges_count(v, _g);
        reg(v, missing_edges);
    }

public: // picking
    vertex_descriptor pick(unsigned fill)
    {
        return *_fill[fill].begin();
    }
    // pick a minimum fill vertex within fill range [lower, upper]
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower=0, unsigned upper=-1) const
    {
        while(_fill[lower].empty()){
            ++lower;
            assert(lower != upper+1);
        }
        vertex_descriptor min_nv = *_fill[lower].begin();
        return std::make_pair(min_nv, lower);
    }
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower, unsigned upper, bool erase)
    {
        vertex_descriptor p = pick_min(lower,upper);
        if(erase){ untested();
            unlink(p.first,p.second);
        }
        return p;
    }

    size_t num_nodes() const{ untested();
        unsigned N=0;
        for(const_iterator i=_fill.begin(); i!=_fill.end(); ++i) { itested();
            N+=i->size();
        }
        return N;
    }

    bag_type const& operator[](size_t x) const
    {
        return _fill[x];
    }
    size_t size() const
    {
        return _fill.size();
    }

public:
    bag_type& operator[](size_t x)
    {
        return _fill[x];
    }

//private: // later.
    container_type _fill;
    std::vector<int> _vals;

private:
    const G_t& _g;
}; // FILL

} //namespace misc


#endif //guard

// vim:ts=8:sw=4:et
