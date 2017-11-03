// Felix Salfelder, 2017
//
// (c) 2017 Felix Salfelder
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
#ifndef TREEDEC_GRAPH_UTIL_HPP
#define TREEDEC_GRAPH_UTIL_HPP

#include "marker_util.hpp"

namespace treedec {

namespace util {

namespace impl {

template<class G, class M>
class edge_map {
public:
	typedef typename boost::graph_traits<G>::vertex_descriptor vd;
	typedef typename M::value_type value_vertex_type;
public: // construct
	edge_map( G const& g, M const& m) : _m(m), _g(g) {}

	std::pair<value_vertex_type, value_vertex_type>
		operator()(typename boost::graph_traits<G>::edge_descriptor e) const {
			auto p=std::make_pair(boost::get(_m, boost::source(e, _g)),
					boost::get(_m, boost::target(e, _g)));
			trace2("edgmap", p.first, p.second);
			return p;
		}

private:
//	boost::iterator_property_map<value_vertex_type*,
//	                             VertexIndexMap,
//										  vd, value_vertex_type&> vd_map;
	M const& _m;
	G const& _g;
};

} // impl

template<class G, class M>
impl::edge_map<G, M> make_edge_map(M const& m, G const& g)
{
	return impl::edge_map<G, M>(g, m);
}

namespace detail{

// template<class S>
// struct visited_mask{
// static_assert(sizeof(S)==0);
// };

template<class S>
struct visited_mask{
	visited_mask(const visited_mask& o) : _s(o._s) {}
	visited_mask(S& s) : _s(s) {}
	bool operator()(unsigned x) const{
		assert(x<size());
		return _s[x];
	}
	bool operator[](unsigned x) const{
		return operator()(x);
	}
	void visit(unsigned x){
		assert(x<size());
		_s[x] = true;
	}
	size_t size() const{
		return _s.size();
	}
	S& _s;
};

} // detail

inline
detail::visited_mask<std::vector<BOOL> > make_incidence_mask(std::vector<BOOL>& v)
{
	return detail::visited_mask<std::vector<BOOL> >(v);
}

} //util

// count number of edges missing in 1-neighborhood of v
template <typename G_t>
inline size_t count_missing_edges(
        const typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const &G)
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
} // count_missing_edges

// same but use marker.
template <typename G_t, class MARKER>
inline size_t count_missing_edges(
        const typename boost::graph_traits<G_t>::vertex_descriptor v,
		  MARKER& marker, G_t const &g)
{
	size_t missing_edges = 0;

	auto p=boost::adjacent_vertices(v, g);
	for(; p.first!=p.second; ++p.first){
//		trace2("visit", v, *p.first);
		marker.clear();
		mark_neighbours(marker, *p.first, g);

		auto q=adjacent_vertices(v, g);
		for(; q.first!=q.second; ++q.first){
			if(*q.first>=*p.first){
				// skip. TODO: more efficient skip
			}else if(marker.is_marked(*q.first)){
				// done, reachable from *p.first
			}else{
//				trace2("found", *p.first, *q.first);
				++missing_edges;
			}
		}
	}
//	trace2("counted_missing_edges w/marker", v, missing_edges);
	return missing_edges;
} // count_missing_edges

} // treedec

#endif
