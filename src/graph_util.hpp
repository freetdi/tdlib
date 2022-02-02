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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
#ifndef TREEDEC_GRAPH_UTIL_HPP
#define TREEDEC_GRAPH_UTIL_HPP

#include <boost/graph/adjacency_list.hpp>
#include "graph.hpp"
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

template<class T>
inline
detail::visited_mask<std::vector<T> > make_incidence_mask(std::vector<T>& v)
{
	return detail::visited_mask<std::vector<T> >(v);
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
            }else{
				}
        }
    }
    return missing_edges;
} // count_missing_edges


template<class G>
class MissingEdgeCounter{
	typedef typename boost::graph_traits<G>::edges_size_type count_t;
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_t;
public:
	explicit MissingEdgeCounter(G const& g) : _g(g) { }

public:
	count_t cme(vertex_t v){
		return count_missing_edges(v, _g);
	}

private:
	G const& _g;
};

// return value: true - vertex removed
//              false - vertex isolated
template <typename G>
bool eliminate_vertex(
        const typename boost::graph_traits<G>::vertex_descriptor v, G& g)
{
#ifndef NDEBUG
	auto check = boost::num_vertices(g);
#endif

	auto adjv = boost::adjacent_vertices(v, g);


	adjv = boost::adjacent_vertices(v, g);
	make_clique(adjv, g);

	typedef typename boost::graph_traits<G>::directed_category Cat;
	if(boost::detail::is_directed(Cat())){ untested();
		for(;adjv.first!=adjv.second; ++adjv.first){ untested();
			boost::remove_edge(*adjv.first, v, g);
		}
	}else{
	}
	boost::clear_vertex(v, g);	

	boost::remove_vertex(v, g);	

	assert(check-1 == boost::num_vertices(g));
	return true;
}

} // treedec

#endif
