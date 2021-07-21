// Felix Salfelder, 2017
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
//
// induced subgraph
#ifndef TREEDEC_INDUCED_SUBGRAPH_HPP
#define TREEDEC_INDUCED_SUBGRAPH_HPP

#include <boost/iterator/filter_iterator.hpp>

// graph const&+member_pred const&=induced subgraph.
//
// this is just used for vertex iteration right now.

namespace treedec{

template<class G, class M, class D>
class INDUCED_SUBGRAPH_1{
public:
	typedef void vertex_bundled; // incomplete
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<G>::vertices_size_type vertices_size_type;
	typedef typename boost::graph_traits<G>::adjacency_iterator all_adj_it;
	typedef typename boost::graph_traits<G>::vertex_iterator all_vert_it;
	typedef G wrapped_type;
#if 0 // nonsense?
	struct mask_wrapper{ // BUG, not here.
		mask_wrapper(M const & m) : _m(m) {}
		bool operator[](vertex_descriptor v) const{
			// skip masked vertices.
			return maskinvert ^ _m[v];
		}
		M const& _m;
	};
#endif
	typedef boost::filter_iterator<M, all_adj_it> adjacency_iterator;
	typedef boost::filter_iterator<M, all_vert_it> vertex_iterator;
	typedef std::pair<adjacency_iterator, adjacency_iterator> adjacency_range;
	typedef std::pair<vertex_iterator, vertex_iterator> vertex_range;
public:
	INDUCED_SUBGRAPH_1(G const& g, M const& m, D const& d)
		: _g(g), _m(m), _d(d)
	{
	}
public:
	vertices_size_type num_vertices() const{
        return boost::num_vertices(_g);
	}
	vertex_range vertices() const{
        auto p=boost::vertices(_g);

        vertex_iterator fb(_m, p.first, p.second);
        vertex_iterator fe(_m, p.second, p.second);

		  return vertex_range(fb, fe);
	}
	adjacency_range adjacent_vertices(vertex_descriptor v) const{
		auto p=boost::adjacent_vertices(v, _g);

		adjacency_iterator fb(_m, p.first, p.second);
		adjacency_iterator fe(_m, p.second, p.second);

		return adjacency_range(fb, fe);
	}
	vertices_size_type out_degree(vertex_descriptor v) const{
        return _d[v]; // bug. idmap??
	}


	// wrapped_type& operator*(){ return _g; }
	wrapped_type const& operator*() const{ return _g; }
private:
	G const& _g;
	M const _m; // for now.
	D const _d; // for now.
};

} //treedec

namespace boost{

template<class G, class M, class D>
struct graph_traits<treedec::INDUCED_SUBGRAPH_1<G, M, D> >
  : graph_traits<G>
{
};

template<class G, class M, class D, class T>
struct property_map<treedec::INDUCED_SUBGRAPH_1<G, M, D>, T> {
	typedef typename property_map<
		typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::wrapped_type, T>::const_type const_type;
	typedef typename property_map<
		typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::wrapped_type, T>::type type;
};

template<class G, class M, class D>
typename property_map<typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::wrapped_type,
                      boost::vertex_index_t>::const_type
get(vertex_index_t, treedec::INDUCED_SUBGRAPH_1<G, M, D> const& g)
{
	return get(vertex_index, *g);
}

template<class G, class M, class D>
typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::vertices_size_type
num_vertices(treedec::INDUCED_SUBGRAPH_1<G, M, D> const& g)
{
	return g.num_vertices();
}

template<class G, class M, class D>
typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::vertex_range
vertices(treedec::INDUCED_SUBGRAPH_1<G, M, D> const& g)
{
	return g.vertices();
}

template<class G, class M, class D>
typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::adjacency_range
adjacent_vertices( typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::vertex_descriptor v,
		 treedec::INDUCED_SUBGRAPH_1<G, M, D> const& g)
{
	return g.adjacent_vertices(v);
}

template<class G, class M, class D>
typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::vertices_size_type
out_degree( typename treedec::INDUCED_SUBGRAPH_1<G, M, D>::vertex_descriptor v,
		 treedec::INDUCED_SUBGRAPH_1<G, M, D> const& g)
{
	return g.out_degree(v);
}

} // boost

#endif
