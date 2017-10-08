#pragma once

namespace treedec {

namespace util {

namespace impl {

template<class G, class M>
class edge_map {
public:
//	typedef typename boost::property_map<G, boost::vertex_index_t>::const_type VertexIndexMap;
	typedef typename boost::graph_traits<G>::vertex_descriptor vd;
//	typedef typename Map_type::value_type va<G>::vertex_descriptor vd;

	edge_map( G const& g, M const& m) : _m(m), _g(g) {}
//			boost::iterator_property_map < value_vertex_type *,
//			VertexIndexMap, value_vertex_type, value_vertex_type&> vd_map, const
//			G& _g): vd_map(vd_map), _g(_g) {}
 
	typedef typename M::value_type value_vertex_type;

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

}

detail::visited_mask<std::vector<BOOL> > make_incidence_mask(std::vector<BOOL>& v)
{
	return detail::visited_mask<std::vector<BOOL> >(v);
}

} //util

}
