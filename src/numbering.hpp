#pragma once

#include <boost/graph/graph_traits.hpp>
// inspired by boost::mindegree
// (but more graph-centric)

namespace treedec{

namespace draft{

template<class G_t>
class NUMBERING_1 {
public:
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::type idmap_type;
    typedef boost::graph_traits<G_t> GraphTraits;
	 typedef typename GraphTraits::vertices_size_type value_type;
	 typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	 typedef std::vector<value_type> container_type;

/* TODO: cleanup
private: // don't use.
	 NUMBERING_1(){ unreachable(); }
	 NUMBERING_1( const NUMBERING_1&){ unreachable(); }
*/
private:
	NUMBERING_1(){ unreachable(); }
	NUMBERING_1(const NUMBERING_1&){ unreachable(); }
public:
	NUMBERING_1(G_t const& g, idmap_type idm)
	    : _current(0), _data(boost::num_vertices(g)),
	      _idmap(idm)
	{
		increment();
	}
	NUMBERING_1(G_t const& g)
	    : _current(0), _data(boost::num_vertices(g)),
	      _idmap(idmap_type(g, boost::vertex_index) )
	{
		trace1("NUMBERING_1", _data.size());
		if(boost::num_vertices(g)){
		}else{
		}
		increment();
	}
#if 0 // later?
	Numbering(G_t const& g, idmap_type i)
		: _current(1), _data(boost::num_vertices(g)), _idmap(idmap_type() ) {}
#endif

	void put(vertex_descriptor v) {
		auto id=get(_idmap, v);
		trace3("NUMBERING", v, id, _current);
		if(_data[id]){
			trace3("eek, already numbered", v, id, _current);
		}else{
		}
		_data[id] = _current;
	}
	void increment(value_type i=1) {
	  	_current -= i;
	}
	bool is_numbered(vertex_descriptor v) const {
		assert(get(_idmap, v) < _data.size());
		return _data[get(_idmap, v)] != 0;
	}
	bool is_not_numbered(vertex_descriptor v) const {
		assert(v < _data.size());
		assert(get(_idmap, v) < _data.size());
		return !is_numbered(v);
	}
	bool is_before(vertex_descriptor v, vertex_descriptor w) const {
		return _data[get(_idmap, v)] > _data[get(_idmap, w)];
	}
	value_type get_position(vertex_descriptor v) const{
		trace2("NUMBERING hack", v, sizeof(value_type));
		// HACK. cleanup later. NUMBERING_2?
#ifndef NDEBUG
		auto id=get(_idmap, v);
		assert(_data[id]);
#endif
		return -1 - _data[get(_idmap, v)];
	}
private:
	value_type _current;
	container_type _data;
	idmap_type _idmap;
};

} // draft

} // treedec
