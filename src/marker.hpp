// Felix Salfelder, 2016-2017
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
#ifndef TREEDEC_MARKER_HPP
#define TREEDEC_MARKER_HPP

// inspired by boost/graph/md/marker
// could just include if there was a seperate header...

namespace treedec{

namespace draft{

template<class U, class key_type /* idmap... */>
class sMARKER{
	BOOST_STATIC_ASSERT( std::numeric_limits<U>::is_integer
	                 && !std::numeric_limits<U>::is_signed);
	BOOST_STATIC_ASSERT( std::numeric_limits<key_type>::is_integer
	                 && !std::numeric_limits<key_type>::is_signed);
public: // duplicate in induced_supergraph.hpp...
	template<class G>
	class predicate_remove_tagged_edges {
		typedef typename boost::graph_traits< G >::vertex_descriptor vertex_t;
		typedef typename boost::graph_traits< G >::edge_descriptor edge_t;

	public:
		predicate_remove_tagged_edges(G const& g, sMARKER const& marker)
		  : _g(&g)
		  , _marker(&marker) {
		}

		bool operator()(edge_t e) const {
			vertex_t dist = boost::target(e, *_g);
			if (_marker->is_marked(dist))
				return true;
			return false;
		}

	private:
		G const* _g;
		sMARKER const* _marker;
	};
	 template<class G>
	 predicate_remove_tagged_edges<G> make_predicate(G const& g){
		 return predicate_remove_tagged_edges<G> (g, *this);
	 }
private: // hide
	sMARKER(){ unreachable(); }
	sMARKER(const sMARKER&){ unreachable(); }
public:
	sMARKER(size_t howmany)
		: _tide(1 /*so we can unmark*/),
		  _tags(howmany)
	{
		assert(sizeof(U));

	}
	~sMARKER(){
	}
	void clear(){
		if(_tide==std::numeric_limits<U>::max()){ untested();
			reset();
		}else{
			++_tide;
		}
		for(unsigned i=0; i<_tags.size(); ++i){
			assert(!is_marked(i));
		}
	}
	void mark(key_type x){
		assert(x<_tags.size());
		_tags[x] = _tide;
	}
	void unmark(key_type x){
		_tags[x] = _tide-1;
	}
	bool is_marked(key_type x) const{
		return(_tags[x] == _tide);
	}
	bool operator()(key_type x) const{
		return(_tags[x] == _tide);
	}
	bool is_done(key_type) const {
		return false;
	}
private:
	void reset(){ untested();
		std::fill(_tags.begin(), _tags.end(), 0);
		_tide = 1;
	}

private:
	U _tide;
	std::vector<U> _tags;
};

} // draft

// multicolor marker from boost/minimum_degree
template < class SignedInteger, class Vertex, class VertexIndexMap >
class Marker {
	static_assert(std::is_signed<SignedInteger>::value, "...");
	typedef SignedInteger value_type;
	typedef typename std::vector< value_type >::size_type size_type;

	static value_type done() {
		return (std::numeric_limits< value_type >::max)() / 2;
	}

public:
	template<class G>
	class predicate_tagged_edges {
		typedef typename boost::graph_traits< G >::vertex_descriptor vertex_t;
		typedef typename boost::graph_traits< G >::edge_descriptor edge_t;

	public:
		predicate_tagged_edges(G const& _g, Marker const& _marker)
			: g(&_g), marker(_marker)
		{
		}

		bool operator()(edge_t e) const
		{
			auto t = boost::target(e, *g);
			auto s = boost::source(e, *g);
			if (marker.is_tagged(t)){
				trace2("p2", s, t);
				return true;
			}else{
				trace2("!p2", t, marker.is_tagged(t));
				return false;
			}
		}

	private:
		G const* g;
		Marker const& marker;
	}; // predicate_tagged_edges
	template<class G>
	predicate_tagged_edges<G> make_predicate(G const& g){
		return predicate_tagged_edges<G> (g, *this);
	}

public:
	Marker(size_type _num, VertexIndexMap index_map)
	  : _tag(1 - (std::numeric_limits< value_type >::max)())
	  , _multiple_tag(1 - (std::numeric_limits< value_type >::max)())
	  , _extra_tag(std::numeric_limits< value_type >::max())
	  , data(_num, -(std::numeric_limits< value_type >::max)())
	  , id(index_map) {
	}

	void mark_done(Vertex node) {
		data[get(id, node)] = done();
	}
	bool is_done(Vertex node) const {
		return data[get(id, node)] == done();
	}
	void mark_tagged(Vertex node) {
		data[get(id, node)] = _tag;
	}
	void mark(Vertex node) {
		mark_tagged(node);
	}
	void unmark(Vertex node) {
		data[get(id, node)] = 1 - (std::numeric_limits< value_type >::max)();
	}
	void mark_multiple_tagged(Vertex node) {
		data[get(id, node)] = _multiple_tag;
	}
	bool is_marked(Vertex node) const {
	  	return data[get(id, node)] >= _tag;
	}
	bool is_extra(Vertex v) {
		return data[get(id, v)] >= _extra_tag;
	}
	bool is_tagged(Vertex node) const {
	  	return is_marked(node);
	}
	bool is_not_tagged(Vertex node) const {
		return data[get(id, node)] < _tag;
	}
	bool is_multiple_tagged(Vertex node) const { itested();
		return data[get(id, node)] >= _multiple_tag;
	}
	void assert_mclear() const {
#ifndef NDEBUG
		for(unsigned i=0; i<data.size(); ++i){
//			trace4("DEBUG clear", i, is_done(i), is_marked(i), is_multiple_tagged(i));
			assert(!is_multiple_tagged(i) || is_done(i));
		}
#endif
	}
	void assert_extra() const {
		assert(_extra_tag > _tag);
	}
	void assert_clear() const {
#ifndef NDEBUG
		for(unsigned i=0; i<data.size(); ++i){
//			trace4("DEBUG clear", i, is_done(i), is_marked(i), is_multiple_tagged(i));
			assert(!is_marked(i) || is_done(i));
		}
#endif
		assert(_multiple_tag <= _tag);
//		assert(_extra_tag <= _tag);
	}
	void clear() {
		increment_tag();
	}
	void increment_tag() {
		const size_type num = data.size();
		++_tag;
//		trace2("increment", _tag, _multiple_tag);
		if (_tag >= done()) { untested();
			_tag = 1 - (std::numeric_limits< value_type >::max)();
			for (size_type i = 0; i < num; ++i){ untested();
				if (data[i] < done()){ untested();
					data[i] = -(std::numeric_limits< value_type >::max)();
				}else{ untested();
				}
			}
		}else{
		}
	}
	void set_extra_tag(value_type e) {
		_extra_tag = _tag + e;
		if (_extra_tag >= done()) { untested();
			incomplete();
		}else{
		}
	}
	void mark_extra(Vertex v) {
		data[get(id, v)] = _extra_tag;
	}
	void set_multiple_tag(value_type mdeg0) {
		const size_type num = data.size();
		_multiple_tag = _tag + mdeg0;

		if (_multiple_tag >= done()) { untested();
			incomplete();
			_tag = 1 - (std::numeric_limits< value_type >::max)();

			for (size_type i = 0; i < num; i++){ untested();
				if (data[i] < done()){ untested();
					incomplete();
					// data[i] -= max - _multiple_tag or so.
					data[i] = -(std::numeric_limits< value_type >::max)();
				}else{ untested();
				}
			}

			_multiple_tag = _tag + mdeg0;
		}else{
		}
	}
	void set_tag_as_multiple_tag() {
		if(_multiple_tag>=_tag){
		}else{
			trace2("mt", _tag, _multiple_tag);
		}
		_tag = _multiple_tag;
	}
	value_type tag() const {
		return _tag;
	}
protected:
	value_type _tag;
	value_type _multiple_tag;
	value_type _extra_tag;
	std::vector< value_type > data;
	VertexIndexMap id;
};

} // treedec

#endif
