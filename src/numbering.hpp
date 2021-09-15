// Felix Salfelder 2017, 2021
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
#ifndef TREEDEC_NUMBERING_H
#define TREEDEC_NUMBERING_H

#include <boost/graph/graph_traits.hpp>
// inspired by boost::mindegree
// (but more graph-centric)

namespace treedec{

namespace draft{

template<class G_t, class T=typename boost::graph_traits<G_t>::vertices_size_type>
class NUMBERING_1 {
public:
	static_assert(std::is_unsigned<T>::value);
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::type idmap_type;
    typedef boost::graph_traits<G_t> GraphTraits;
//	 typedef typename GraphTraits::vertices_size_type value_type;
	 typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	 typedef T value_type;
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
	// resoncsturct numbering from order and sns
	template<class C, class S>
	NUMBERING_1(G_t const& /*g*/, idmap_type idm, C const& data, S const& sns)
	    : _current(0),
		   _data(data.size()),
	      _idmap(idm)
	{
		for(unsigned i=0; i<data.size(); ++i){
			if(sns[i]>0){
				_data[i] = - data[i] - 1;
			}else{
				_data[i] = sns[i] - 1;
			}
			trace3("reconst numbering", i, data[i], _data[i]);
		}
	}
	NUMBERING_1(G_t const& g, idmap_type idm)
	    : _current(0),
		   _data(boost::num_vertices(g)),
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
		}else{ untested();
		}
		increment();
	}
#if 0 // later?
	Numbering(G_t const& g, idmap_type i)
		: _current(1), _data(boost::num_vertices(g)), _idmap(idmap_type() ) {}
#endif

	void put(vertex_descriptor v) {
		auto id = get(_idmap, v);
		trace3("NUMBERING", v, id, value_type(-_current));
		if(_data[id]){
			trace3("eek, already numbered", v, id, value_type(-_current));
			assert(false);
		}else{
		}
		_data[id] = _current;
	}
	void unput(vertex_descriptor v) {
		auto id = get(_idmap, v);
		_data[id] = 0;
	}
	void put_(vertex_descriptor v) {
		auto id = get(_idmap, v);
		_data[id] = _current;
	}
	bool all_done(value_type i = 0) const {
		trace3("all_done", i, _data.size(), value_type(-_current));
		trace2("all_done", -_current, value_type(value_type(_data.size()) + _current));
//		@#@all_done  i=0  _data.size()=9  _current=65526
//		// why is value_type(_data.size()) + _current not of type value_type?
 // "num + i > max_num;"
		return i + value_type(-_current) > _data.size();

		//  return -current + i > size(); }
	}
	void init(value_type i=1) {
		// hack. should be -1 already.
	  	_current = -i;
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
	bool is_after(vertex_descriptor v, vertex_descriptor w) const {
		return is_before(w, v);
	}

	// is_before(any_numbered, any_not_numbered) -> true
	// otherwise: -5 is before -10 ...
	bool is_before(vertex_descriptor v, vertex_descriptor w) const {
		return value_type(_data[get(_idmap, v)]) > value_type(_data[get(_idmap, w)]);
	}
	size_t total() const{
		return value_type(-1) - _current;
	}
	void indistinguishable(vertex_descriptor i, vertex_descriptor j) {
		_data[get(_idmap, i)] = -(get(_idmap, j) + 1 /*offset?*/);
	}
	value_type get_position(vertex_descriptor v) const{
//		trace2("NUMBERING hack", v, sizeof(value_type));
		// HACK. cleanup later. NUMBERING_2?
		auto id = get(_idmap, v);
#ifndef NDEBUG
	//	assert(_data[id]);
#endif
		return value_type(-1) - _data[id];
	}

	template<class O, class S>
	void get_ordering(O& ord, S& supernode_by_idx) {

		typedef typename O::value_type ovt;
		static_assert(std::is_signed<ovt>::value); // BUG
		static_assert(std::is_signed<typename S::value_type>::value); // BUG

		// collect the permutation info
		auto n = _data.size();
		typedef long diff_t;
		unsigned long i;

		// stash numbering in ord
		for (i = 0; i < n; ++i) {
			trace3("bp in", i, value_type(-_data[i]), supernode_by_idx[i]);
			diff_t size = supernode_by_idx[i];
			if ( size < 0 ) { untested();
				// reachable? (second call?)
			}else if ( size == 0 ) {
				ord[i] = - ( ovt(std::numeric_limits<value_type>::max()) - _data[i] ) - 1;
				supernode_by_idx[i] = ord[i] + 1;
				trace1("bp record sup", std::numeric_limits<value_type>::max());
				trace4("bp record sup", i, ord[i], _data[i], sizeof(value_type));
				trace1("bp record sup", supernode_by_idx[i]);
			} else {
				// not a supernode. now ord > 0.
            _data[i] = value_type(- _data[i] - 1); // final value.
            ord[i] = _data[i] + 1;
			}
		}

		// fill gaps in _data.
		for (i = 0; i < n; ++i) {
			if ( ord[i] > 0 ){
				// all set.
				trace2("all set", i, ord[i]);
				continue;
			}else{
				trace2("find", i, ord[i]);
			}
			long parent = i;
			while ( ord[parent] < 0 ) {
				parent = - ord[parent] - 1;
			}

			long root = parent;
			assert(ord[root]>=0);
			value_type num = ord[root];
			_data[i] = num;
			++ord[root]; // increment number in parent.
			// so next node gets this one...

			// redirect children.
			parent = i;
			while (ord[parent] < 0) {
				ord[parent] = - root - 1;
				parent = - ord[parent] - 1;
			}
		}

		// copy back into ord
		for (i = 0; i < n; i++) {
			value_type num = _data[i];
			ord[num] = i;
			_data[i] = -num-1;
			trace2("bp data?", i, num);
		}
		for (unsigned num=0; num < n; num++) {
			trace4("bp", num, _data[num], ord[num], supernode_by_idx[ord[num]]);
		}
	} // build_permutation()
private:
	value_type _current;
public: // BUG
	container_type _data;
	idmap_type _idmap;
};

} // draft

template<class G, class O>
draft::NUMBERING_1<G, typename O::value_type> make_numbering(G const& g, O const& o)
{ untested();
	auto n = draft::NUMBERING_1<G, typename O::value_type>(g, o);
	return n;
}

} // treedec

#endif
