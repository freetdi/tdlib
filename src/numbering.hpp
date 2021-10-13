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
#include "graph_traits.hpp"
// #include <boost/property_map/property_map.hpp> ?
// inspired by boost::mindegree
// (but more graph-centric)

namespace treedec{

namespace draft{

template<class G_t, class T=typename boost::graph_traits<G_t>::vertices_size_type>
class NUMBERING_1 {
public:
	enum ORD_FLAGS{
	  oDefault,
	  oFlatten,
	  oTree,
     oAuto
	};
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
public:
	NUMBERING_1(){ untested(); }
private:
	NUMBERING_1(const NUMBERING_1&){ unreachable(); }
public:
	// reconstruct numbering from order and sns
	template<class C, class S>
	NUMBERING_1(G_t const& /*g*/, idmap_type idm, C const& data, S const& sns)
	  : _current(0)
	  , _data(data.size())
	  , _idmap(idm)
	  , _number(data.size())
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
	  : _current(0)
	  , _data(boost::num_vertices(g))
     , _idmap(idm)
	  , _number(boost::num_vertices(g))
	{
		increment();
	}
	NUMBERING_1(G_t const& g)
	  : _current(0)
	  , _data(boost::num_vertices(g))
	  , _idmap(idmap_type(g, boost::vertex_index) )
	  , _number(boost::num_vertices(g))
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

	template<class S>
	vertex_descriptor find_parent(vertex_descriptor v, S const& sns) const{ itested();
		while(sns[v] <= 0){ untested();
			trace3("fp", v, short(_data[v]), sns[v]);
			v = - _data[v] - 1;
		}
		return v;
	}

	void resize(size_t v) {
		untested();
		incomplete();
		_data.resize(v);
	}
	long operator[](vertex_descriptor v) const {
		auto id = get(_idmap, v);
		return value_type(-1) - _data[id];
		// return _data[v];
	}

	size_t size() const {
		return _data.size();
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
		auto idx = get(_idmap, i);
		_data[idx] = -value_type(get(_idmap, j) + 1) /*offset?*/;

		assert(idx<_number.size());
		_number[idx] = _current; // keep track of introduction...
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
	void get_ordering(O& ord, S& sns, ORD_FLAGS o=oDefault) { untested();
		switch(o){
		case oDefault:
			return get_ordering(ord, sns, _ord_mode);
		case oFlatten:
			untested();
			return get_ordering_(ord, sns);
		case oTree:
			untested();
			return get_ordering2(ord, sns);
		case oAuto:
			incomplete();
		}
	}
private:
	template<class O, class S>
	void get_ordering_(O& ord, S& supernode_by_idx) {

		typedef typename O::value_type ovt;
		static_assert(std::is_signed<ovt>::value); // BUG
		static_assert(std::is_signed<typename S::value_type>::value); // BUG

		// collect the permutation info
		auto n = _data.size();
		typedef long diff_t;
		unsigned long i;

		// stash numbering in ord
		// destroy negative supernode entries. (why?)
		// forget 0 size nodes.
		for (i = 0; i < n; ++i) {
				auto pp = - ( ovt(std::numeric_limits<value_type>::max()) - _data[i] ) - 1;
				pp = (pp>0)?(-1):pp;
			trace4("bp in", i, value_type(-_data[i]), supernode_by_idx[i], pp);
			diff_t size = supernode_by_idx[i];
		///	if ( size < 0 ) { untested();
		///		// don't touch supernode? (good?)
		///		// reachable? (second call?)
		///	}else
			if ( size <= 0 ) {
				// the link to the parent is in _data...
				// (should it be in supernode?)
				trace4("bp record sup?", i, ord[i], _data[i], sizeof(value_type));
				ord[i] = - ( ovt(std::numeric_limits<value_type>::max()) - _data[i] ) - 1;
				supernode_by_idx[i] = ord[i] + 1;
				trace4("bp record sup", i, ord[i], _data[i], sizeof(value_type));
			} else {
				// not a supernode. now ord > 0.
            _data[i] = value_type(- _data[i] - 1); // final value.
            ord[i] = _data[i] + 1;
				trace4("bp...", i, ord[i], _data[i], sizeof(value_type));
			}
		}
		for (i = 0; i < n; ++i) {
			if ( ord[i] > 0 ){
				// all set.
				trace2("bp ord set", i, ord[i]);
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
			trace2("bp find", i, parent);

			long root = parent;
			assert(ord[root]>=0);
			value_type num = ord[root];
			_data[i] = num;
			++ord[root]; // increment number in parent.
			// so next node gets this one...

			// redirect children.
			parent = i;
			auto op = ord[parent];

			while (op < 0) {
				trace2("bp set", parent, -root-1);
				ord[parent] = - root - 1;
				parent = - op - 1;
				op = ord[parent];
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
	} // get_ordering()

	template<class O, class S>
	void get_ordering2(O& ord, S& supernode_by_idx) {
		auto& sns = supernode_by_idx;

		typedef typename O::value_type ovt;
		static_assert(std::is_signed<ovt>::value); // BUG
		static_assert(std::is_signed<typename S::value_type>::value); // BUG

		// collect the permutation info
		auto n = _data.size();
		typedef long diff_t;
		unsigned long i;
		trace1("get_ordering2", n);

		// collect root nodes.
		for (i = 0; i < n; ++i) {
			assert(ord[i] == 0);
			auto pp = - _data[i] - 1;
			pp = (sns[i]<=0)?pp:-1;
			trace4("bp in", i, value_type(-_data[i]), sns[i], pp);
			diff_t size = sns[i];

			// sns > 0: root node
			// sns = 0: leaf node
			// sns < 0: inner node
			if ( size >= 1 ) {
            _data[i] = value_type(- _data[i] - 1); // final position of node i
            ord[i] = _data[i] + 1; // position of next block underneath.
				// at least 1. so ord>0 indicates "allocated"
			}else if ( size == 0 ) {
				// leaf node, controlled by some other node.
				ord[i] = - ( ovt(std::numeric_limits<value_type>::max()) - _data[i] ) - 1; // parent
				assert(ord[i]<0);
				_data[i] = 1 - sns[i]; // memory size
			} else {
				assert(size!=1);
				ord[i] = - ( ovt(std::numeric_limits<value_type>::max()) - _data[i] ) - 1; // parent of i
				assert(ord[i]<0);
				_data[i] = 1 - sns[i]; // memory size
			}
		}

		for (i = 0; i < n; ++i) {
			if ( ord[i] >= 0 ){
				// all set.
				trace3("bp ord set", i, ord[i], sns[i]);
			}else{
				trace3("bp ord is parent", i, ord[i], sns[i]);
			}
		}
		// fill gaps in _data.
		for (i = 0; i < n; ++i) {
			if ( ord[i] >= 0 ){
				trace3("final _data set", i, ord[i], _data[i]);
				continue;
			}else{
				auto imm_parent = - ord[i] - 1;
				trace2("bp find", i, imm_parent);
			}

			// find root. leave a trace in ord.??
			//  but dont change ord[root].
			unsigned long chld = i;
			assert(ord[chld]<0);

			auto p = - ord[chld] - 1;
			while ( ord[p] < 0 ) { untested();
				auto nup = ord[p]; // go up.
				trace3("trace", i, p, chld);
				ord[p] = chld; // leave trace
				chld = p;
				p = - nup - 1 ;
			}

			assert(ord[p]>=0);
			trace5("bp found root and position", i, p, sns[p], ord[p], chld);
			
			// go down
			while(true){
				auto s = _data[chld];
				auto next =  ord[chld];
				assert(_data[chld] < n || long(chld)==long(i));

				_data[chld] = ord[p];

				if (chld==i){
					ord[chld] = ord[p] + 1;
					ord[p] += s;
					// assert(s==sns[i]);
					trace2("bp leaf", i, sns[i]);
					break;
				}else{ untested();
					// assert(_data[chld] == sns[chld]); // no.
					trace5("bp noleaf", chld, p, i, sns[chld], ord[p]);
					ord[chld] = ord[p] + 1;
					ord[p] += s;
				}

				p = chld;
				chld = next;
			}

		}

		// copy back into ord
		for (i = 0; i < n; i++) { itested();
			value_type num = _data[i];
			ord[num] = i;
//			trace3("bp data?", i, _data[i], num);
			_data[i] = -num-1;
		}
		for (unsigned num=0; num < n; num++) { itested();
			trace4("bp go2 out", num, _data[num], ord[num], supernode_by_idx[ord[num]]);
		}
	} // get_ordering2()

public:
// 	bool is_root(vertex_descriptor v) const{
// 		return sns...
// 	}
	void set_mode_tree(){
		_ord_mode = oTree;
	}
	bool is_mode_tree() const{
		return _ord_mode == oTree;
	}

private:
	value_type _current;
	ORD_FLAGS _ord_mode{oFlatten};
	container_type _data;
	idmap_type _idmap;
private:
	// container_type _parent; // ?
	container_type _number;
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
