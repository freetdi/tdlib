// Felix Salfelder, 2017, 2021
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
// induced supergraph and supernode iteration
//
// | Includes code from boost/graph/minimum_degree_ordering.hpp.
// | =======================================================================
// |  Copyright 1997-2001 University of Notre Dame.
// |  Authors: Lie-Quan Lee, Jeremy Siek
// |
// |  Distributed under the Boost Software License, Version 1.0. (See
// |  accompanying file LICENSE_1_0.txt or copy at
// |  http://www.boost.org/LICENSE_1_0.txt)
// | =======================================================================
//

#ifndef TREEDEC_INDUCED_SUPERGRAPH_HPP
#define TREEDEC_INDUCED_SUPERGRAPH_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/graph/vertex_and_edge_range.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/detail/numeric_traits.hpp> // for integer_traits

// #include "treedec_traits.hpp"
#include "trace.hpp"
#include "marker.hpp"
#include "treedec_traits.hpp"
#include "treedec.hpp"

#include <boost/graph/graph_utility.hpp>

namespace treedec{

// further up?
template<class Iter>
struct iter_pair_range : std::pair<Iter,Iter> {
    iter_pair_range(Iter const& x, Iter const& y)
    : std::pair<Iter,Iter>(x, y) {
	 }
    iter_pair_range(std::pair<Iter,Iter> const& x)
    : std::pair<Iter,Iter>(x) { untested();
	 }
    iter_pair_range(std::pair<Iter,Iter>&& x)
    : std::pair<Iter,Iter>(std::move(x)) {
	 }
    Iter const& begin() const {return this->first;}
    Iter const& end() const {return this->second;}
    Iter& begin() {return this->first;}
    Iter& end() {return this->second;}
};

namespace detail{

// from boost/graph/mdo.hpp
// move to bits?
template < class SignedInteger >
class Stacks {
	static_assert(std::is_signed<SignedInteger>::value);
	typedef SignedInteger value_type; // why signed?
	typedef typename std::vector< value_type >::size_type size_type;

public:
	Stacks(size_type n) : data(n) {}

	//: stack
	class stack {
		typedef typename std::vector< value_type >::iterator Iterator;
		class const_iterator{
		friend class stack;
		protected:
			const_iterator(Iterator const& d, value_type c)
			  : _data(d), _current(c) { untested();
			}
		public:
			value_type operator*() const { untested();
				return _current;
			}
			const_iterator operator++() { untested();
				_current = _data[_current];
				return *this;
			}
			bool operator!=(const_iterator const&o) const { untested();
				return(_current!=o._current);
			}
			bool operator==(const_iterator const&o) const { untested();
				return(_current==o._current);
			}
		private:
			Iterator const& _data;
			value_type _current;
		};

	public:
		stack(Iterator _data, const value_type& head)
			: data(_data), current(head) { untested();
		}

		// did not use default argument here to avoid internal compiler
		// error in g++.
		stack(Iterator _data)
			: data(_data), current(-(std::numeric_limits< value_type >::max)()) {
		}

		const_iterator begin() const { untested();
			return const_iterator(data, current);
		}
		const_iterator end() const { untested();
			return const_iterator(data, -(std::numeric_limits< value_type >::max)());
		}
		void pop() {
			assert(!empty());
			current = data[current];
		}
		void push(value_type v) {
			data[v] = current;
			current = v;
		}
		bool empty() {
			return current == -(std::numeric_limits< value_type >::max)();
		}
		value_type& top() { return current; }

	private:
		Iterator data;
		value_type current;
	};

	// To return a stack object
	stack make_stack() { return stack(data.begin()); }

protected:
	std::vector< value_type > data;
};

template<class S>
class sn_iter {
public: // types
	typedef iter_pair_range<sn_iter> sn_range;
	typedef typename S::numbering_type numbering_type;
	typedef typename S::wrapped_type G;
	typedef typename S::marker_type marker_type;
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
	typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
	typedef typename std::tuple<adjacency_iterator, adjacency_iterator, vertex_descriptor> adj_range;

public: // construct
	explicit sn_iter(vertex_descriptor v, S const& s)
	  : _s(s) {
		assert(!_s.numbering().is_numbered(v));

		auto x = boost::adjacent_vertices(v, _s.g());
		_current = x.first;
		_end = x.second;

		update_leaf();
		next();
	}
	sn_iter(adjacency_iterator v, S const& s)
		: sn_iter(std::make_pair(v, v), s){ untested();
	}

public: // ops
	vertex_descriptor operator*() const {
		return current();
	}
	sn_iter& operator++() {
		assert(_current!=_end);

		if(_s.is_numbered_(*_current)){
			assert(_leaf != _leaf_end);
			++_leaf;
		}else{
			++_current;
			update_leaf();
		}

		next();
		return *this;
	}
	bool operator==(const sn_iter&e) const {
		if(e.atend()){
			return atend();
		}else{ untested();
			incomplete();
		}
		return false;
	}
	bool operator!=(const sn_iter& e) const {
		return !operator==(e);
	}
private:
	bool atend() const{
		return _current == _end;
	}

public: // construct

	// needed for end?
	sn_iter(S const& s) : _s(s){
	  _end = _current;
	  assert(_end == _current); // who knows...
	}
public: // copy
	sn_iter(S const& s, int /*dummy*/)
	  : _s(s) { untested();
	}
//	sn_iter(sn_iter const& p)
//	  :_q(*(new scratch_type(p._q))), _qq(&_q),
//	   _g(p._g), _numbering(p._numbering) { untested();
//		incomplete();
//	}
	sn_iter(sn_iter const&& p)
	  : _current(std::move(p._current))
	  , _end(std::move(p._end))
	  , _leaf(std::move(p._leaf))
	  , _leaf_end(std::move(p._leaf_end))
	  , _s(p._s) {
	}

public: // assign
	sn_iter& operator=(const sn_iter& other) { untested();
		(void) other;
		unreachable();
		assert(0);
		assert(&_s == &other._s);
		// _current = std::move(p._current);
		// _end = std::move(p._end);
		// _leaf = std::move(p._leaf);
		// _leaf_end = std::move(p._leaf_end);
		return *this;
	}
	sn_iter& operator=(const sn_iter&& p) { untested();
		_current = std::move(p._current);
		_end = std::move(p._end);
		_leaf = std::move(p._leaf);
		_leaf_end = std::move(p._leaf_end);
		assert(&_s == &p._s);
		return *this;
	}

private:
	void update_leaf(){
		while(_current!=_end){
			if(*_current>=_s._supernode_size.size()){ untested();
				break;
			}else if(_s._supernode_size[*_current]>0){
				break;
			}else{
				++_current;
			}
		}

		if(_current==_end){
		}else if(_s.is_numbered_(*_current)){
			if(_marker){ untested();
				_marker->mark(*_current);
			}else{
			}
			auto x = boost::adjacent_vertices(*_current, _s.g());
			_leaf = x.first;
			_leaf_end = x.second;
		}else{
			assert(_s.find_parent(*_current) == *_current);
		}
	}
	vertex_descriptor current() const {
		assert(_current!=_end);

		if(!_s.is_numbered__(*_current)){
			return *_current;
		}else{
			assert(_s.find_parent(*_current) == *_current);
			assert(_leaf!=_leaf_end);
//			trace1("leaf", *_leaf);
			return *_leaf;
		}
	}

	// sn_iter::
	void next() {
		while(true){
			while(_current!=_end){
				if(_marker && _marker->is_marked(*_current)){ untested();
					++_current;
				}else if(*_current>=_s._supernode_size.size()){ untested();
					break;
				}else if(_s._supernode_size[*_current]>0){
					// filter out indistinguishables?
					break;
				}else{ untested();
					++_current;
				}
			}
			if(_current == _end){
				break;
			}

			if(_marker && _marker->is_marked(*_current)){ untested();
				++_current;
			}else if(!_s.is_numbered_(*_current)){
				// found it.
				break;
			}else if(*_current != _s.find_parent(*_current)){ untested();
				// can't happen because non-parents are numbered.
				assert(false);
			}else if(_leaf != _leaf_end ) {
				if(!_s.is_numbered__(*_leaf)){
					return;
				}else if(_s.supernode_size(*_leaf)<=0){
					++_leaf;
				}else if(*_leaf != _s.find_parent(*_leaf)){ untested();
					// somebody has lumped it together. not sure when.
					if(_s.is_before_(*_current, _s.find_parent(*_leaf))){ untested();
						return;
					}else{ untested();
						++_leaf;
					}
				}else if(!_s.is_before_(*_current, *_leaf)){
					trace2("b4", *_current, *_leaf); //  *_current=6  *_leaf=7
					return;
				}else{ untested();
					trace2("b42", *_current, *_leaf); //  *_current=6  *_leaf=7
					++_leaf;
				}
			}else{
				++_current;
				update_leaf();

#if 0 // wrong?
				if(_current==_end){ untested();
				}else if(_s.is_numbered__(*_current)){ untested();
					 assert(_s.find_parent(*_current) == *_current);
					 auto x = boost::adjacent_vertices(*_current, _s.g());
					 _leaf = x.first;
					 _leaf_end = x.second;
				}else{ untested();
				}
#endif
			}
		}
	}

public: // destroy
	~sn_iter() {
	}

private: // state
#ifndef NDEBUG
public:
#endif
	adjacency_iterator _current;
	adjacency_iterator _end;
	adjacency_iterator _leaf;
	adjacency_iterator _leaf_end;
	S const& _s;
	marker_type* _marker{nullptr}; // incomplete.
}; // sn_iter


template<class S, class N>
class bag_iter {
public: // types
	typedef iter_pair_range<bag_iter> sn_range;
	typedef N numbering_type;
	typedef typename S::wrapped_type G;
	typedef typename S::marker_type marker_type;
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
	typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
	typedef typename std::tuple<adjacency_iterator, adjacency_iterator, vertex_descriptor> adj_range;
	typedef typename std::deque<adj_range> scratch_type; // stack?

public: // construct
	explicit bag_iter(vertex_descriptor v, S const& s, N const& num, scratch_type* c, marker_type* m=nullptr)
	  : _q(c?(*c):(*(new scratch_type(0))))
	  , _qq(&_q)
	  , _s(s)
	  , _numbering(num)
	  , _marker(m) {
		assert(c);
		if(c){
			_qq = NULL;
			_q.resize(0);
		}else{ untested();
		}

		_base_vertex = v;
		push_back_range(v);
		dfs();
	}
private:
	bag_iter(adjacency_iterator v, S const& s, N const& num, scratch_type* c) :
		bag_iter(std::make_pair(v, v), s, num, c){ untested();
	}

public: // ops
	vertex_descriptor operator*() const {
		return current();
	}
	bag_iter& operator++() {
		assert(!_q.empty());
		assert(std::get<0>(front_range()) != std::get<1>(front_range()));
		++std::get<0>(front_range());
		dfs();
		return *this;
	}
	bool operator==(const bag_iter&e) const {
		if(_q.empty()){
			return e._q.empty();
		}else if(e._q.empty()){
			return false;
		}else{ untested();
			return(get<0>(_q.front())==get<1>(e._q.front()));
		}
	}
	bool operator!=(const bag_iter& e) const {
		return !operator==(e);
	}

public: // construct

	// needed for end?
	bag_iter(S const& s, N const& n, scratch_type* c=NULL)
	  : _q(c?(*c):(*(new scratch_type(0)))), _qq(&_q)
	  , _s(s)
	  , _base_vertex()
	  , _numbering(n) {
			if(c){ untested();
				_qq = NULL;
				_q.resize(0);
			}else{
			}
		}
public: // copy
//	bag_iter(S const& s, int /*dummy*/)
//	  :_q(*(new scratch_type())), _qq(&_q),
//		_s(s) { untested();
//	}
//	bag_iter(bag_iter const& p)
//	  :_q(*(new scratch_type(p._q))), _qq(&_q),
//	   _g(p._g), _numbering(p._numbering) { untested();
//		incomplete();
//	}
	bag_iter(bag_iter const&& p)
	  : _q(p._q), _qq(p._qq),
	    _s(p._s),
	    _base_vertex(p._base_vertex),
	    _numbering(p._numbering) {
		p._qq = nullptr;
	}

public: // assign
	bag_iter& operator=(const bag_iter& other) { untested();
		(void) other;
		unreachable();
		assert(0);
		assert(&_s == &other._s);
		assert(&_numbering == &other._numbering);
		return *this;
	}
	bag_iter& operator=(const bag_iter&& p) { untested();
		_q = std::move(p._q);
		assert(&_s == &p._s);
		p._qq = nullptr; // assert?
		return *this;
	}

private:
	vertex_descriptor base() const {
		return _base_vertex;
	}
	vertex_descriptor current_base() const {
		assert(!_q.empty());
		assert(get<0>(_q.front()) != get<1>(_q.front()));

		return get<2>(front_range());
	}
	vertex_descriptor current() const {
		assert(!_q.empty());
		assert(get<0>(_q.front()) != get<1>(_q.front()));

		auto vd = *get<0>(front_range());
		return vd;
	}
	void push_front_range(vertex_descriptor v) { untested();
		auto x = boost::adjacent_vertices(v, _s.g());

		if(x.first==x.second){ untested();
		}else{ untested();
			_q.push_front(std::make_tuple(x.first, x.second, v));
		}
	}
	void push_back_range(vertex_descriptor v) {
		auto x = boost::adjacent_vertices(v, _s.g());

		if(x.first==x.second){
		}else{
			_q.push_back(std::make_tuple(x.first, x.second, v));
		}
	}

	// bag_iter::
	void dfs() {
		while(true){
			if(_q.empty()){
				trace0("dfs empty...");
				break;
			}else if(get<0>(_q.front()) == get<1>(_q.front())){
				trace0("dfs pop...");
				_q.pop_front();
#if 0
			}else if(_s.supernode_size(current()) <= 0
					&& _s.find_parent(current()) == base()
			//		&& _s.find_parent(current()) == base()
					){ untested();
				trace2("push parent", current(), base());
				push_back_range(current());
				//++get<0>(_q.front());
				break;
#endif
			}else if(_marker && _marker->is_marked(current())){ untested();
				++get<0>(_q.front());
			}else if(  current_base() != current() &&
					_numbering.is_before(_s.find_parent(current()), _s.find_parent(current_base()))
					// !_numbering.is_before(_s.find_parent(current_base()), _s.find_parent(current()))
			      &&  _numbering.is_before(_s.find_parent(current()), base()) 
			      ){
				auto c = current();
//				c= _s.find_parent(c);
				trace4("skip??", c, base(), current_base(), _q.size());
				trace1("push1", c);

				//if(base()==5844){ untested();
				//std::cerr << "XXX push " << c << "\n";
				//}
				if(!_marker){
				}else if(_marker->is_marked(c)){ untested();
				}else{
					_marker->mark(c);
				}
				push_back_range(c);

				// c can't be in this bag, because it is already gone.
				++get<0>(_q.front());
			}else if(current() == current_base()){ untested();
				trace1("skip1", current());
				++get<0>(_q.front());
			}else if(current() == base()){
				trace1("skip2", current());
				++get<0>(_q.front());
			}else if( _numbering.is_before(_s.find_parent(current()), base()) ) {
				trace1("skip3", current());
				// assert(_s.supernode_size(current())>0);
				++get<0>(_q.front());
			}else{
				// got it?
				if(_marker){
//					_marker->mark(current());
				}else{
				}
				break;
			}
		}
	}

private:
	typename scratch_type::value_type const& front_range() const {
		assert(_q.begin() != _q.end());
		return _q.front();
	}
	typename scratch_type::value_type& front_range() {
		assert(_q.begin() != _q.end());
		return _q.front();
	}

public: // destroy
	~bag_iter()
	{
		if(_qq){
			assert(_qq==&_q);
			delete(_qq);
		}else{
		}
	}

private: // state
#ifndef NDEBUG
public:
#endif
	scratch_type& _q;
	mutable scratch_type* _qq; // bool?
	S const& _s;
	vertex_descriptor _base_vertex;
	numbering_type const& _numbering; // _g.numbering?
	marker_type* _marker{nullptr};
}; // bag_iter

}// detail

// move into Supergraph? private scratch...
template<class G>
typename detail::sn_iter<G>::sn_range
   make_sn_range( typename G::vertex_descriptor v,
			         G const& g)
{
	return std::make_pair(
	  detail::sn_iter<G>(v, g),
	  detail::sn_iter<G>(g));
}

// move into Supergraph? private scratch...
template<class G, class N>
typename detail::bag_iter<G, N>::sn_range
   make_bag_range( typename G::vertex_descriptor v,
			         G const& g,
			         N const& n,
                  typename detail::bag_iter<G, N>::scratch_type* s=nullptr,
						typename G::marker_type* m=nullptr)
{
	return std::make_pair(
	  detail::bag_iter<G, N>(v, g, n, s, m),
	  detail::bag_iter<G, N>(g, n));
}

template < class SignedInteger, class Vertex, class VertexIndexMap >
class degreelists_marker {
public:
	typedef SignedInteger value_type;
	typedef typename std::vector< value_type >::size_type size_type;
	degreelists_marker(size_type n, VertexIndexMap id)
	  : marks(n, 0), id(id) {
	}
	void mark_need_update(Vertex i) {
		marks[get(id, i)] = 1;
	}
	bool need_update(Vertex i) const {
		return marks[get(id, i)] == 1;
	}
	bool outmatched_or_done(Vertex i) const {
		return marks[get(id, i)] == -1;
	}
	void mark(Vertex i) { untested();
		marks[get(id, i)] = -1;
	}
	void unmark(Vertex i) {
		marks[get(id, i)] = 0;
	}

private:
	std::vector< value_type > marks;
	VertexIndexMap id;
}; // degreelist_marker


// predicateRemoveEdge1
template < class Graph, class MarkerP, class NumberD, class Stack,
  class VertexIndexMap >
class remove_and_collect {
	typedef typename boost::graph_traits< Graph >::vertex_descriptor vertex_t;
	typedef typename boost::graph_traits< Graph >::edge_descriptor edge_t;

public:
	remove_and_collect(Graph& _g, MarkerP& _marker, NumberD const& _numbering,
			Stack& n_e, VertexIndexMap id)
		: g(&_g)
		  , marker(&_marker)
		  , numbering(_numbering)
		  , neighbor_elements(&n_e)
		  , id(id)
	{
	}

	bool operator()(edge_t e)
	{
		vertex_t dist = boost::target(e, *g);
		if (marker->is_tagged(dist)){ untested();
			trace1("delete", dist);
			return true;
		}else{
		}
		marker->mark_tagged(dist);
		if (numbering.is_numbered(dist)) {
			neighbor_elements->push(get(id, dist));
			trace1("delete2", dist);
			return true;
		}else{
#ifdef COUNT
			++_cnt;
			assert( dist != boost::source(e, *g));
#endif
		}
		return false;
	}

#ifdef COUNT
	size_t cnt() const{ untested();
		return _cnt;
	}
#endif

private:
	Graph* g;
	MarkerP* marker;
	NumberD const& numbering;
	Stack* neighbor_elements;
	VertexIndexMap id;
#ifdef COUNT
	size_t _cnt{0};
#endif
};

template<class G, class N, class D>
class Supergraph {
public: // types
	typedef Supergraph S;
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
//	typedef boost::identity_property_map VertexIndexMap;
	typedef typename boost::property_map<G, boost::vertex_index_t>::type VertexIndexMap;
private: // types
	typedef typename std::vector<int>::size_type size_type;
	typedef typename boost::detail::integer_traits< size_type >::difference_type diff_t;
	typedef degreelists_marker< diff_t, vertex_descriptor, VertexIndexMap >
            DegreeListsMarker;
public:
	typedef typename boost::graph_traits<G>::vertices_size_type vertices_size_type;
	typedef typename boost::graph_traits<G>::edges_size_type edges_size_type;
//	typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
	typedef Marker< diff_t, vertex_descriptor, VertexIndexMap > marker_type;
	typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
	typedef N numbering_type;
	typedef typename std::vector<diff_t> sns_type;

	template<class B>
	class is_marked_collect{
	public:
		explicit is_marked_collect(marker_type const& m, G const& g, numbering_type const& n, B& b)
			: _mark(m), _num(n), _buf(b), _g(g) { untested();
		}
		bool operator()(edge_descriptor e){ untested();
			auto v = boost::target(e, _g);
			trace1("delete?", v);
			if(_mark.is_marked_collect(v)){ untested();
				if(_num.is_numbered(v)){ untested();
					_buf.push_back(v);
				}else{ untested();
				}
				++_howmany;
				return true;
			}else{ untested();
				return false;
			}
		}
		size_t howmany() const{return _howmany;}
	private:
		marker_type const& _mark;
		numbering_type const& _num;
		B& _buf;
		G const& _g;
		size_t _howmany{0};
	};

	class isN{
	public:
		explicit isN(N const& n, G const& g) : _n(n), _g(g) {}
		bool operator()(vertex_descriptor v) const{
			return !_n.is_numbered(v);
		}
		bool operator()(edge_descriptor e) const{ untested();
			auto v = boost::target(e, _g);
			return !_n.is_numbered(v);
		}
	private:
		N const& _n;
		G const& _g;
	};
	typedef void vertex_bundled; // incomplete
	typedef typename boost::graph_traits<G>::adjacency_iterator all_adj_it;
	typedef typename boost::graph_traits<G>::vertex_iterator all_vert_it;
	typedef G wrapped_type;
	typedef typename detail::sn_iter<S> adjacency_iterator;
	typedef typename detail::bag_iter<S, N> bag_iterator;
	typedef typename bag_iterator::scratch_type scratch_type;
	typedef iter_pair_range<adjacency_iterator> adjacency_range;
	typedef iter_pair_range<bag_iterator> bag_range;
	typedef boost::filter_iterator<isN, all_vert_it> vertex_iterator;
	typedef iter_pair_range<vertex_iterator> vertex_range;
	typedef detail::Stacks< diff_t > Workspace;

private:
//	typedef boost::iterator_property_map< vertex_descriptor*, boost::identity_property_map,
//			  vertex_descriptor, vertex_descriptor& >
//				  IndexVertexMap;
	typedef boost::identity_property_map IndexVertexMap;


	Supergraph(Supergraph const& s) = delete;
public: // construct
	template<class SNS = sns_type>
	explicit Supergraph(G& g, N& num, D& d, SNS const* sns=nullptr )
	 : _g(g),
	   _numbering(num),
		_d(d),
		//_marker(boost::num_vertices(g)),
	   _marker(boost::num_vertices(g), boost::get(boost::vertex_index, _g)),
	   _work_space(boost::num_vertices(g)),
	 //  _marker(boost::num_vertices(g)),
	   _degree_lists_marker(boost::num_vertices(g), boost::get(boost::vertex_index, _g))	{

		if(sns){ untested();
			_supernode_size = sns_type(sns->begin(), sns->end());
		}else{
			_supernode_size.resize(boost::num_vertices(g), 1);
		}

		// _d.resize(num_vertices());
		if(!_d.size()){ untested();
			// trace no degrees in sg.
		}else{
			assert(num_vertices()<=_d.size());
			for(auto v : vertices()){
				auto dv = boost::out_degree(v, _g);
				_d[v] = dv;
			}
		}
	}

public:
	vertex_range vertices() const{
		auto p = boost::vertices(_g);
		isN n(_numbering, _g);

		vertex_iterator fb(n, p.first, p.second);
		vertex_iterator fe(n, p.second, p.second);

		return vertex_range(fb, fe);
	}
	vertices_size_type num_vertices() const{
		return boost::num_vertices(_g) - _numbering.total();
	}
	edges_size_type num_edges() const{
		assert(!_numbering.total()); // this only works during initialisation...
		return boost::num_edges(_g);
	}
	//auto bag_vertices2(vertex_descriptor v) const{ untested();
	//		return boost::adjacent_vertices(v, _g);
	//}
	bag_range bag_vertices(vertex_descriptor v, marker_type* m=nullptr) const{
//		trace2("bv", v, _numbering.get_position(v));
		assert(_numbering.is_numbered(v));

		assert(_stack.empty());
		return make_bag_range(v, *this, _numbering, &_stack, m);
	}
	adjacency_range adjacent_vertices(vertex_descriptor v) const{
		assert(!_numbering.is_numbered(v));
		assert(_stack.empty());
		return make_sn_range(v, *this);
	}
	vertices_size_type out_degree(vertex_descriptor v) const{ itested();
		assert(!_numbering.is_numbered(v));
		assert(!_degree_lists_marker.need_update(v));
		assert(v < _d.size());
		return _d[v]; // bug. idmap??
	}
	vertices_size_type bagsize(vertex_descriptor v) const{ itested();
		assert(_numbering.is_numbered(v));
		assert(v < _d.size());
		return _d[v] + 1; // bug. idmap??
	}
	vertex_descriptor find_parent(vertex_descriptor v_) const{ itested();
		if(v_ < _supernode_size.size()){ itested();
			long v = v_;
		//	return _numbering.find_parent(v);
			while(_supernode_size[v] <= 0){ itested();
				trace3("parent", v, _supernode_size[v], _numbering.get_position(v));
				// v = _numbering.get_position(v);
				v = - _supernode_size[v];

				if(v==0 && _supernode_size[v]==0){
					break; // HACK
				}else{
				}
				//assert(v == - _supernode_size[v_]); // really?
				assert(v>=0);
				trace3("parent2", v, _supernode_size[v], _numbering.get_position(v));
			}
			return v;
		}else{ untested();
			// no supernodes.
			return v_;
		}
	}
	bool is_numbered__(vertex_descriptor a) const{
		// .. not numbered and not pointing to a superelt.
		return _numbering.is_numbered(a);
	}
	bool is_numbered_(vertex_descriptor a) const{
		if(a == find_parent(a)){
		}else{ untested();
			trace2("is_numbered_", a, _supernode_size[a]);
			assert(0 && "is not parent");
		}
		return _numbering.is_numbered(a);
	}
	bool is_before_(vertex_descriptor a, vertex_descriptor b) const{
		assert(a == find_parent(a));
		assert(b == find_parent(b));
		return _numbering.is_before(a, b);
	}
	bool is_before(vertex_descriptor a, vertex_descriptor b) const{
		auto aa = find_parent(a);
		auto bb = find_parent(b);
		trace4("is_before", a, aa, b, bb);
		return _numbering.is_before(aa, bb);
	}
	bool is_after(vertex_descriptor a, vertex_descriptor b) const{ untested();
		auto aa = find_parent(a);
		auto bb = find_parent(b);
		return _numbering.is_after(aa, bb);
	}

	unsigned eliminate_vertex(vertex_descriptor c) {
		// switch degree.

		_marker.clear();
		_numbering.put(c);
		_marker.mark(c);
		cleanup_mark_n(c); // marks neighbors of c?
		eliminate_n(c);

		trace2("eliminated", c, _supernode_size[c]);
		_numbering.increment(_supernode_size[c]);
		assert(_numbering.is_numbered(c));

		if(_d.size()){
			size_type deg = _d[c];
			update_degrees(c, deg);
		}else{ untested();
		}

		return boost::out_degree(c, _g);
	}

	// the version that uses the built-in marker and avoids duplicates.
	template<class MM = marker_type>
	size_t count_missing_edges(vertex_descriptor c, MM* mm=nullptr) const {

		if(mm){
		}else{
			mm = (MM*)&_marker;
		}
		auto& marker=*mm;

		trace1("sg cme", c);
		marker.clear();
		marker.set_multiple_tag(0);
		marker.assert_clear();
		marker.assert_mclear();

		auto buf = _work_space.make_stack();

		size_t degv = _supernode_size[c]-1;

		assert(!marker.is_done(c));
		marker.mark(c); // BUG
		assert(!_numbering.is_numbered(c));
#ifdef DO_TRACE_
		{ untested();
			auto pp = boost::adjacent_vertices(c, _g);
			for(; pp.first!=pp.second; ++pp.first){ untested();
				auto n = *pp.first;
				trace2("cme g neighbour", c, n);
			}
		}
#endif
		auto pp = adjacent_vertices(c);
		for(; pp.first!=pp.second; ++pp.first){
			auto n = *pp.first;
			assert(!is_numbered(n));
			if(_supernode_size[n]<=0){ untested();
				assert(marker.is_done(n));
			}else if(!marker.is_marked(n)) {
				assert(n != c);
				buf.push(n);
				marker.mark(n);
				assert(_supernode_size[n]);
				degv += _supernode_size[n];
			}else if(marker.is_done(n)) { untested();
				assert(_supernode_size[n]<=0);
			}else{ itested();
				// dup.
			}
		}

#ifdef DO_TRACE_
//		marker.set_multiple_tag(0);
		trace2("cme21 ", c, degv);
		for(auto b : buf){ untested();
			trace2("cme stack", c, b);
			//auto q2 = boost::adjacent_vertices(b, _g);
			//for(; q2.first!=q2.second; ++q2.first){ untested();
			//	trace1("cme -->", *q2.first);
			//}
		}
#endif

		// maximum missing edges at each neighbour
		size_t d = degv - (_supernode_size[c] -1);

		diff_t m = 0;
		size_t e = 0; // count nonedges

		// count edges between neighbours of c
		size_t debug_e = (_supernode_size[c]-2) * (_supernode_size[c]-1);

		while (!buf.empty()) {
			marker.set_multiple_tag(++m);
			// element absorb
			size_type nn = buf.top();

			// maximum missing edges at nn
			diff_t deg = d - _supernode_size[nn];
			size_t nnedg = _supernode_size[nn] * (_supernode_size[nn]-1);
			nnedg += 2* _supernode_size[nn] * (_supernode_size[c]-1);

			assert(nn!=c);
			auto q = adjacent_vertices(nn);
			marker.mark_multiple_tagged(nn);
			marker.mark_multiple_tagged(c);
			for(; q.first!=q.second; ++q.first){
				auto w = *q.first;
				if(_supernode_size[w]<=0){ untested();
					assert(marker.is_done(w));
				}else if(marker.is_done(w)){ untested();
					assert(_supernode_size[w]<=0);
				}else if(marker.is_multiple_tagged(w)){
		//			trace4("cme mmt", c, nn, w, deg);
				}else if(marker.is_tagged(w)){
					assert(c!=w);
					assert(nn!=w);
					marker.mark_multiple_tagged(w);
//					trace5("cme got", c, nn, w, deg, _supernode_size[w]);
					nnedg += _supernode_size[w] * _supernode_size[nn];
					deg -= _supernode_size[w]; //  * _supernode_size[nn];
				}else{
		//			trace5("cme miss", c, nn, w, deg, _supernode_size[w]);
				}
			} // NC stack item loop

//			trace5("cme subsum", c, nn, _supernode_size[nn], nnedg, deg);
			debug_e += nnedg;

			assert(_supernode_size[nn]>0);
			e += deg * _supernode_size[nn];
			buf.pop();
		} // stack

		marker.set_tag_as_multiple_tag();

		assert(!(debug_e%2));
		auto debug_me = (degv*(degv-1) - debug_e)/2;

		trace5("got me", e, c, degv, debug_e, debug_me);
		assert(!(e%2));
		e/=2;
		assert(e == debug_me);
		return e;
	}

public:
	void cleanup_mark_n(vertex_descriptor c) {
		auto element_neighbor = _work_space.make_stack();

		// Create two function objects for edge removal
		typedef typename Workspace::stack WorkStack;
		remove_and_collect<G, marker_type, numbering_type, WorkStack, VertexIndexMap>
				p(_g, _marker, _numbering, element_neighbor, get(boost::vertex_index, _g));

		// Reconstruct the adjacent node list, push element neighbor in a
		// List.
		trace1("rm", c);
		remove_out_edge_if(c, p, _g);
		// during removal element neighbors are collected.
		//
#ifdef COUNT // later
		auto cnt = p.cnt();
#endif

		while (!element_neighbor.empty()) {
			// element absorb
			size_type e_id = element_neighbor.top();
			vertex_descriptor element = get(_index_vertex_map, e_id);
			auto ii = boost::adjacent_vertices(element, _g);
			for(; ii.first!=ii.second; ++ii.first) {
				vertex_descriptor i_node = *ii.first;
				if (_marker.is_tagged(i_node)){
				}else if( !_numbering.is_numbered(i_node)) {
					_marker.mark_tagged(i_node);

					// add node to newly created clique
					boost::add_edge(c, i_node, _g);
				}else{ untested();
				}
			}
			element_neighbor.pop();
		}
	} // cleanup_n
private:
	void eliminate_n(vertex_descriptor c) {
		auto p2 = _marker.make_predicate(_g);

		auto vv = boost::adjacent_vertices(c, _g);
		for(; vv.first!=vv.second; ++vv.first) {
			vertex_descriptor v_node = *vv.first;
			if (_degree_lists_marker.need_update(v_node)){ untested();
			}else if(!_degree_lists_marker.outmatched_or_done(v_node)) {
				// degreelists.remove(v_node);
			}else{ untested();
			}
			// update out edges of v_node
			trace1("update", v_node);

			// TODO: which edges have been added?
			// c is also marked (reinsert below **?)
			boost::remove_out_edge_if(v_node, p2, _g);

			if (0 && boost::out_degree(v_node, _g) == 0) { untested();
				// indistinguishable nodes
				_supernode_size[c] += _supernode_size[v_node];
				_supernode_size[v_node] = 0;
				trace2("indistinguishable", v_node, c);
				_numbering.indistinguishable(v_node, c);
				_marker.mark_done(v_node);
				_degree_lists_marker.mark(v_node);
			} else {
			  	// not indistinguishable nodes
				//
				// what if already v_node->c **?
				boost::add_edge(v_node, c, _g);
				_degree_lists_marker.mark_need_update(v_node);
			}
		}
	} // eliminate_n()
	bool eliminate_old(vertex_descriptor c) { untested();
		assert(!_numbering.is_numbered(c));

		_marker.clear();
		auto rg = boost::adjacent_vertices(c, _g);
		auto k = 0;
		for (; rg.first!=rg.second; ++rg.first){ untested();
			auto n = *rg.first;
			assert(n!=c);
			_marker.mark(n);
			trace1("mark", n);
			++k;
		}

		auto rs = adjacent_vertices(c);
		k = 0;
		for (; rs.first!=rs.second; ++rs.first){ untested();
			auto n = *rs.first;
			assert(n!=c);
			_marker.mark(n);
			trace1("mark", n);
			++k;
		}

		trace2("marked", c, k);
		auto dc = boost::out_degree(c, _g);
		assert(k==dc);

		typedef std::vector<vertex_descriptor> B;
		B buf;

		auto r = adjacent_vertices(c);
		for (; r.first!=r.second; ++r.first){ untested();
			auto n = *r.first;
			trace2("cleanup", c, n);
			buf.clear();
			is_marked_collect<B> m(_marker, _g, _numbering, buf);
			boost::remove_out_edge_if(n, m, _g);

			if(buf.size()){ untested();
				boost::add_edge(n, c, _g);
			}else{ untested();
			}
			trace3("degs", n, _d[n], m.howmany());
			_d[n] -= m.howmany();
			for(auto b : buf){ untested();
				assert(_numbering.is_numbered(b));
//				_d[n] -= _d[b] - 1;
			}
			--_d[n]; // disconnect center
			_d[n] += dc - 1; // connect clique.
		}

		_d[c] = dc; // assert?

		_numbering.put(c);
		_numbering.increment();
		return true; // vertex gone.
	}

	// single bmd::update
	void update_degrees(vertex_descriptor current, size_type& min_degree) {
		trace1("update", current);
		assert(_numbering.is_numbered(current));
		size_type min_degree0 = min_degree + delta + 1;

		size_type deg, deg0 = 0;

		_marker.set_multiple_tag(min_degree0);
		auto q2list = _work_space.make_stack();
		auto qxlist = _work_space.make_stack();

		auto ci = boost::adjacent_vertices(current, _g);
		for (; ci.first != ci.second; ++ci.first) {
			auto i_node = *ci.first;
			trace2("update", i_node, _supernode_size[i_node]);
			const size_type i_id = get(boost::vertex_index, _g, i_node);
			if (_supernode_size[i_node] >= 0) {
				deg0 += _supernode_size[i_node];
				_marker.mark_multiple_tagged(i_node);

				if (!_degree_lists_marker.need_update(i_node)) { untested();
					// not getting here? not outmatching yet?
				}else if (boost::out_degree(i_node, _g) == 2){
					q2list.push(i_id);
				}else{
					qxlist.push(i_id);
				}
			}else{ untested();
				// not getting here because only numbered nodes have zero supernode size?
			}
		}

		// neighbors of current that need update.
		while (!q2list.empty()) {
			const size_type u_id = q2list.top();
			trace2("q2", u_id, _numbering.is_numbered(u_id));

			// neighbor of the eliminated node, or adjacent supernode
			auto u_node = get(_index_vertex_map, u_id);

			// if u_id is outmatched by others, no need to update degree
			if (_degree_lists_marker.outmatched_or_done(u_node)) { untested();
				trace1("outmatched", u_node);
				q2list.pop();
				continue;
			}else{
			}

			_marker.increment_tag();
			deg = deg0;

			auto nu = boost::adjacent_vertices(u_node, _g).first;
			auto neighbor = *nu;
			if (neighbor == u_node) { untested();
				// why would this be reachable?
				++nu;
				neighbor = *nu;
			}else{
			}

			// walk cliques? anchored at neighbor
			if (_numbering.is_numbered(neighbor)) {
				trace2("walk", u_id, neighbor);
				auto ii = boost::adjacent_vertices(neighbor, _g);
				for(; ii.first!=ii.second; ++ii.first) {
					auto i_node = *ii.first;

					if (i_node == u_node){
					}else if ( _supernode_size[i_node] <= 0){ untested();
					}else if (!_marker.is_tagged(i_node)) {
						_marker.mark_tagged(i_node);
						deg += _supernode_size[i_node];
					}else if (!_degree_lists_marker.need_update(i_node)) {
						// nothing to do.
					}else if (0 && boost::out_degree(i_node, _g) == 2) { untested();
						// is indistinguishable
					 	_supernode_size[u_node] += _supernode_size[i_node];
					 	_supernode_size[i_node] = 0;
					 	trace2("indistinguishable2", i_node, u_node);
					 	_numbering.indistinguishable(i_node, u_node);
					 	_marker.mark_done(i_node);
						_degree_lists_marker.mark(i_node);
					} else {
						// is outmatched
						trace1("outmatch?", i_node);
						// _degree_lists_marker.mark(i_node);
					}
				}
			}else{
				deg += _supernode_size[neighbor];
			}

			deg -= _supernode_size[u_node];
			_d[u_node] = deg; // update degree
			trace2("new degree2", u_node, deg);
			// degreelists[deg].push(u_node);
			// u_id has been pushed back into degreelists
			_degree_lists_marker.unmark(u_node);
			if (min_degree > deg){
				min_degree = deg;
			}else{
			}
			q2list.pop();
		} // while (!q2list.empty())

		while (!qxlist.empty()) {
			const size_type u_id = qxlist.top();
			// neighbor of the eliminated node, or adjacent supernode
			auto u_node = get(_index_vertex_map, u_id);

			// if u_id is outmatched by others, no need to update degree
			if (_degree_lists_marker.outmatched_or_done(u_node)) { untested();
				trace1("outmatched", u_node);
				qxlist.pop();
				continue;
			}else{
			}
			_marker.increment_tag();
			deg = deg0;
			auto ii = boost::adjacent_vertices(u_node, _g);
			for(; ii.first!=ii.second; ++ii.first) {
				auto i_node = *ii.first;
				if (_marker.is_tagged(i_node)){ untested();
					continue;
				}else{
				}
				_marker.mark_tagged(i_node);

				if (_numbering.is_numbered(i_node)) {
					auto jj = boost::adjacent_vertices(i_node, _g);
					for(; jj.first!=jj.second; ++jj.first) {
						auto j_node = *jj.first;
						if (_marker.is_not_tagged(j_node)) { itested();
							_marker.mark_tagged(j_node);
							deg += _supernode_size[j_node];
						}else{
						}
					}
				} else{
					deg += _supernode_size[i_node];
				}
			} // for adjacent vertices of u_node
			deg -= _supernode_size[u_node];
			trace2("new degree", u_node, deg);
			assert(u_node < _d.size());
			_d[u_node] = deg;
			// degreelists[deg].push(u_node);
			// u_id has been pushed back into degreelists
			_degree_lists_marker.unmark(u_node);
			if (min_degree > deg){
				min_degree = deg;
			}else{
			}
			qxlist.pop();
		} // while (!qxlist.empty())

		_marker.set_tag_as_multiple_tag();
//		llist.pop();

	} // update()

public:
	long supernode_size(vertex_descriptor i) const { itested();
		assert(i < _supernode_size.size());
		return _supernode_size[i];
	}
	G const& g() const {
		return _g;
	}
	N const& numbering() const{
		return _numbering;
	}
	wrapped_type const& operator*() const{
		return _g;
	}
	bool is_numbered(vertex_descriptor v) const{
		return _numbering.is_numbered(v);
	}

private:
public: // BUG
	G& _g;
private:
	N& _numbering;
	D& _d; // for now.
public: // BUG
	sns_type _supernode_size; // for now.
private:
	mutable scratch_type _stack;
public: // BUG
	mutable marker_type _marker;
	mutable Workspace _work_space; // move to greedy_base?
private:

	int delta{0};
	IndexVertexMap _index_vertex_map;
	DegreeListsMarker _degree_lists_marker; // not here.
}; // Supergraph

template<class G, class M, class D>
void eliminate_vertex( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D>& g)
{
	g.eliminate_vertex(v);
	// return true?
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertices_size_type
count_missing_edges( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g)
{
	return g.count_missing_edges(v);
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertices_size_type
count_missing_edges( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g, M&)
{ untested();
	// use built-in marker.
	// assert(&m==&g._marker)?
	return g.count_missing_edges(v);
}

template<class G, class M, class D, class MM>
typename treedec::Supergraph<G, M, D>::vertices_size_type
count_missing_edges( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		MM& m, treedec::Supergraph<G, M, D> const& g)
{ untested();
	// use built-in marker.
	// assert(&m==&g._marker)?
	return g.count_missing_edges(v, &m);
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertices_size_type
count_missing_edges( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		M&, treedec::Supergraph<G, M, D>& g)
{ untested();
	// use built-in marker.
	// assert(&m==&g._marker)?
	return g.count_missing_edges(v);
}

template<class G, class N, class D, class S>
Supergraph<G, N, D> const make_supergraph(G const& g, N const& n, D const& d, S const& s)
{ untested();
	G& g_ = const_cast<G&>(g);
	N& n_ = const_cast<N&>(n);
	D& d_ = const_cast<D&>(d);
	return Supergraph<G, N, D>(g_, n_, d_, &s);
}

namespace draft{

template<class V, class N, class O, class S, class P>
void visit_bag(V v, N const& num, O const& ordering, S const& sns, P& visitor)
{
//	auto const& num = s.numbering();
// auto s = g.supernode_size(v);
	auto s = sns[v];
//
	auto p = num.get_position(v);

	bool descend = visitor(v);
	if(s<0){
		s = -s;
	}else{
		--s;
	}

	if(descend){
		trace2("descend", v, s);
		for(auto h=1; h<=s; ){
			auto w = ordering[p+h];
			visit_bag(w, num, ordering, sns, visitor);
			auto cs = sns[w];
			assert(cs<=0);
			h -= cs;
			++h;
		}
	}else{
	}
}

template <typename G, typename O, class T, class S>
void tree_from_sg(G &s, O const& o, T& t, size_t bagsize, S const& sns)
{
//    typedef typename boost::property_map<G, boost::vertex_index_t>::type::value_type vertex_index_type;
    size_t num_vert = o.size();
    auto const& num = s.numbering();

    if(num_vert == 0){ untested();
        boost::add_vertex(t);
    }else{
        assert(num_vert == o.size());
        O iOlocal;

        std::vector<unsigned> nb; // which positions in o correspond to bags?
        std::vector<long> nbi(num_vert, -1ul);
        std::vector<long> nbe(num_vert, -1ul);

        for(unsigned i=0; i<num_vert; ++i){
            int oi = o[i];
            long w = sns[oi];

            if(w > 0){
                trace6("bag prep", w, i, oi, sns[oi], o.size(), nb.size());
                nbi[oi] = nb.size();

                for(int y=0; y<w; ++y){
						 assert(i+y<o.size());
						 nbi[o[i+y]] = nb.size();
                }

                nbe[nb.size()] = oi;
                nb.push_back(i);
            }else{
                trace3("nobag", i, oi, sns[oi]);
            }
        }

        auto nodes_left = o.size();

        trace2("loop done", nodes_left, num_vert);
        assert(! boost::num_vertices(t));

        std::vector<unsigned> edges(nb.size()-1, -1u);

			//boost::add_vertex(t);
		  t = T(nb.size());
        treedec::set_bagsize(t, bagsize);

        for(unsigned i = 0; i < num_vert; ++i){
			  s._marker.unmark(i);
		  }

        for(unsigned i = 0; i < nb.size(); ++i){
            s._marker.clear();

            int oi = o[nb[i]];
            //assert(i+1 == boost::num_vertices(t));
            auto& b = boost::get(treedec::bag_t(), t, i);

            trace3("bag contents start...", i, oi, sns[oi]);
            // TODO: don't use o;
            // TODO: include with s.bag_vertices(oi)?
#if 0
            unsigned pos = nb[i];
            auto w = sns[oi];
            while(w>0){ untested();
                trace3("... opush?", i, pos, o[pos]);
                assert(pos < o.size());
                push(b, o[pos]);
                s._marker.mark(o[pos]);
                ++pos;
                --w;
                w = 0;
            }
#else
            push(b, oi);
#endif

            trace4("==== bag contents rest...", i, oi, sns[oi], num.get_position(oi));
            trace1("==== bag contents rest...", num.is_numbered(oi));
            auto bb = s.bag_vertices(oi, &s._marker);
            for(; bb.first!=bb.second; ++bb.first){
					long v = *bb.first;
//					if(oi==5844){ untested();
//						std::cerr << "c DEBUG " << oi << " " << v << " " << num.get_position(v) << "\n";
//					}
                trace2("... bag contents", v, sns[v]);
                trace1("... bag contents", num.get_position(v));
                assert(sns[v] <= 0 || !num.is_before(v, oi));

#if 0
                if( s._marker.is_done(v)){ untested();
						 //BUG
						 s._marker.unmark(v);
						 trace1("undone marker", v);
					 }else{ untested();
					 }
#endif
                if( s._marker.is_marked(v)){
						 trace2("marked", i, v);
                    // incomplete();
                }else{
                    s._marker.mark(v);
                    auto pv = s.find_parent(v);

                    assert(!num.is_before(pv, oi));
                    // assert(num.is_before(oi, pv)); no.

                    trace4("... push", oi, v, sns[v], pv);
                    push(b, v);

                    assert(v>=0);
                    assert(pv<sns.size());
                    assert(pv<nbi.size());
                    assert(nbi[pv]<long(nb.size()));
                    assert(sns[pv]>=0);

                    if(i>=edges.size()){
                        trace2("noneed.(.", i, edges.size());
                        // no need for edge.
                    }else if( nbi[v] >= long(nb.size())){ untested();
                        trace2("not a bag...", v, nbi[v]);
                         // not a bag
//                    }else if(sns[v] <= 0){ untested();
//                        trace3("edg no cand", i, v, sns[v]);
                    }else if(edges[i] == -1u){
							  if( i!=nbi[v]) {
								  edges[i] = nbi[v];
									assert(edges[i] < nb.size());
							  }else{
							  }
                        trace4("edg init", i, v, nbi[v], nbi[pv]);
                        assert(edges[i] != i);
                    }else if(s.is_before(pv, nbe[edges[i]])){
                        trace3("edg upd", i, v, nbi[v]);
                        trace4("edg upd", nbi[pv], pv,  nbe[edges[i]], edges[i]);
                        assert(pv<nbi.size());
                        if(v>=long(nbi.size())){ untested();
                        }else if(i==nbi[v]){ untested();
                        }else if(nbi[v] < long(nb.size())){
                            edges[i] = nbi[v];
									trace2("edg again??", i, edges[i]);
                        }else{ untested();
                        }
                        assert(edges[i] != i);
                        assert(edges[i] < nb.size());
                    }else if(!s.is_before(nbe[edges[i]], pv)){
                        trace2("draw", nbe[edges[i]], pv);
                    }else{
                        trace6("miss", nb.size(), i, v, edges[i], nb[edges[i]], o[nb[edges[i]]]);
                    }
                }
            }
            if(i>=edges.size()){
				}else if(edges[i]==-1u){
					trace1("BUG", i);
					unreachable();
					// assert(false);
				}else{
                trace1("set edg", edges[i]);
                trace5("set edg", i, oi, nbe[i], edges[i], nbe[edges[i]]);
                assert(edges[i] < nb.size());
                auto v = nbe[edges[i]];
                auto w = nbe[i];

                // omitted during mass elimination.
                // tmp hack
                if(boost::edge(v, w, *s._g).second){ untested();
                }else{
                    trace2("extra edg", v, w);
                    boost::add_edge(v, w, *s._g);
                }
            }

        } // bag loop

        // invert edge direction?
        for(unsigned i = 0; i < edges.size(); ++i){
            trace3("edg", i, nb[i], edges[i]);
            boost::add_edge(i, edges[i], t);
        }

//        trace4("edg", boost::num_edges(t), boost::num_vertices(t), nb.size(), edges.size());
//        assert(boost::num_edges(t) +1 == boost::num_vertices(t));

        for(unsigned i = 0; i < num_vert-nodes_left; i++){ untested();
            auto& b=boost::get(treedec::bag_t(), t, i);
//            assert(b.size());
            treedec::sort(b); // bug in check_tree_decomp, need sorted bags...
        }
    }
} // tree_from_sg

} // draft
} //treedec

namespace boost{

template<class G, class M, class D>
void clear_vertex( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g)
{ untested();
	incomplete();
}

template<class G, class M, class D>
void remove_vertex( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g)
{ untested();
	incomplete();
}

template<class G, class M, class D>
struct graph_traits<treedec::Supergraph<G, M, D> >
  : graph_traits<G> { //
	  typedef treedec::Supergraph<G, M, D> type;
	  typedef typename type::adjacency_range adjacency_range;
};

template<class G, class M, class D, class T>
struct property_map<treedec::Supergraph<G, M, D>, T> {
	typedef typename property_map<
		typename treedec::Supergraph<G, M, D>::wrapped_type, T>::const_type const_type;
	typedef typename property_map<
		typename treedec::Supergraph<G, M, D>::wrapped_type, T>::type type;
};

template<class G, class M, class D>
typename property_map<typename treedec::Supergraph<G, M, D>::wrapped_type,
                      boost::vertex_index_t>::const_type
get(vertex_index_t, treedec::Supergraph<G, M, D> const& g)
{
	return get(vertex_index, *g);
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertices_size_type
num_vertices(treedec::Supergraph<G, M, D> const& g)
{
	return g.num_vertices();
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::edges_size_type
num_edges(treedec::Supergraph<G, M, D> const& g)
{
	return g.num_edges();
}

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertex_range
vertices(treedec::Supergraph<G, M, D> const& g)
{
	return g.vertices();
}

#if 1
template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::adjacency_range
adjacent_vertices( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g)
{
	return g.adjacent_vertices(v);
}
#else
template<class G, class M, class D, class V>
typename treedec::Supergraph<G, M, D>::adjacency_range
adjacent_vertices( V v,
		 treedec::Supergraph<G, M, D> const& g)
{ untested();
	return g.adjacent_vertices(v);
}
#endif

template<class G, class M, class D>
typename treedec::Supergraph<G, M, D>::vertices_size_type
out_degree( typename treedec::Supergraph<G, M, D>::vertex_descriptor v,
		 treedec::Supergraph<G, M, D> const& g)
{
	return g.out_degree(v);
}

} // boost

#endif
