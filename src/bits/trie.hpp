/*
 * Copyright (C) 2017 Felix Salfelder
 * Authors: Felix Salfelder
 *
 * This file is part of "freetdi", the free tree decomposition intiative
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 *
 * a trie based container (inspired by Hisao Tamakis tw-exact implementation)
 */
#ifndef TREEDEC_TRIE_HPP
#define TREEDEC_TRIE_HPP
/*--------------------------------------------------------------------------*/
// #include "trace.hpp"
#include "triealloc.hpp"

#undef lassert
#define lassert(x)

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace detail{
/*--------------------------------------------------------------------------*/
struct some_true {
	bool operator()(unsigned) const {return true;}
};
/*--------------------------------------------------------------------------*/
template<class V>
bool eVal(some_true const*, V const&)
{
	return true;
}
/*--------------------------------------------------------------------------*/
template<class P, class V>
bool eVal(P const* p, V const& v)
{
	return (*p)(v);
}
/*--------------------------------------------------------------------------*/
template<class V>
bool eVal(void const*, V const&)
{
	return true;
}
/*--------------------------------------------------------------------------*/
} // detail
/*--------------------------------------------------------------------------*/
// a map-like container using sequences as keys.
template<class key_type, class value_type, class Alloc=std::allocator<char> >
class TRIE{
public: // types
	typedef key_type T;
	typedef typename T::value_type char_type;
	typedef typename T::value_type V;
private:
	struct node_t { //
		value_type block;
		struct node_t *left; // hmm try unisgned?
		struct node_t *right;
		char_type v;
	};
	typedef node_t node_t;
	typedef node_t const* NODEp;
public: // iter types
	static constexpr unsigned node_size=sizeof(node_t);
public: // iter types
	struct iterNode { //
		iterNode() {}
		iterNode(NODEp a, char_type b) : cur(a), v(b){
			lassert(b<depth());
		}
		NODEp cur;
		char_type v;
		operator bool() const{return cur;}
		NODEp operator->() const{return cur;}
	};
	typedef iterNode* range_scratch_type; /// hmm
	struct end_iterator{
	};
	// TODO: scratch should probably be a template arg.
	// void for local scratch?
	class const_iterator{
	protected: // construct
	public: // bug?
		template<typename PRED>
		const_iterator(node_t const* n, TRIE const& t, PRED const* p,
				range_scratch_type scratch=NULL)
		    : _stack(scratch?scratch:t.new_range_scratch()),
		      _seek(0), _own_stack(!scratch)
		{
			if(n->left){
				lassert(n->right);
			}else{
			}
			lassert(n);
			lassert(t.depth());
			{
				push(iterNode(n, 0));
				skip(p);

			//	trace1("const_iterator done", _stack[0]);
			}

			//trace1("bgin", _seek);
			if(!empty()){
				lassert(top().v==top()->block->component.back());
			}
		}
	public: // construct/move
		const_iterator(TRIE const& t, NODEp* scratch=NULL, bool end=false)
		    : _stack(scratch?scratch:make_range_scratch(t)),
		      _seek(_stack), _own_stack(!scratch)
		{ untested();
			incomplete();
		}
#if 0
		const_iterator(const const_iterator& i)
		    : _stack(i._stack),
		      _seek(i._seek),
		      _own_stack(false)
		{
			// inefficient. don't use.
				//_seek=0;
				// lassert(_stack[0]);
		}
#endif
		const_iterator(const const_iterator&& i)
		    : _stack(i._stack),
		      _seek(i._seek),
		      _own_stack(i._own_stack)
		{
			i._own_stack=false;
		}
		~const_iterator(){
			if(_own_stack){
				delete[] _stack;
			}else{
			}
		}
	private:
		// check if this is a candidate during iteration.
		bool check_top(){
			if(empty()){
				return true;
			}
			lassert(top().v<n);

			if(!top()->block){
//				trace0("no block here");
				return false;
			}else{
				lassert(_pred(top().v));
				return true;
			}
		}
		template<class P>
		void skip(P const* p){
		  	while(!check_top()){
				iterNode t=top();
				pop();
				expand(t, p);
				// lassert(top().v<n); no, if empty
			}
		}
	public:
		bool empty(){
			return !_seek;
		}
	private: // stack
		void push(iterNode p){
			_stack[_seek]=p;
			++_seek;
		}
		void pop(){
			lassert(!empty());
			--_seek;
		}
		iterNode top() const{
			lassert(_seek);
			return(_stack[_seek-1]);
		}
//		const_iterator_& operator=(const const_iterator_&& i){ untested();
//			incomplete();
//		}
	public: // ops
		value_type operator*() const{ untested();
			/// hmm better return pair?
			lassert(_stack[_seek-1]);
			auto block=_stack[_seek-1];
			lassert(block);
			lassert(cbset::contains(block->component, back()));

			return block;
		}
		node_t const* operator->() const{

			// incomplete();
			lassert(_seek);
//			lassert(_stack[_seek-1]);

			return _stack[_seek-1].cur;
		}
		bool operator!=(const const_iterator& p) const{ untested();
			incomplete();
         //trace2("!=", _seek, p._seek);
			lassert(p._stack);
			// lassert(p._stack[0]);
			return _seek || p._seek;
		}
		bool operator!=(const end_iterator&) const{
         //trace2("!=", _seek, _stack[0]);
			return _seek;
		}
		bool operator<(const const_iterator& p) const{ untested();
			incomplete();
//			lassert(_seek); no. could be empty
//			lassert(p._seek); no. could be end...
         //trace2("!=", _seek, p._seek);
			lassert(p._stack);
			lassert(p._stack[0]);
			return _stack[_seek]->v < p._stack[p._seek]->v;
		}

		// replace node (on stack) by its successors.
		template<class P=detail::some_true>
		void expand(iterNode t, P const* p){
			auto cur=t;

			if(cur->right){
				lassert(cur->v<depth());
				if(!p || detail::eVal(p, cur->v)){
					push(iterNode(cur->right, cur->v)); // for later...
					lassert(top().v<depth());
				}else{
				}
			}
			if(!t->left){
//			}else if(t->left->v==n && t->left->block){ untested();
//				// dangling left dead end.
//				//	incomplete, actually
//				push(iterNode(t->left, t.v));
//				lassert(top().v<depth());
			}else{
				push(iterNode(t->left, t.v));
				lassert(top()->block || top()->right);
			}
		}
		const_iterator& operator++(){
			iterNode tmp=top();
			pop();

			expand(tmp, (detail::some_true*)NULL);

			if(empty()){
			}else{
				skip((detail::some_true*)NULL);
			}

			return *this;
		}
		template<class P=detail::some_true>
		const_iterator& inc(P const* p=NULL){
			lassert(!empty());

			iterNode tmp=top();
			pop();

			expand(tmp, p);

			if(empty()){
			}else{
				skip(p);
			}

			return *this;
		}
	public: // extra
		unsigned back() const{
			lassert(_seek);
			lassert(_stack[_seek-1]);
			lassert(_stack[_seek-1].cur->block);
			lassert(top()->block->component.size());

			lassert(top().v==top()->block->component.back());
			return(top().v);
		}
	public: // speed hacks (doesn't work)
		range_scratch_type scratch_tail(){
			return &_stack[_seek];
		}
	private:
		range_scratch_type _stack;
		unsigned _seek;
		mutable bool _own_stack;
	}; // const_iterator
//	typedef const_iterator_<some_true> const_iterator;
public: //iter
	const_iterator begin(range_scratch_type r=NULL) const
	{
		return const_iterator(_root, *this, (void*)NULL, r);
	}
	template<typename PRED>
	const_iterator begin(PRED const* pred, range_scratch_type r=NULL) const
	{
		return const_iterator(_root, *this, pred, r);
	}
	static end_iterator end(){
		return end_iterator();
	}
protected:
	unsigned depth() const{
		return _depth;
	}

public:
	TRIE(unsigned depth, Alloc const& alloc=Alloc())
		: _allocator(alloc),
		  _root(NULL),
		  _depth(depth),
	     _size(0)
	{
			// incomplete();
		std::cerr << "incomplete ../../src/bits/trie.hpp:336:TRIE\n";
		// _root = new_node(-1u);
	}
	TRIE(Alloc const& alloc=Alloc())
		: _allocator(alloc),
		  _root(NULL),
		  _depth(0),
	     _size(0)
	{
			incomplete();
		// _root = new_node(-1u);
	}
	~TRIE(){
// TODO: deallocate
//		delete area;
	}
	void clear(){
	   _size = 0;
		assert(_depth);
		_root = new_node(-1u);
	}
	value_type& operator[](T const& sequence);
	size_t size() const{return _size;}
private:
	void insert(value_type const& block);
	node_t* new_node();
	node_t* new_node(char_type v, node_t* left=NULL, node_t* right=NULL);
private:
public:	 // incomplete.
	Alloc const& _allocator;
	node_t* _root;
	size_t _depth;
	size_t _size;
public:
	range_scratch_type new_range_scratch(unsigned howmany=1) const{
		return new iterNode[howmany*depth()];
	}
	void delete_range_scratch(range_scratch_type* x) const{ untested();
		delete[] x;
	}
}; // TRIE
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
template<class V, class B, class Alloc>
inline typename TRIE<V, B, Alloc>::node_t*
TRIE<V, B, Alloc>::new_node()
{
  Alloc a;
  node_t* newnode = (node_t*) const_cast<Alloc&>(_allocator).allocate(sizeof(node_t));
  lassert(newnode);
  newnode->v = -1u; // hmm.
  newnode->left = NULL;
  newnode->right = NULL;
  newnode->block = NULL;

  __builtin_prefetch (newnode+1, 1);
  return newnode;
}
template<class V, class B, class Alloc>
inline typename TRIE<V, B, Alloc>::node_t*
TRIE<V, B, Alloc>::new_node(V v, node_t* left, node_t* right)
{
  Alloc a;
  node_t* newnode = new_node();
  lassert(newnode);
  newnode->v = v;
  newnode->left = left;
  newnode->right = right;

  return newnode;
}

//put block into the tree rooted at v.
// v is a neighbor
// TODO: hint?
// hmmm block->component is the key,
//      block->otherstuff is the value?
template<class V, class B, class Alloc>
inline void TRIE<V, B, Alloc>::insert(B const& block)
{
	incomplete();
	operator[](block->component) = block;
}
/*--------------------------------------------------------------------------*/
// TODO. not tried to access existing item
template<class T, class B, class Alloc>
inline B& TRIE<T, B, Alloc>::operator[](const T& seq)
{
	node_t root;
	root.v = -1u;
	root.block = NULL;
	root.right = _root;
	root.left = NULL;

	node_t* par=&root;
	node_t* p=par->right;
	lassert(p);
	for (auto const&w : seq) {
		//trace2("ins loop", w, p->v);
		while(p->v<w){
			if(p->left){
				par = p;
				p = p->left; // not in component
				lassert(p->right);
			}else{
				par = p;
				p->left = new_node(-1u, NULL, NULL);
				p = p->left;
			}
		}

		if(!p->right){
 			lassert(p->v==-1u);
			p->v = w;
			par = p;
			lassert(!p->right);
			lassert(!p->left);
			p->right = new_node(-1u, NULL, NULL);
			p = p->right;
		}else

		if(p->v == w){
			//trace2("gotit", w, p->v);
			par = p;
			lassert(p->right);
			p = p->right;
		}else{
			// somewhere in between
			assert(p);
			assert(p->right);
			assert(p->v>w);

			node_t* tn = new_node(w, p, new_node(-1u, NULL, NULL));
			assert(par);
			if(p == par->left){
				assert(par!=&root);
				par->left = tn;
			}else{
				lassert(p==par->right);
				// assert(par!=&root);
				par->right = tn;
			}

			par = tn;
			p = tn->right;
			lassert(tn->right);
		}
	}
	_root = root.right;
	assert(par);
	assert(!p->block);
	++_size;
	return p->block;
}
/*--------------------------------------------------------------------------*/
template<class TRIE>
typename TRIE::range_scratch_type make_range_scratch(
		TRIE const& t, unsigned howmany=1)
{
	return t.new_range_scratch(howmany);
}
/*--------------------------------------------------------------------------*/
// iterate over all sets in pred^-1(true)
template<class TRIE, class PRED=detail::some_true>
std::pair<typename TRIE::const_iterator, typename TRIE::end_iterator>
make_range(TRIE const& t, PRED const* p=NULL,
		typename TRIE::range_scratch_type scr=NULL)
{
	return std::make_pair(t.begin(p, scr), t.end());
}
/*--------------------------------------------------------------------------*/
#if 0 // later
template<class TRIE, class value_type>
std::pair<TRIE::iter, TRIE::iter> make_range(TRIE t, begin, end)
{ untested();
	incomplete();
}
#endif
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#endif
