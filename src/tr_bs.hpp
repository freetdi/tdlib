// Felix Salfelder, 2021
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

// This is derived from MIT licensed java code from the PACE submission
// Copyright (c) 2017, Hisao Tamaki

#ifndef TR_BS_H
#define TR_BS_H
#include <bit>
#include <bitset>

// #define STORE_SIZE

#define MIN_NODES 2
namespace detaiL{
/*--------------------------------------------------------------------------*/
struct some_true {
	bool operator()(unsigned) const {return true;}
};
/*--------------------------------------------------------------------------*/
template<class V>
bool eVal(some_true const*, V const&)
{ untested();
	return true;
}
/*--------------------------------------------------------------------------*/
template<class P, class V>
bool eVal(P const* p, V const& v)
{ untested();
	return (*p)(v);
}
/*--------------------------------------------------------------------------*/
template<class V>
bool eVal(void const*, V const&)
{ untested();
	return true;
}
/*--------------------------------------------------------------------------*/
} // detaiL
/*--------------------------------------------------------------------------*/


static const std::string spaces64 =
    "                                                                ";

template<class cfg_myset>
class KEYXS{
	KEYXS(KEYXS const&) = delete;
public:
	explicit KEYXS(cfg_myset const& s) : _s(s){
		// trace1("XS", _s);
	}

	size_t size() const{
		return (cfg_myset::max_element + 1)/64;
	}
	ulong const& operator[](int i) const{
//		trace3("XS", _s, i, _s.chunk(i));
		return _s.chunk(i);
	}
private:
	cfg_myset const& _s;
};

template<class key_type, class value_t,int MAX_CHILDREN_SIZE = 512>
class BlockSieve{
	typedef std::vector<uint64_t> long_array;
	typedef key_type cfg_myset;
public:
	typedef std::vector< value_t > set_out_hack;
	static size_t get_size(cfg_myset const* s)
	{ itested();
		assert(s);
		return s->size();
	}
	static size_t get_size(cfg_myset const& s)
	{ untested();
		return s.size();
	}

private:
	template<class label_t, class label_ft>
	class Node_ {
		static constexpr unsigned max_width = sizeof(label_t)*8;
		friend class BlockSieve;

	public:
//		using NodeBase::add;
		int index() const{
			return _index;
		}

		ulong getMask() const{ itested();
			assert(_width);
			// 0000000000011111111111111100000000000
			//           ^ntz+width     ^ntz       ^lsb=0
			ulong a = ulong(-1) << _ntz;
			ulong b = ulong(-1) >> (64 - _width - _ntz);
			return a & b;
		}

		bool isLeaf() const {
			return (_index == 0) && _ntz == 0;
		}

		bool isFirstInInterval() const {
			return _ntz == 0; // _ntz + _width == 64;
		}

		bool isLastInInterval() const{ untested();
			return _ntz + _width == 64;
		}
	public:
		explicit Node_() {
			_data._values = nullptr;
			_data._children = nullptr;
		}
		explicit Node_(int index, int width, int ntz)
			: _index(index), _width(width), _ntz(ntz) {

			_data._values = nullptr;
			_data._children = nullptr;
			assert(width <= int(8*sizeof(label_t)));
		}
		Node_(const Node_&p)
			: _index(p._index), _width(p._width), _ntz(p._ntz) { untested();
			assert(!p.size());
			assert(!p._data._values);
			_data._values = nullptr;
			_data._children = nullptr;
		}
		Node_(Node_&&p) noexcept
			: _index(p._index), _width(p._width), _ntz(p._ntz), _size(p.size),
			_labels(p._labels),
			_data(p._data)
		{ untested();

			_size = p._size;
			p._data = nullptr;
			p._labels = nullptr;
			p._size = 0;
		}
		//Node_(const Node_&&p) :
		//	NodeBase(p),
		//	_children(p._children)
		//{ untested();

		//}
		Node_& operator=(Node_& x) = delete;
		Node_& operator=(Node_&& p) noexcept {
			_index = p._index;
			_width = p._width;
			_ntz = p._ntz;
#ifdef STORE_SIZE
			_cardinalities = std::move(p._cardinalities);
#endif
			_labels = p._labels;
			_size = p._size;
				assert(size() <= MAX_CHILDREN_SIZE);

			if(isLeaf()){
				_data._values = p._data._values;
				p._data._values = nullptr;
			}else{
				_data._children = p._data._children;
				p._data._children = nullptr;
			}
			p._labels = nullptr;
			p._size = 0;

			return *this;
		}
		~Node_(){
			if(isLeaf()){
				if(_data._values){
					free( _data._values );
				}else{
				}
			}else{
				if(_data._children){
					//delete[] _data._children;
					free(_data._children);
				}else{
				}
			}
			free(_labels);
		}

		int cardinality(unsigned i) const{
			assert(i<size());
#ifdef STORE_SIZE
			assert(i<int(_cardinalities.size()));
			assert(get_size(value(i)) == size_t(_cardinalities[i]));
			return _cardinalities[i];
#else
			return(get_size(value(i)));
#endif
		}

		size_t size() const {
			return _size;
		}
	private:
		value_t& add_value(unsigned i) {
//			assert(_children.size() == 0);
#ifdef STORE_SIZE
			assert(_values.size() == _cardinalities.size());
#endif
        for(unsigned j=size() - 1; j > i; --j){
				value(j) = value(j - 1);
			}

			//		  trace2("insert", i, value);
			//        _values[i] = new cfg_myset(value); // eek.
//			_values[i] = value; // eek.

#ifdef STORE_SIZE
			_cardinalities.resize(_cardinalities.size() + 1);
			for(int j = _cardinalities.size() - 1; j - 1 >= i; j--){ untested();
				_cardinalities[j] = _cardinalities[j - 1];
			}

			assert(i<int(_cardinalities.size()));
			_cardinalities[i] = cardinality;
#endif
			// assert(_values[i]);
//			assert(!_children.size());
			return value(i);
		}
		// really??
//		value_type& add(ulong label, value_type const& value){ untested();
//			return add(label, value, get_size(value));
//		}

		void add_child(Node_* child, unsigned i) {
			assert(!isLeaf());
			assert(_data._children);

			for(unsigned j = size()-1; j > i; j--){
				_data._children[j] = std::move(_data._children[j - 1]);
			}
			_data._children[i] = std::move(*child);
		}
		int add(ulong label, Node_* child) {
			assert(!isLeaf());

			auto p = add_label(label);
			int i = p.second;
			node_t* n = p.first;
			//		 trace2("add", _children.size(), i);
			
			assert(child->size() <= MAX_CHILDREN_SIZE);
			n->add_child(child, i);

			return i;
		}

		Node_ const& child(unsigned i) const{
			assert(!isLeaf());
			assert(i<size());
			return _data._children[i];
		}
		value_t const& value(unsigned i) const{
			assert(isLeaf());
			assert(i<size());
			return _data._values[i];
		}
		value_t& value(unsigned i){
			assert(isLeaf());
			assert(i<size());
			assert(i<MAX_CHILDREN_SIZE);
			assert(_data._values);
			return _data._values[i];
		}
		int find_mid(int& p) const{
			int lo = 1;
			int hi = _width-1;
			std::vector<ulong> l(size()); //  = new ulong[leng];

			int mid = (lo+hi)/2;
			while(lo != mid) {
				trace3("====try", lo, mid, hi);
				assert(lo<=mid);
				assert(hi>=mid);
				int leftntz = _ntz + mid;
				// ulong m = ulong(-1) >> ( 64 - hi );

				p = 0;
				for(size_t i = 0; i < size(); i++){
					ulong l = label(i);
					trace2("full", std::bitset<64>(l), _ntz);
				}

				ulong prev = label(0);
				prev = prev >> leftntz;
				int cnt=1;
				assert(size());
				for(size_t i = 1; i < size(); i++){
					ulong l = label(i) >> leftntz;
					if(l==prev){
					}else{
						prev = l;
						++cnt;
					}
				}

				p = cnt;

				trace5("next?", p, lo, mid, hi, size());
				if(2*p==int(size())){ untested();
					break;
				}else if(2*p<int(size())){
					trace5("too few", p, lo, mid, hi, size());
					assert(hi!=mid);
					// need more high bits -> lower mid.
					hi = mid;
				}else{ untested();
					trace5("too many", p, lo, mid, hi, size());
					assert(lo!=mid);
					lo = mid+1;
				}

				// similar to old condition
				if(p < 2){
					// more space.
			//	}else if(2*p > MAX_CHILDREN_SIZE){ untested();
			//		 break;
				}else if(p < MAX_CHILDREN_SIZE){
					 break;
				}else{ untested();
				}
				assert(mid != (lo+hi)/2);

				mid = (lo+hi)/2;
				trace4("next?", p, lo, mid, hi);
			}
			trace3("found", mid, p, size());

			return mid;
		}
		value_t& put(ulong bits, KEYXS<cfg_myset> const& longs, int i)  {
			node_t* node = this;

			value_t* ret = nullptr;
			size_t n = size();
			if(n == MAX_CHILDREN_SIZE){ untested();
			}else{
			}
			assert(size() + 1 <= MAX_CHILDREN_SIZE);

			if(isLeaf()){
				ret = &add(bits);
				assert(size() <= MAX_CHILDREN_SIZE);
			}else if(isFirstInInterval()){
				auto p = newPath(i - 1, longs);
				ret = p.second;
				add(bits, p.first);
				assert(size() <= MAX_CHILDREN_SIZE);
			}else{
				trace2("put mask", i, std::bitset<64>(getMask()));
				trace2("put mask     ", i, std::bitset<64>(bits));

				assert(_ntz);
				node_t* header = newNode( i, _ntz, 0);
				if(!header->isLeaf()){ itested(); // 10x10
					auto p = newPath(i - 1, longs);
					ret = p.second;
					header->add(bits, p.first);
				} else{
					ret = &header->add(bits);
				}

				assert(!isLeaf());
				add(bits, header);
			}
			assert( n+1 == node->size()); // ??

			assert(ret);
			return *ret;
		}
		ulong label(int i) const{
			assert(size_t(i)<size());
			return ulong(_labels[i]) << _ntz;
		}

	private:
		Node_& child(unsigned i){
			assert(!isLeaf());
			assert(i<size());
			assert(i<MAX_CHILDREN_SIZE);
			assert(_data._children);
			return _data._children[i];
		}

	public:
      int indexOf(ulong label){
			ulong mask = ulong(-1) >> ( 64 - _width );
			auto ml = (label >> _ntz) & mask;
			assert(std::popcount(ml)<=int(sizeof(label_t)*8));

			auto i = std::lower_bound(_labels, _labels + _size, ml);

			int ret;

			if(i == _labels + _size){
				ret = - size() - 1;
			}else if(*i == ml){
				ret = i - _labels;
			}else{
				ret = - (i - _labels + 1);
				auto insert_pos = - ret - 1;
				assert(ml<_labels[insert_pos]);
			}

			return ret;
      }

	public:
		value_t& add(ulong label) {
			assert(isLeaf());

			auto p = add_label(label);
			int i = p.second;
			node_t* n = p.first;

			auto& ret = n->add_value(i);
			return ret;
		}
		std::pair<Node_*, int> add_label(ulong l){

			assert(_width <= int(8*sizeof(label_t)));
			ulong mask = ulong(-1) >> ( 64 - _width );
			auto ml = (l >> _ntz) & mask;

			// make sure it is not there yet.
			for(unsigned i=0; i<size(); ++i){
				assert(_labels[i]!=ml);
			}

			int m = indexOf(l);
			assert(m<0); // label not there yet.
			unsigned i = - m - 1;

			assert(i==size() || ml<_labels[i]);

			//_labels.resize(_labels.size() + 1); //  = Arrays.copyOf(_labels, _labels.length + 1);
			unsigned old_cap = alloc_cap(_size);
			unsigned cap = alloc_cap(_size+1);
			if(!isLeaf()){
				trace1("noleaf", cap);
			}
			alloc_all(cap, old_cap);
			assert(_labels);
			++_size;
			trace3("size", cap, old_cap, size());
			assert(_size <= cap);

			for(unsigned j=size() - 1; j >= i + 1; j--){
				_labels[j] = _labels[j-1];
			}
			// _labels[i] = (label & getMask()) >> ntz;
			//
			assert(i<size());
			assert(_labels);
			_labels[i] = (l >> _ntz) & mask;
			// _labels[i] = (label >> ntz) & (getMask() >> ntz);
#ifndef NDEBUG
			//     for(int ii = 1; ii<_labels.size(); ++ii){ untested();
			//		  assert(_labels[ii-1] < _labels[ii]);
			//	  }
#endif

			return std::make_pair(this, i);
		}

		unsigned alloc_cap(unsigned s) const {
			unsigned c;
			if(s==0){
				c = 0;
			}else if(s<MIN_NODES){
				c = MIN_NODES;
			}else{
			 c = std::bit_ceil(s);
			}
			trace2("alloc_cap" , s, c);
			return c;
		}
		void alloc_all(unsigned cap, unsigned old_cap){
			alloc(cap, old_cap, _labels);
			assert(_labels);
			if(isLeaf()){
				trace2("leaf?", cap, old_cap);
				alloc(cap, old_cap, _data._values);
			}else{
				trace2("chld?", cap, old_cap);
				alloc(cap, old_cap, _data._children);
			}
		}
	private:

		template<class T>
		void alloc(unsigned new_cap, unsigned old_cap, T*& data){
			trace4("alloc?", new_cap, old_cap, sizeof(T), sizeof(T*));
			{
			  // new_cap = 1, 3, 5, 9, 17 ..

			  if(old_cap==0){
				  assert(!data);
					data = (T*)calloc(sizeof(T), MIN_NODES);
// 				  memset((void*)_data._children, 0, sizeof(Node_*)*2);
               for(int i=0; i<MIN_NODES; ++i){
						new (data+i) T();
					}
			  }else if(new_cap<=MIN_NODES){
			  }else{
				  T* old_data = data;
#if 0
				  data = (T*)malloc(sizeof(T) * new_cap);
				  unsigned i=0;
				  for(; i<new_cap; ++i){
					  new (data+i) T();
				  }
				  for(i=0; i<old_cap; ++i){
					  data[i] = std::move(old_data[i]);
				  }
				  free(old_data);
#else // seems faster. works for POD
				  data = (T*)realloc((void*)data, sizeof(T) * new_cap);
				  unsigned i=old_cap;
				  memset((void*)(data + old_cap),  0, sizeof(T)*(new_cap-old_cap));
#endif
			  }
		  }
		}


	private:
		std::ostream& dump(std::ostream& ps, std::string const& indent, int last) const{ untested();
        for(unsigned i=0; i<size(); ++i){ untested();
          ps << indent << _labels[i];
          if(!isLeaf()){ untested();
//				 assert(_children[i]);
            _data._children[i].dump(ps, indent + spaces64, last);
          }else{ untested();
			 }
        }
		  return ps;
      }

		// Node_::
      void filterSuperblocks( int tbs, int margin,
				KEYXS<cfg_myset> const& longs,
				KEYXS<cfg_myset> const& neighbors,
				int intersects, set_out_hack& list) const{
			assert(size()<=MAX_CHILDREN_SIZE);
//        long mask = getMask();
		  ulong mask_ = ulong(-1) >> ( 64 - _width );
        label_ft bits = 0;
        if(_index < longs.size()){
          bits = label_t((longs[_index] >> _ntz) & mask_);
        }else{ untested();
			  assert(false);
		  }

        label_ft neighb = 0;
        if(_index < neighbors.size()){
          // neighb = ((neighbors[_index] & mask) >> _ntz);
          neighb = (neighbors[_index] >> _ntz) & mask_;
        }else{ untested();
			  assert(false);
		  }

        if(isLeaf()){
          for(int i = size() - 1; i >= 0; --i){
            label_ft label = _labels[i];
            if(bits > label){
              break;
            }else if((bits & ~label) == 0){
              int intersects1 = intersects + std::popcount(label_t(label & neighb));
              if(intersects1 + cardinality(i) <= tbs){
					  list.push_back(value(i)); // HACK
              }else{
				  }
            }
          }
        } else{
          for(int i = size() - 1; i >= 0; i--){
            label_ft label = _labels[i];
            if(bits > label){
              break;
            }else if((bits & ~label) == 0){
              assert(label_t(label & neighb) == label_ft(label & neighb));
              int intersects1 = intersects + std::popcount(label_t(label & neighb));
              if(intersects1 <= margin){
                child(i).filterSuperblocks(tbs, margin, longs, neighbors, intersects1, list);
              }else{
				  }
            }else{
				}
          }
        }
      }

		// Node_::
		void filterSubblocks( int tbs, int margin,
				KEYXS<cfg_myset> const& longs,
				KEYXS<cfg_myset> const& neighbors,
				int intersects, set_out_hack& list) const { untested();

			assert(size()<=MAX_CHILDREN_SIZE);
			//        long mask = getMask();
			ulong mask_ = ulong(-1) >> ( 64 - _width );

			label_ft bits = 0;
			if(_index < longs.size()){ untested();
				bits = (longs[_index] >> _ntz ) & mask_;
			}else{ untested();
			}

			label_ft neighb = 0;
			if(_index < neighbors.size()){ untested();
				neighb = (neighbors[_index] >> _ntz) & mask_;
			}else{ untested();
			}

			if(isLeaf()){ untested();
				for(unsigned i=0; i < size(); ++i){ untested();
					label_ft label = _labels[i];
					if(bits < label){ untested();
						break;
					}else if((~bits & label) == 0){ untested();
						int intersects1 = intersects + std::popcount(label & neighb);
						if(intersects1 + cardinality(i) <= tbs){ untested();
							list.push_back(value(i));
						}else{ untested();
						}
					}else{ untested();
					}
				}
			}else{ untested();
				for(unsigned i=0; i < size(); ++i){ untested();
					label_ft label = _labels[i];
					if(bits < label) { untested();
						break;
					}else{ untested();
					}
					if((~bits & label) == 0){ untested();
						int intersects1 = intersects + std::popcount(label & neighb);
						if(intersects1 <= margin){ untested();
							child(i).filterSubblocks( tbs, margin, longs, neighbors, intersects1, list);
						}else{ untested();
						}
					}else{ untested();
					}
				}
			}
		}

	private:
		void resizeWidth(Node_* node);
		uint16_t _index{0};
		uint16_t _width{0};
		uint16_t _ntz{0};
		uint16_t _size{0};
	public: // HACK
		typedef Node_<uint64_t, uint_fast64_t> Node64;
		label_t* _labels{nullptr};

		union { //
		  Node64* _children;
		  value_t* _values;
		} _data{nullptr};
#ifdef STORE_SIZE
		std::vector<int> _cardinalities;
#endif
	}; // Node_

	typedef Node_<uint64_t, uint_fast64_t> node_t;
	typedef typename key_type::value_type char_type;

	typedef node_t const* NODEp;

public: // iter types
	struct iterNode {
		iterNode() {}
		iterNode(NODEp a, char_type b) : cur(a), v(b){ untested();
//			lassert(b<depth());
		}
		NODEp cur;
		char_type v;
		operator bool() const{return cur;}
		NODEp operator->() const{return cur;}
	};
	typedef iterNode* range_scratch_type; /// hmm
	struct end_iterator{
	};
	range_scratch_type new_range_scratch(unsigned howmany=1) const{ untested();
		return new iterNode[howmany*depth()];
	}
	unsigned depth() const{ untested();
		incomplete();
		return 128;
	}
	// TODO: scratch should probably be a template arg.
	// void for local scratch?
	class const_iterator{
	protected: // construct
	public: // bug?
		template<typename PRED>
		const_iterator(node_t const* n, BlockSieve const& t, PRED const* p,
				range_scratch_type scratch=NULL)
		    : _stack(scratch?scratch:t.new_range_scratch()),
		      _seek(0), _own_stack(!scratch)
		{ untested();
			lassert(t.depth());
			{ untested();
				push(iterNode(n, 0));
				skip(p);

			//	trace1("const_iterator done", _stack[0]);
			}

			//trace1("bgin", _seek);
		//	if(!empty()){ untested();
		//		lassert(top().v==top()->block->component.back());
		//	}
		}
	public: // construct/move
		const_iterator(BlockSieve const& t, NODEp* scratch=NULL, bool end=false)
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
		{ untested();
			// inefficient. don't use.
				//_seek=0;
				// lassert(_stack[0]);
		}
#endif
		const_iterator(const const_iterator&& i)
		    : _stack(i._stack),
		      _seek(i._seek),
		      _own_stack(i._own_stack)
		{ untested();
			i._own_stack=false;
		}
		~const_iterator(){ untested();
			if(_own_stack){ untested();
				delete[] _stack;
			}else{ untested();
			}
		}
	private:
		// check if this is a candidate during iteration.
		bool check_top(){ untested();
#if 0 // incomplete
			if(empty()){ untested();
				return true;
			}else{ untested();
			}

			if(!top()->block){ untested();
//				trace0("no block here");
				return false;
			}else{ untested();
				lassert(_pred(top().v));
				return true;
			}
#endif
			return false;
		}
		template<class P>
		void skip(P const* p){ untested();
		  	while(!check_top()){ untested();
				iterNode t=top();
				pop();
				expand(t, p);
				// lassert(top().v<n); no, if empty
			}
		}
	public:
		bool empty(){ untested();
			return !_seek;
		}
	private: // stack
		void push(iterNode p){ untested();
			_stack[_seek]=p;
			++_seek;
		}
		void pop(){ untested();
			lassert(!empty());
			--_seek;
		}
		iterNode top() const{ untested();
			lassert(_seek);
			return(_stack[_seek-1]);
		}
//		const_iterator_& operator=(const const_iterator_&& i){ untested();
//			incomplete();
//		}
	public: // ops
		value_t const& operator*() const{ untested();
			/// hmm better return pair?
			lassert(_stack[_seek-1]);
			iterNode const& block = _stack[_seek-1];
			lassert(block);

			assert(block.cur);
			return *new cfg_myset(); // block;
		}
		node_t const* operator->() const{ untested();

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
		bool operator!=(const end_iterator&) const{ untested();
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
		template<class P=detaiL::some_true>
		void expand(iterNode t, P const* p){ untested();
			incomplete();
			(void)p;
			(void)t;
#if 0
			auto cur=t;

			if(cur->right){ untested();
				lassert(cur->v<depth());
				if(!p || detaiL::eVal(p, cur->v)){ untested();
					push(iterNode(cur->right, cur->v)); // for later...
					lassert(top().v<depth());
				}else{ untested();
				}
			}
			if(!t->left){ untested();
//			}else if(t->left->v==n && t->left->block){ untested();
//				// dangling left dead end.
//				//	incomplete, actually
//				push(iterNode(t->left, t.v));
//				lassert(top().v<depth());
			}else{ untested();
				push(iterNode(t->left, t.v));
				lassert(top()->block || top()->right);
			}
#endif
		}
		const_iterator& operator++(){ untested();
			iterNode tmp=top();
			pop();

			expand(tmp, (detaiL::some_true*)NULL);

			if(empty()){ untested();
			}else{ untested();
				skip((detaiL::some_true*)NULL);
			}

			return *this;
		}
		template<class P=detaiL::some_true>
		const_iterator& inc(P const* p=NULL){ untested();
			lassert(!empty());

			iterNode tmp=top();
			pop();

			expand(tmp, p);

			if(empty()){ untested();
			}else{ untested();
				skip(p);
			}

			return *this;
		}
	public: // extra
		unsigned back() const{ untested();
			lassert(_seek);
			lassert(_stack[_seek-1]);
			lassert(_stack[_seek-1].cur->block);
			lassert(top()->block->component.size());

			lassert(top().v==top()->block->component.back());
			return(top().v);
		}
	public: // speed hacks (doesn't work)
		range_scratch_type scratch_tail(){ untested();
			return &_stack[_seek];
		}
	private:
		range_scratch_type _stack;
		unsigned _seek;
		mutable bool _own_stack;
	}; // const_iterator

public:
	explicit BlockSieve();

	~BlockSieve() {
//		delete _root;
	}

	void set_n(int n){
		// incomplete(); // clear?
		_last = (n - 1) / 64; // KEYXS<cfg_myset>...
		assert(!_root.size());
		_root._index = _last;
	}
public: //iter
	const_iterator begin(range_scratch_type r=NULL) const { untested();
		return const_iterator(&_root, *this, (void*)NULL, r);
	}
	template<typename PRED>
	const_iterator begin(PRED const* pred, range_scratch_type r=NULL) const
	{ untested();
		return const_iterator(_root, *this, pred, r);
	}
	static end_iterator end(){ untested();
		return end_iterator();
	}

public:
	value_t& operator[](key_type const& bs);

private:
  static node_t* newNode(int index, int width, int ntz){
	  assert(width);
      return new node_t(index, width, ntz);
#if 0
    if(width > 32){ untested();
      return new Node64(index, width, ntz);
    } else if(width > 16){ untested();
      return new Node32(index, width, ntz);
    } else if(width > 8){ untested();
      return new Node16(index, width, ntz);
    } else{ untested();
      return new Node8(index, width, ntz);
    }
#endif
  }

  static std::pair<node_t*, value_t*> newPath(int index, KEYXS<cfg_myset> const& longs){
    auto node = new node_t(index, node_t::max_width, 0);

    ulong bits = 0;
    if(index < 0){ untested();
		 assert(false);
	 }else if(index < int(longs.size())){
      bits = longs[index];
    }else{ untested();
	 }

	 value_t* ref = nullptr;
    if(index == 0){
      ref = &node->add(bits);
    }else{ untested();
		auto p = newPath(index - 1, longs);
      node->add(bits, p.first);
		ref = p.second;
    }

	 assert(node->size()==1);
    return std::make_pair(node, ref);
  }

public:
  void collectSuperblocks( int tbs, int margin,
		  cfg_myset const& component, cfg_myset const& neighbors,
		  set_out_hack& list) const{
		KEYXS<cfg_myset> c(component);
		KEYXS<cfg_myset> n(neighbors);
    _root.filterSuperblocks(tbs, margin, c, n, 0, list);
  }


#if 0  // not used
  public void collectSubblocks(
      XBitSet component, XBitSet neighbors, ArrayList< XBitSet > list){ untested();
    _root.filterSubblocks(component.toLongArray(), 
        neighbors.toLongArray(), 0, list);
  }
#endif

  size_t size() const{
    return _size;
  }
		 void filterSuperblocks( int tbs, int margin,
				KEYXS<cfg_myset> const& longs,
				KEYXS<cfg_myset> const& neighbors,
				int intersects, set_out_hack& list) const;

		 void filterSubblocks( int tbs, int margin,
				KEYXS<cfg_myset> const& longs,
				KEYXS<cfg_myset> const& neighbors,
				int intersects, set_out_hack& list) const;
private:

  std::ostream& dump(std::ostream& ps){ untested();
    return _root.dump(ps, "");
  }

private:
   int _n{0};
   int _last{0};
   int _tbs{0};
   int _size{0};

	node_t _root;
}; // BlockSieve

template<class key_type, class value_t, int MAX_CHILDREN_SIZE>
BlockSieve<key_type, value_t ,MAX_CHILDREN_SIZE>::BlockSieve()
  : _root(0, node_t::max_width, 0)
{
	_last = 0;
	trace1("BlockSieve", sizeof(node_t));
}

template<class key_type, class value_t, int MAX_CHILDREN_SIZE>
inline value_t&
		 BlockSieve<key_type, value_t,MAX_CHILDREN_SIZE>::operator[](key_type const& bs)
{
	KEYXS longs(bs);
	node_t* current = &_root;

	int i = _last;
	ulong bits = 0;
	for(;;){
		bits = 0;
		assert(i < int(longs.size()));
		if(i < int(longs.size())){
			bits = longs[i];
		}else{ untested();
		}
		int j = current->indexOf(bits);
		assert(current->size() <= MAX_CHILDREN_SIZE);

		if(j < 0){
//			trace2("not there", current->size(), current->getMask());
			// bits are not there...
			if(current->size() + 1 > MAX_CHILDREN_SIZE){
				current->resizeWidth(current); // (, MAX_CHILDREN_SIZE - 1);
				assert(i == current->index()); // ?

				assert(current->size());

			}else{
				break;
			}
		}else if(current->isLeaf()){ untested();
			trace2("leaf", current->size(), current->getMask());
			return current->value(j);
		}else{
			// descend
			// parent = node;

			current = &current->child(j);
			//trace1("descend", std::bitset<64>(current->getMask()));
			//trace2("descend", i, current->isLeaf());
			//trace2("descend", i, current->size());
			//trace2("descend", i, current->_children.size());
			//trace2("descend", i, current->_values.size());
			//				node = node->_children[j];
			i = current->index();
		}
	}

	value_t* ret = nullptr;
//	trace2("nput", i, std::bitset<64>(bits));
	ret = &current->put(bits, longs, i);
	assert(current->size() <= MAX_CHILDREN_SIZE);

	++_size;

	assert(current->size()<=MAX_CHILDREN_SIZE);
	assert(ret);
	return *ret;
}

template<class key_type, class value_t, int MAX_CHILDREN_SIZE>
	template<class label_t, class label_ft>
inline void // BlockSieve<key_type, value_t, MAX_CHILDREN_SIZE>::NodeBase*
       BlockSieve<key_type, value_t , MAX_CHILDREN_SIZE>::Node_<label_t, label_ft>::resizeWidth(
       BlockSieve<key_type, value_t, MAX_CHILDREN_SIZE>::Node_<label_t, label_ft>* node)
{
	assert(node==this);
	assert(node->size());
	int w = node->_width;
	ulong oldmask = getMask();
	assert(oldmask);
	int old_size = size();
	assert(old_size);
	bool was_leaf = node->isLeaf();
	std::vector<ulong> old_labels;
	Node_* old_children;
	value_t* old_values;

	if(isLeaf()){
		old_values = _data._values;
		_data._values = nullptr;
	}else{
		old_children = _data._children;
		_data._children = nullptr;
	}
	_data._children = nullptr;

#ifdef STORE_SIZE
	auto oldcard = _cardinalities; // move?
#endif

	int p;
	int mid = find_mid(p);
	trace2("found mid", mid, w);

	old_labels.push_back(label(0));
	for(int i=1; i < old_size; ++i) {
		assert(node->label(i)>node->label(i-1));
		old_labels.push_back(label(i));
	}

	ulong m = ulong(-1) >> ( 64 - w );

	trace2("resize", std::bitset<64>(m), w);
	node = nullptr;

	int hi = w;
	m = ulong(-1) >> ( 64 - hi );

	int leftntz = _ntz + mid;
	int leftwidth = _width - mid;
	ulong leftmask = ulong(-1) << leftntz;

	int rightntz = _ntz;
	int rightwidth = mid;
	ulong rightmask = ulong(-1) >> ( 64 - (_ntz+rightwidth) );
	trace1("resize  ", std::bitset<64>(oldmask));
	trace1("resize", std::bitset<64>(rightmask));
	trace1("resize ", std::bitset<64>(leftmask));

	ulong collect = (label(0) >> leftntz);
	trace1("collect push1", std::bitset<64>(collect));

	// clear stuff ===========
	assert(old_labels[0] >> leftntz == collect);
#ifdef STORE_SIZE
	_cardinalities.clear();
#endif
	free(_labels);
	_labels = nullptr;
	_size = 0;
//	_cap = 0;

	_width = leftwidth;
	_ntz = leftntz;
	trace1("resize", std::bitset<64>(getMask()));
	node_t* n1 = this;
	/// =========================

	auto ccl = collect << leftntz;
	assert(rightwidth);
	node_t* newc = newNode(_index, rightwidth, rightntz);
//	auto newc = dynamic_cast<node_t*>(newc_);
//	assert(newc);
	n1->add(ccl, newc);

	for(unsigned i=1; i < unsigned(old_size); ++i){
		ulong upper_bits = old_labels[i] >> leftntz;
		if(upper_bits == collect){
			trace1("collect nopush", std::bitset<64>(upper_bits));
		}else{
			collect = upper_bits;
			auto ccl = collect << leftntz;
			assert(rightwidth);
			node_t* newc = newNode(_index, rightwidth, rightntz);
			n1->add(ccl, newc);
			trace2("collect push", _index, rightntz);
			//trace2("collect push", rightwidth, std::bitset<64>(upper_bits));
		}
	}

	assert(std::popcount(m) == w);

	// m = m >> ntz;
	auto cc = n1->_data._children;
	assert(cc);
	int j = 0;

	trace2("distributing leaves", p, old_size);
	for(int i=0; i < old_size; ++i){
		ulong label = old_labels[i]; //  << oldntz;
		assert(label==(label&oldmask));
		// collect?
		ulong rightlabel = rightmask & label;
		ulong leftlabel = leftmask & label;

		if(leftlabel != n1->label(j)){
			++cc;
			++j;
		}else{
		}
		assert(leftlabel == n1->label(j));

		if(!was_leaf){
			assert(i<old_size);
			cc->add(rightlabel, &old_children[i]);
			assert(!cc->isLeaf());
		}else{
#ifdef STORE_SIZE
			cc->add(rightlabel, old_values[i], oldcard[i]);
#else
			cc->add(rightlabel) = std::move(old_values[i]);
#endif
			assert(cc->isLeaf());
		}
		assert(cc->size() <= MAX_CHILDREN_SIZE);
	}
	trace2("resize post distrib", size(), old_size);

	if(was_leaf){
		free(old_values);
	}else{
		free(old_children);
	}
//	free(old_labels);
}
#endif

