/* Copyright (C) 2017 Felix Salfelder
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
 * Hisao Tamaki tree decomposition kernel. C++ rewrite.
 */
/*--------------------------------------------------------------------------*/
#ifndef TREEDEC_EXACT_TA_HPP
#define TREEDEC_EXACT_TA_HPP

#include "bits/trie.hpp"
#include "bits/predicates.hpp"
#include "graph.hpp"
#include "graph_util.hpp"

#include <gala/boost.h>
#include "graph_gala.hpp"

// #include "status.hpp"
#include "exact_base.hpp"
#include "trace.hpp"
#include <limits>
#ifdef STDHASH
#include <unordered_map>
#endif

#include <boost/graph/copy.hpp>

#undef tassert
#define tassert(x) assert(x)

#ifdef STATUS
long unsigned comp_size;
long unsigned delta_size;
long unsigned neigh_size;
#endif

// FIXME
#define NB_MAX (1l << 20)
#define HASH_FACTOR 4
#define TRIE_FACTOR 50

#define EXTA_t template<class G, template<class GG, class ...> class CFGT>
#define EXTA_a G, CFGT

#define stcnt(x)

namespace treedec{

// cleanup later.
template<class B, class T>
void merge(B&b, T const& s)
{
	for (auto const&v : s) {
		treedec::insert(b, v);
	}
}

namespace detail{

template<class S>
struct incidence_mask{
	incidence_mask(S& s) : _s(s) {}
	bool operator()(unsigned x) const{
		bool ret=_s.find(x)==_s.end();
		assert(ret == !_s.contains(x));
		trace3("imask", x, ret, _s);
		return ret;
	}
	bool operator[](unsigned x) const{
		return operator()(x);
	}
	void visit (unsigned x){
		_s.erase(x);
	}
	size_t size(){ return -1u; }
	S& _s;
};

}

template<class S>
detail::incidence_mask<S> make_incidence_mask(S& s)
{
	return detail::incidence_mask<S>(s);
}

/// treedec::algo::default_config?
template<class GraphType>
struct ta_config_default : treedec::algo::default_config<GraphType> {
	typedef typename boost::graph_traits<GraphType>::vertices_size_type vst;
	static constexpr unsigned max_vertex_index=std::numeric_limits<vst>::max();
};
/*--------------------------------------------------------------------------*/
template<class T, class S>
void assign_delta(T& t, S& s){
	clear(t);
	for(auto x : s){
		cbset::insert(t, x);
	}
}
/*--------------------------------------------------------------------------*/
namespace hack{
template<class T>
struct help{
	static void* next(T const& sth)
	{
		return (void*)(&sth + 1);
	}
};
#ifndef NOSTUFFD
template<unsigned w, class c>
struct help<  cbset::BSET_DYNAMIC<w, c> >{
	static void* next(cbset::BSET_DYNAMIC<w, c> const& sth)
	{ itested();
		return sth.next();
	}
};
#endif
template<class T>
void* next(T const& sth)
{
	return help<T>::next(sth);
}

} // hack
/*--------------------------------------------------------------------------*/
namespace bits{

// preallocated fixed vector.
// substitute for std::vector
template<class T>
class fvec{
public:
	typedef T value_type;
public:
	class const_iterator{
	public:
		const_iterator(fvec const& s, unsigned p=0)
		    : _s(s), _pos(p) {
		}
	public:
		void operator++(){ ++_pos; }
		bool operator!=(const const_iterator& o){
		  	return _pos!=o._pos;
		}
		T const& operator*(){return _s._data[_pos]; }
	private:
		const fvec& _s;
		unsigned _pos;
	};
public:
	fvec(unsigned n)
	    : _data(new T[n]),
		   _pos(0){
	}
	~fvec(){
		delete[] _data;
	}
	void clear(){
		_pos = 0;
	}
	unsigned size() const{
		return _pos;
	}
	void pop_back(){ untested();
		--_pos;
	}
	void push_back(unsigned x){
		_data[_pos++] = x;
	}
	value_type& back(){ untested();
		assert(_pos);
		return _data[_pos-1];
	}
	const_iterator begin()const{
		return const_iterator(*this);
	}
	const_iterator end()const{
		return const_iterator(*this, _pos);
	}

private:
	T* _data;
	unsigned _pos;
};
}
/*--------------------------------------------------------------------------*/
template<unsigned K, class CHUNK_T>
using bsd=cbset::BSET_DYNAMIC<K, CHUNK_T, cbset::nohowmany_t, cbset::nooffset_t, cbset::nosize_t>;
// using bsd=cbset::BSET_DYNAMIC<K, CHUNK_T, cbset::nohowmany_t, cbset::nooffset_t, unsigned>;

/*--------------------------------------------------------------------------*/
template<class G, template<class G_, class ...> class CFGT=ta_config_default>
class exact_ta{ //
public: // types
	typedef CFGT<G> CFG;
	constexpr static unsigned L=CFG::max_vertex_index;
#ifdef DYNAMIC
	typedef uint8_t CHUNK_T; // pick from config?
	constexpr static unsigned K=unsigned(L/8/sizeof(CHUNK_T)+1); // that many chunks
	typedef cbset::BSET_DYNAMIC<K, CHUNK_T, uint8_t, uint8_t, uint8_t /*BUG*/> T;
	typedef cbset::BSET_DYNAMIC<K, CHUNK_T, uint8_t, uint8_t, cbset::nosize_t> Tns;
#else
 	typedef uint64_t CHUNK_T; // pick from config?
	constexpr static unsigned K=unsigned(L/8/sizeof(CHUNK_T)+1); // that many chunks
	typedef bsd<K, CHUNK_T> T;
	typedef bsd<K, CHUNK_T> Tns;
#endif
	typedef typename T::value_type vertex_t;
	template<class A, class...>
	using myset=T;
	typedef gala::graph<myset, std::vector, unsigned> graph_type;

#ifdef STDHASH
	// incomplete
	typedef std::unordered_map<T, long> stlHashMap;
#endif
	class BLOCK {
	public:
		T component;
	private:
		T _neighbors;
		Tns _delta;
	public:
		Tns& delta_hack(){ untested();
		  	return _delta;
		}
	public: // stuffed access
		T const& neighbours()const{ return *(T*)hack::next(component); }
		Tns const& delta()const{ return *(Tns*)hack::next(neighbours()); }

		template<class DS>
		void stuff(T const& c, T& n, DS const& d){
			tassert(component==c); // still
			// hmm efficient? trim here?
			cbset::trim(n);
			auto nsize=uintptr_t(hack::next(n))-uintptr_t(&n);
			memcpy(hack::next(component), &n, nsize);

			Tns* nd=new (hack::next(neighbours())) Tns();
			assign_delta(*nd, d);
		}
		BLOCK* next() {
			return (BLOCK*)hack::next(delta());
		}
//		unsigned stuffedsize() const { untested();
//			return sizeof(3+3+2+component.howmany() +
//					neighbours().howmany() + delta().howmany());
//		}
	};
	typedef struct { //
		BLOCK const* bi;
	} ENTRY;

   static constexpr unsigned node_size=TRIE<T, BLOCK*>::node_size;
	typedef TRIE<T, BLOCK*, TRIE_SHARED_AREA<node_size> > trie_t;
	typedef typename trie_t::const_iterator NODE;

#ifndef STDHASH
typedef ENTRY* hashMap;
#else
typedef stlHashMap hashMap;
#endif
public:
	exact_ta()
	    : _blockmem(NULL), _top_block(NULL)
	{ untested();
		//incomplete(); // FIXME
		//_trie(boost::num_vertices(g))
		allocate();
		cbset::fullSet(all, n());
		cbset::clear(_emptyblock.component);
	}
	~exact_ta(){
		deallocate();
		::free(_blockmem);
	}
	template<class oG, class M=boost::identity_property_map>
	exact_ta(oG const& g, M const& m=boost::identity_property_map() /* =vertex_index_map?! */);
public:
	template<class T>
	bool try_it(T& t, unsigned x){ untested();
		if(_bag_size){ untested();
			incomplete();
			// retrying?!
			assert(x=_bag_size+1);
			assert(!solution);
			//	retry_decompose(x);
		}else{
		}

		if(try_decompose(x)){ untested();
			make_td(t);
			return true;
		}else{ untested();
			return false;
		}
	}
#if 0
	TD& get_treedec()
	{ untested();
		make_td(_t);
		return _t;
	}
#endif
#if 0
	void get_treedec(TT t, map)
	{ untested();
		make_td(t);
	}
#endif

private:
	size_t n() const{
		return _g.num_vertices();
	}
	void allocate();
	void deallocate();
	void clear();
	void clear_tries();
	void process(BLOCK *b);
	template<class N, class D>
	void registerBlock(N const& c, N /*const*/& onb, D const& delta);
	void try_combine(NODE *node, vertex_t v, T const& c, T const& neighb, NODE* from,
			unsigned last);
	void extendBy(NODE *node, vertex_t v, T const& c, T const& neighb, NODE* from,
			unsigned last=-1u);
	void extendByIterative(NODE *node, vertex_t v, T const& c, T const& neighb, NODE* from);
	bool is_valid(vertex_t v) const{ return v!=-1u && v<n(); }

private: // graph access. ideally forward to treedec
	const T& out_edges(vertex_t v) const{
		return _g.out_edges(v);
	}
	void merge_neighbors(T& s, vertex_t v) const{
		cbset::unionWith(s, out_edges(v));
	}
	T make_closed_neigh(T const& t) const{
		T tmp(t);
		treedec::close_neighbourhood(tmp, _g);
		return tmp;
	}
	T make_open_neigh(T const& t) const{
		T tmp(t);
		treedec::open_neighbourhood(tmp, _g);
		return tmp;
	}
	T make_saturation(T const& t) const{
		T tmp(t);
		treedec::saturate(tmp, _g);
		return tmp;
	}
	std::pair<bool, vertex_t> is_saturated(T const& t) const{
		T onb(make_open_neigh(t));
		return is_saturated(t, onb, onb);
	}
	std::pair<bool, vertex_t>
	    is_saturated(T const& vertices,
	                 T const& onb, vertex_t t=-1) const
	{
		assert(onb==make_open_neigh(vertices));
		return is_saturated(vertices, onb, onb, t);
	}
	std::pair<bool, vertex_t> is_saturated(
		T const& vertices,
		T const& onb,
		T const& cand,
		vertex_t t=-1u) const
	{
	  assert(onb==make_open_neigh(vertices));
	  T cnb=cbset::union_(vertices, onb);
	  assert(cnb==make_closed_neigh(vertices));
	  assert(isSubset(cand, onb));

#if 0 //  STRANGEHACK, not here.
	  // always do stuff, if v is an absorbable
	  if(t==-1u){
	  }else if(cbset::isSubset(_neighbours[t], cnb)) {
		  return std::make_pair(false, t);
	  }
#endif
	  (void) t;

	  for (auto const&v : cand) {
		  // if(t!=v){ untested(); << = slow
		  // }else
		  if(cbset::isSubset(out_edges(v), cnb)) {
			  return std::make_pair(false, v);
		  }
	  }
	  return std::make_pair(true, -1u);
	}
private: // rewrite
	void q_base_sets();
	void q_base_set(vertex_t v);
	typedef typename trie_t::const_iterator trie_citer;
	void extendByNew(trie_citer const& node, // ??
	                 vertex_t v,
	                 T const& c, T const& neighb,
	                 unsigned e=-1);
	void extendByNew(trie_citer const& node, // ??
	                 vertex_t v,
	                 BLOCK const& b,
	                 unsigned e=-1)
	{
		return extendByNew(node, v, b.component, b.neighbours(), e);
	}
	template<class NI>
	void try_combine_new(NI const& node,
			vertex_t v,
			T const& c,
			T const& neighb);
	template<class NI>
	void try_extend_union(NI const& node,
			T const& uc,
			T const& un,
			vertex_t v,
			T const* nnn=NULL);
	template<class S=T>
	void try_extend_by_vertex(T const& uc,
	                          T const& un,
	                          vertex_t v,
	                          S const* comm=NULL);
private: // hash
	void free(ENTRY*& x);
	void alloc(ENTRY*& x);
	void clearh(ENTRY*& x);
private:
	template<class D, class S=T>
	bool resaturate(T& c,
	                T const& onb, vertex_t v, T& cand, D& delta,
	                S const* comm=NULL);
	void registerForVertex(vertex_t v, BLOCK *block);
public:
	void do_it(unsigned n=2);
	template<class T>
	void do_it(T&, unsigned& bagsize);
	template<class t>
	void make_td(t&) const;
private:
	bool try_decompose(unsigned bs); // try_it?!
	template<class TREEDEC_>
	unsigned make_td(BLOCK const* b, TREEDEC_* td) const;
private: // i/o
	graph_type _g;
	unsigned _trieMax; // obsolete
   TRIE_SHARED_AREA<node_size> _shared_trie_area;
	std::vector<trie_t> _trie;
	typename trie_t::range_scratch_type _range_scratch;
private: //
	unsigned _bag_size;
	bits::fvec<unsigned> _delta;
private: // blockstack
	BLOCK* _blockmem;
	BLOCK* _top_block;
	BLOCK* lastblock;
	hashMap hashTable;
	BLOCK const* solution;
	BLOCK _emptyblock;
	BLOCK& top_block(){
		assert(_blockmem);
		assert(_top_block);
		return *_top_block;
	}
private: // hmm
	T all;
	T _empty;
	T _singleton;
	size_t _nHash;
private: // here?
	BLOCK const* getEnd(ENTRY const*)
	{ untested();
	  return NULL;
	}
	template<class S>
	BLOCK const*& getHashSpot(ENTRY* x, S component)
	{
		unsigned long h = cbset::hash(component) % _nHash;
//		trace1("gh", h);
		while (x[h].bi) {
//			trace1("eq?", component);
			if(x[h].bi->component == component){
				stcnt(ST_coll);
				break;
			}else{ untested();
				h = (h + 1) % _nHash;
			}
		}
		return x[h].bi;
	}
	// move to hash header?
	// use getHashSpot.
	template<class S>
	ENTRY const* getHash(ENTRY const* x, S component) const
	{
	  trace1("getHash", component);
	  unsigned long h = cbset::hash(component) % _nHash;
//	  trace1("gh", component);
//	  trace1("gh", h);
//	  trace1("gh", nHash);
//	  assert(h<x->size());
	  while (x[h].bi) { // }
		 trace1("eq?", component);
		 if(x[h].bi->component == component){
			return &x[h];
		 }else{ itested();
			h = (h + 1) % _nHash;
		 }
	  }
	  return NULL;
	}
	template<class S>
	ENTRY* getHash(ENTRY* x, S component)
	{ untested();
	  unsigned long h = cbset::hash(component) % _nHash;
	  trace1("gh", _nHash);
	  trace1("gh", component);
	  trace1("gh", h);
	  while (x[h].bi) { // }
		 trace1("eq?", component);
		 if(x[h].bi->component == component){ untested();
			return &x[h];
		 }else{ itested(); // 4x6
			h = (h + 1) % _nHash;
		 }
	  }
	  return NULL;
	}

	template<class KEY>
	void putHash(ENTRY* hp, KEY const& component, long bi, BLOCK* blocks);
	/// ...
	BLOCK const* getBlock(ENTRY const* x, T c) const
	{
		trace1("getblock1", c);
		auto hash=getHash(x, c);
		assert(hash);
		BLOCK const* b(hash->bi);
		(void)b;
		assert(b->neighbours()==make_open_neigh(b->component));
		trace1("gh", b->component);
		return hash->bi;
	}
	BLOCK* getBlock(ENTRY const* x, T c)
	{ untested();
		trace1("getblock2", c);
		auto hash=getHash(x, c);
		assert(hash);
		assert(hash->bi);
		return hash->bi;
	}
#ifdef STDHASH
	template<class X, class S>
	static void alloc(std::unordered_map<X, S>& m)
	{ untested();
	}
#endif
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EXTA_t
template<class oG, class M>
exact_ta<EXTA_a>::exact_ta(oG const& g, M const& m)
    : _g(),
      _trie(boost::num_vertices(g),
            trie_t(boost::num_vertices(g), _shared_trie_area)),
      _range_scratch(make_range_scratch(_trie[0]
#ifdef MORESCRATCH
					, _g.degree()
#endif
					)),
      _bag_size(0),
		_delta(boost::num_vertices(g)),
      _blockmem(NULL),
      _top_block(NULL)
{

	auto nv=boost::num_vertices(g);
	auto ne=boost::num_edges(g);

	//treedec::util::edge_map EM(m);
	auto EM=treedec::util::make_edge_map(m, g);
	auto g_edges=boost::edges(g);

	auto edgbegin=boost::make_transform_iterator(g_edges.first, EM);
	auto edgend=boost::make_transform_iterator(g_edges.second, EM);

	CFG::message(0, "graph size %d of %d\n", nv, L+1);
	CFG::message(0, "using %d chunks of size %d bit\n", K, 8*sizeof(CHUNK_T));

	_g = graph_type(edgbegin, edgend, nv, ne);
	assert(nv==boost::num_vertices(_g));
	assert(ne==boost::num_edges(_g));

#ifdef DYNAMIC
	CFG::message(0, "dynamic mode\n");
#else
	CFG::message(0, "nondynamic mode\n");
#endif
	allocate();
	cbset::fullSet(all, n());
}
/*--------------------------------------------------------------------------*/
namespace hack{
template<class T>
bool is(T const* a, T const* b){ untested();
	return a==b;
}
template<class T, class S>
bool is(T const* a, S const* b){ untested();
	return a->is(b);
}
}
/*--------------------------------------------------------------------------*/
template<class S>
static inline bool contains(S const& s, typename S::value_type i)
{
//  assert(i<::n);
  return s.contains(i);
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class D, class S>
inline bool exact_ta<EXTA_a>::resaturate(
		T& c,
		T const& onb,
		vertex_t v,
		T& newonb,
		D& delta,
		S const* comm)
{
	stcnt(ST_sat);
	/// v is an absorbable
	tassert(_bag_size>=size(onb));
	assert(!delta.size());
	tassert(!cbset::contains(c, v));
	// tassert(c==saturate(c)); not necessary.
	// tassert(onb == openNeighbor(c)); // not yet.
	T cnb = cbset::union_(onb, c);
	tassert(cnb == make_closed_neigh(c));
	//  unionWith(cnb, onb);

	merge_neighbors(cnb, v);

	tassert(cbset::contains(onb, v));
	tassert(!cbset::contains(c, v));
	cbset::add(c, v);
	tassert(equals(cnb, make_closed_neigh(c)));

	// 2-neighs of v can also be candidates.
	// ... if they are neighs of c.
	newonb = cnb; // cbset::union_(onb, neighborSets[v]);

	// we got them already.
	newonb.carve(c);
	// now newonb= onb(c_in + v)

//	assert(equals(newonb, _g.openNeighbor(c))); // no, not trimmed.

	if(cbset::size(newonb) + 1 > _bag_size){
		return false;
	}else{
	}

	// newonb is not trimmed...
	for (auto const&w : newonb) {
		assert(w!=v);
		if (cbset::isSubset(_g.out_edges(w), cnb)) {
			// found absorbable.
			if(!comm){
				stcnt(ST_now);
				// no whitelist. accept all.
			}else if(hack::is(comm, &_singleton)){ untested();
				stcnt(ST_by3);
				return false; // no, comm is false if saturated...
			}else if(comm->size()==1u){ untested();
				unreachable();
				assert(!cbset::contains(*comm, w));
				stcnt(ST_by0);
				return false;
			}else if(cbset::contains(*comm, w)){ untested();
			}else{ untested();
				stcnt(ST_by1);
				return false;
			}
			delta.push_back(w);
			tassert(!cbset::contains(c, w));
		}else{
		}
	}

	assert(delta.size()<=newonb.size());
	remove_sorted_sequence(newonb, delta); // do this later?
	add_sorted_sequence(c, delta);
	delta.push_back(v); // already in c

	tassert(cbset::contains(c,v));
	// tassert(equals(newonb, _g.openNeighbor(c))); no. not trimmed.
	assert(is_saturated(c).first);
	assert(cbset::size(newonb) + delta.size() <= _bag_size);
	return true;
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class KEY>
void exact_ta<EXTA_a>::putHash(ENTRY* hp, KEY const& component, long, BLOCK* block)
{ untested();
	trace1("putHash", component);
  unsigned long h=cbset::hash(component);

  // while (!equals(hp[h].key, empty))
  while (hp[h].bi) { itested(); // 4x6
    if (equals(hp[h].bi->component, component)) { untested();
		 // collision?!
		 unreachable();
      hp[h].bi = block;
      return;
    }else{ itested(); // 4x6
      stcnt(ST_coll);
    }
    h = (h + 1) % _nHash;
	 // incomplete. better prefetch.
    __builtin_prefetch (&hp[h], 1);
  }
  hp[h].bi = block;
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::free(ENTRY*& x)
{
  ::free((void*)x);
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::alloc(ENTRY*& x)
{
  x = (ENTRY *) malloc(_nHash * sizeof(ENTRY));
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::clearh(ENTRY*& x)
{
  memset(x, char(0), _nHash * sizeof(ENTRY));
}
/*--------------------------------------------------------------------------*/
EXTA_t
void exact_ta<EXTA_a>::deallocate()
{

  free(hashTable);
#ifdef VERBOSE
  printf("deallocated\n");
#endif
}
/*--------------------------------------------------------------------------*/
EXTA_t
void exact_ta<EXTA_a>::allocate()
{
	if(n()){
	}else{ untested();
	}

	size_t nbMax = NB_MAX;
	long* trial; 
	while (1) {
		_trieMax = nbMax * TRIE_FACTOR;
		_nHash = nbMax * HASH_FACTOR - 1;
		long total = _trieMax  * node_size + 
			_nHash * sizeof(ENTRY) + 
			nbMax * sizeof(BLOCK) + 
			// n() * sizeof(BAG) + 
			// n() * sizeof(TDEDGE) + 
			10 * n() * sizeof(int) + /* for bag entries */
			30 * n() * sizeof(void*); /* for stacks for extendByIterative */
		trial = (long int*) malloc(total);
		if (trial != NULL){
			break;
		}
		nbMax /= 2;
	}
	::free(trial);

	// BUG, proper alloc
	assert(!_blockmem);
	_blockmem = (BLOCK *) calloc(nbMax, sizeof(BLOCK));
	lastblock = (BLOCK*)( uintptr_t(_blockmem) + (nbMax-1)*sizeof(BLOCK) );

	//_trieArea = (NODE *) malloc(trieMax * sizeof(NODE) );
	_shared_trie_area.reserve(_trieMax);
	alloc(hashTable);
	trace1("allocated", nbMax);
}
/*--------------------------------------------------------------------------*/
EXTA_t
void exact_ta<EXTA_a>::clear_tries()
{
	_shared_trie_area.deallocate();

	for(auto i=_trie.begin(); i!=_trie.end(); ++i){
		i->clear();
	}
}
/*--------------------------------------------------------------------------*/
EXTA_t
void exact_ta<EXTA_a>::clear()
{
	assert(_blockmem);
// 	memset(_blockmem, char(0), nbMax * sizeof(BLOCK));
	_top_block = (BLOCK*)_blockmem;
	clear_tries();
	clearh(hashTable);
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class N, class D>
void exact_ta<EXTA_a>::registerBlock(N const& component, N& onb, D const& delt)
{
//  T const& component
  top_block().component = component;
  auto const& delta=delt;
  tassert(cbset::size(onb) + delta.size() <= _bag_size);

#ifdef DEBUG_NETYET // uses _g...
  itested();
  T components[n()];
  unsigned m=_g.getComponents(component, components);
  (void)m;
  tassert(m==1);
#endif

  trace1("ghs", component);
  BLOCK const*& NH=getHashSpot(hashTable, component);

  if (NH){
//	 trace1("REG0", component);
    stcnt(ST_alr);
    // already queued.
	 return;
  }else if (_top_block > lastblock) { itested();
    fprintf(stderr, "block area exausted\n");
/// bug. throw
    exit(1);
  }else{
//	  trace1("REG1", component);
    stcnt(ST_reg);
#ifdef STATUS
    comp_size+=(component.size());
    delta_size+=(delta.size());
    neigh_size+=(onb.size());
#endif

    unsigned cs=size(component);

    if(cs+_bag_size >= n()){
		 trace1("SOL", n());
		 trace1("SOL", component);
		 trace1("SOL      ", onb);
		 //trace1("SOL    ", delta);
      if(solution){
			incomplete();
		}
      solution = &top_block();
		trace1("found solution", solution->component);
    }else{
    }

    tassert (component.size());
//    tassert(equals(neighbors, openNeighbor(component)));
    tassert(top_block().component == component);
	 NH = &top_block();

	 top_block().stuff(component, onb, delta);
    assert(equals(top_block().neighbours(), onb));
    tassert(equals(top_block().neighbours(), make_open_neigh(component)));
	 _top_block = _top_block->next();
    // nb++;
    __builtin_prefetch (&top_block(), 1);
  }
}

// register comp into neigh tries.
// then trie to extend for all neighbors
EXTA_t
void exact_ta<EXTA_a>::process(BLOCK *b)
{
	tassert(equals(b->neighbours(), make_open_neigh(b->component)));
//  trace1("PROC", b->component);
  tassert(b);
  tassert(!solution);
  assert(b->neighbours()==b->neighbours());
  for(auto const&v : b->neighbours()) {
	  assert(!solution);

	 // trie[v].insert(b);
	  if(solution){ untested();
		  break;
	  }

    T c; // (top_block().component);
    T nn; //=top_block().neighbours_hack();
#if 0
	 std::vector<unsigned> delta;
#endif
	 _delta.clear();
    c = b->component;

#ifdef DEBUG
    T dbugc(c);
    cbset::add(dbugc, v);
    dbugc = make_saturation(dbugc);
    // trace1("", dbugc);
#endif
    tassert(equals(b->neighbours(), make_open_neigh(c)));

    cbset::clear(nn);
    stcnt(ST_proc);
    /// hmm c is already sat
    bool x=resaturate(c, b->neighbours(), v, nn, _delta);
	 assert(!x || is_saturated(c).first);
	 tassert(!solution);

    // pass parent? here "b".
    if (x){
      // tassert(nn==_g.openNeighbor(c)); not trimmed!
		// delta: the verts added by saruration.
      tassert(cbset::size(nn) + _delta.size() <= _bag_size);
      registerBlock(c, nn, _delta);
		if(solution){ untested();
			break;
		}else{
		}
    }else{
		 tassert(!solution);
	 }

    tassert(is_saturated(b->component).first);
    tassert(is_saturated(b->component, b->neighbours()).first);
	 T empty;
	 cbset::clear(empty); // why?
	 // assert(*b==*b);
	 extendByNew(_trie[v].begin() /* ->right */, v, *b);

    if(solution){
      return;
    }else{
	 }
	 _trie[v][b->component] = b;
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::q_base_set(vertex_t v)
{
	T c; // (top_block().component);
	//			Tns delta; // (top_block().delta);
	_delta.clear();
	T onb; // =top_block().neighbours_hack();
	cbset::clear(c);

//	trace1("cl", c);

#ifdef noRESAT // not yet.
	resaturate(c, c, v);
#else
	cbset::add(c, v);
//	trace1("added", c);

	T c1 = make_saturation(c);
//	assert(c1.size()<=::n);
	c = c1;
#endif

	tassert(c1 == make_saturation(c1));
	onb = make_open_neigh(c1);
	// delta = c1;

//	trace1("xx", c1);
	// parent. here nULL?
	// oops?!
	// assert(nn==nn);
	// BUG: check size.
	if (cbset::size(onb) + _delta.size() > _bag_size) { untested();
	}else{
		registerBlock(c, onb, _delta);
		if(solution){ untested();
			// BUG
			// sometimes looking for tw=5 finds a td of width 4. why?!
			return;
		}else{
		}
	}
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::q_base_sets()
{
//	for ( vertices... )
	for (unsigned v=0; v<n() && !solution; v++) {
		if (boost::out_degree(v, _g) < _bag_size) {
			q_base_set(v);
		}else{
		}
	}
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline bool exact_ta<EXTA_a>::try_decompose(unsigned bs)
{
	if(n()>L+1){ untested();
		std::cerr<<"too big: " << n() << "(" << L+1 << ")\n";
		unreachable();
		exit(25);
	}else{
	}

	assert(_blockmem);
	if(_bag_size+1==bs){
		// reuse.
		// fscs are still what they are.
		clear_tries();
	}else{
		clear();
	}
//		clear(); // !
	_bag_size = bs;

	// trace1("TRY", _bag_size);
	// tassert(_top_block==_blockmem); no. subsequent tries

	fprintf(stderr, "try bagsize = %d\n", _bag_size);

	q_base_sets();

	BLOCK* bb=_blockmem;
	while(bb!=&top_block() && !solution){
		// try to extend newly discovered blocks.
		process(bb);
		bb = bb->next();
	}
	return solution;
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class TT>
inline void exact_ta<EXTA_a>::do_it(TT& t, unsigned& bs)
{
	assert(bs>1);

	if(_bag_size){ untested();
		// retrying?
	}else{
	}

	do_it(bs);

	bs = _bag_size;
	make_td(t); // bug. not here.
}
/*--------------------------------------------------------------------------*/
EXTA_t
inline void exact_ta<EXTA_a>::do_it(unsigned bs)
{
	assert(bs>1);
	solution = NULL;
	while (!solution) {
		try_decompose(bs);
		if(!solution){
			++bs;
		}
	}
	assert(solution);
}
/*--------------------------------------------------------------------------*/
#include <boost/graph/graph_traits.hpp>
#include "graph_traits.hpp"
/*--------------------------------------------------------------------------*/
EXTA_t
template<class TREEDEC_t>
inline void exact_ta<EXTA_a>::make_td(TREEDEC_t& td) const
{
  auto component=solution->component;
  trace2("TDTDTDTD", component, cbset::size(component));
  //boost::clear(td);
  assert(boost::num_vertices(td)==0);

  tassert(_bag_size);

  if(n() - cbset::size(component)) {
    unsigned k=boost::add_vertex(td);
    auto& b=boost::get(bag_t(), td, k);
	 T s=cbset::diff(all, component);
    treedec::merge(b, s);
    unsigned j=make_td(solution, &td);
    boost::add_edge(k, j, td);
  }else{ untested();
    // all in one...
    make_td(solution, &td);
  }
  std::cerr << "make_td nvt " << boost::num_vertices(td) << "\n";
  assert(boost::num_edges(td)+1==boost::num_vertices(td));

#ifdef DEBUG
  auto e=boost::edges(td);
  for(;e.first!=e.second;++e.first){
	  std::cerr << boost::source(*e.first, td);
	  std::cerr << ":" <<  boost::target(*e.first, td) << "\n";
  }
#endif

#ifdef STATUS
  std::cout << "==== STATUS ====\n";
  std::cout << "csum " << comp_size << " " << cnt[ST_reg] << "\n";
  std::cout << "cavg " << comp_size / cnt[ST_reg] << "\n";
  std::cout << "dsum " << delta_size << " " << cnt[ST_reg] << "\n";
  std::cout << "davg " << delta_size / cnt[ST_reg] << "\n";
  std::cout << "nsum " << neigh_size << " " << cnt[ST_reg] << "\n";
  std::cout << "navg " << neigh_size / cnt[ST_reg] << "\n";
  std::cout << "==== STATUS ====\n";
  for(unsigned i=0; i<STNUM; ++i){ untested();
	  std::cout << i << ":" << stn[i] << " " << cnt[i] << "\n";
  }
  exit(0);
#endif
}
/*--------------------------------------------------------------------------*/
template<class T, class TREEDEC_>
static unsigned addBag(T s, TREEDEC_* td)
{
  unsigned k=boost::add_vertex(*td);
  auto& b=boost::get(treedec::bag_t(), *td, k);
  treedec::merge(b, s);
  return k;
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class TREEDEC_>
inline unsigned exact_ta<EXTA_a>::make_td(BLOCK const* block, TREEDEC_* td) const
{
  std::vector<BLOCK const*> bStack(n());
  std::vector<int> aStack(n());

  bStack[0]=block;
  aStack[0]=-1;
  int top=0;
  unsigned r=0;

  while (top >= 0) {
    trace1("", top);
    BLOCK const* b = bStack[top];
    int k = aStack[top];
    top--;
    if (cbset::size(b->component) + cbset::size(b->neighbours())
      <= _bag_size) {
      int j = addBag(cbset::union_(b->component, b->neighbours()), td);
      if (k >= 0) {
			boost::add_edge(k, j, *td);
      }else{ untested();
			r = j;
		}
      continue;
    }else{

	 }
    T s=cbset::union_(b->neighbours(), b->delta());

    unsigned j=addBag(s, td);
    if(k >= 0){
      boost::add_edge(k, j, *td);
    }else{
      r = j;
    }
    T base = cbset::diff(b->component, b->delta());
	 cbset::trim(base);
	 trace1("base", base);
	 auto roi(base);
	 assert(roi==base);

#if 1
	 auto vert=boost::vertices(_g);
	 auto mask=make_incidence_mask(base);
//    auto cmps_range=make_components_range(roi.begin(), roi.end(), _g, mask);
    auto cmps_range=make_components_range(vert.first, vert.second, _g, mask);

	 for(; cmps_range.first!=cmps_range.second; ++cmps_range.first){
		 top++;
		 auto comp_range=*cmps_range.first;

		 T cmp_i;
		 assert(cmp_i.empty());
		 unsigned cnt=0;
		 for(; comp_range.first!=comp_range.second; ++comp_range.first){
			 trace1("found", *comp_range.first);
			 cmp_i.add(*comp_range.first);
			 ++cnt;
		 }
		 assert(cnt==cmp_i.size());
		 bStack[top] = getBlock(hashTable, cmp_i);
		 if(!bStack[top]){ untested();
			 unreachable();
			 std::cerr << "something is wrong\n";
		 }else{
		 }
		 assert(bStack[top]);
		 aStack[top] = j;
	 }
#else
	 T components[n()];
    unsigned m=_g.getComponents(base, components);
    for (unsigned i=0; i<m; i++) {
      top++;
      bStack[top] = getBlock(hashTable, components[i]);
		if(!bStack[top]){ untested();
			unreachable();
			std::cerr << "something is wrong\n";
		}
      tassert(bStack[top]);
      aStack[top] = j;
    }
#endif
  }
  return r;
}
/*--------------------------------------------------------------------------*/
// only called by process
EXTA_t
inline void exact_ta<EXTA_a>::extendByNew(
		typename trie_t::const_iterator const& /*node*/, // needed?!
		vertex_t v,
		// BLOCK const& b,
		T const& c, T const& neighb,
		unsigned /*to*/)
{
//	tassert(is_connected(c));
	auto lu(make_lazy_union(c, neighb));
	auto P(::util::make_not_in_set(lu));

#if 0 // ndef NDEBUG
	{ untested();
		auto y=make_range(_trie[v]);
		trace1("loopstart", v);

		unsigned cnt=0;
		for(;y.first!=trie_t::end(); ++y.first){ untested();
			++cnt;
		}
		assert(cnt==_trie[v].size());
	}
#endif

	auto y=make_range(_trie[v], &P, _range_scratch); // P
	//trace1("loopstart", v);

	for(;y.first!=trie_t::end(); y.first.inc(&P)){
		tassert(y.first->block);
//		trace1("loop", y.first->block->component);

		auto const& iternode=y.first;
		tassert(cbset::contains(iternode->block->component, iternode.back()));

		tassert(!cbset::contains(c, iternode.back()));
		tassert(!cbset::contains(neighb, iternode.back()));

		tassert(!cbset::intersects(iternode->block->component, c));
		tassert(!cbset::intersects(iternode->block->component, neighb));

		assert(neighb==neighb);
#if 0
		try_combine_new(y.first, v, c, neighb, &neighb);
#else
		try_combine_new(y.first, v, c, neighb);
#endif
		if(solution){
		  	break;
		}
	}
}
/*--------------------------------------------------------------------------*/
// uc is saturated, or v is an absorbable.
// if comm, then give up if there's an absorbable not in comm
EXTA_t
template<class S>
inline void exact_ta<EXTA_a>::try_extend_by_vertex(
		T const& uc,
		T const& un,
		vertex_t v,
		S const* comm)
{
	T c2; //(top_block().component);
	T nn; // (top_block().neighbours_hack());
	_delta.clear();

	c2 = uc;
	T c2onb = un;

	stcnt(ST_ext);
	// hmm c2 is union...
	bool x=exact_ta<EXTA_a>::resaturate(c2, c2onb, v, nn, _delta, comm);
	if(x){
		registerBlock(c2, nn, _delta);
		if (solution) {
			return;
		}
		tassert(size(c2) == size(uc) + _delta.size());
	}else{
		tassert(!solution);
	}
}
/*--------------------------------------------------------------------------*/
EXTA_t
template<class NI>
inline void exact_ta<EXTA_a>::try_combine_new(
		NI const& node,
		vertex_t v, T const& c,
		T const& neighb)
{
	tassert(node->block);
	tassert(!cbset::contains(node->block->component, v));
	tassert(cbset::size(neighb) <= _bag_size);
	tassert(cbset::size(node->block->neighbours()) < _bag_size);

	tassert(!cbset::intersects(node->block->component, cbset::union_(c, neighb))
			&& "trie search inconsistent for %d\n");

	T un = cbset::union_(neighb, node->block->neighbours());
	assert(un==un);
	if(solution){
		unreachable();
	}

	if (cbset::size(un) > _bag_size) {
		// too many neigbours.
	}else{
		T uc = cbset::union_(c, node->block->component);
		assert(un == make_open_neigh(uc));

		// intersect with previous common neighbours.
		auto abs=is_saturated(uc, un, v);

		if(abs.first){
			// saturated.
			// try to add neighbours
			try_extend_by_vertex(uc, un, v);

			stcnt(ST_rec);
			try_extend_union(node.back(), uc, un, v);
		}else if	(abs.second==v){
			try_extend_by_vertex(uc, un, v);
		}else{
			// found an absorbable. but it's not v.
		}
	}
}

/*--------------------------------------------------------------------------*/
EXTA_t
template<class NI>
inline void exact_ta<EXTA_a>::try_extend_union(
		NI const&nodeback, // needed?!
		T const& uc,
		T const& un,
		vertex_t v,
		T const* //remove?
		)
{
	// avoid sets that intersect cnb(union)
	// only visit sets with lower vertex indexes
	//
	auto LU(make_lazy_union(uc, un));
	auto P(::util::make_not_in_set(LU));
	auto LT(::util::make_lt(nodeback));
	auto cond(::util::make_conj(LT, P));

#ifdef MORESCRATCH
	// incomplete(); might make sense for iterative code
	auto x=make_range(_trie[v], cond, node.scratch_tail());
#else
	auto x=make_range(_trie[v], &cond);
#endif

	for( ; x.first!=trie_t::end(); x.first.inc(&cond)){
		tassert(x.first->block);
		auto const& iternode=x.first;

		tassert(!cbset::contains(uc, iternode.back()));
		tassert(!cbset::contains(un, iternode.back()));

		tassert(!cbset::intersects(iternode->block->component, uc));
		tassert(!cbset::intersects(iternode->block->component, un));
		tassert(x.first->block->component.back() < nodeback);

		try_combine_new(x.first, v, uc, un);
		// hmm v is in the neighborhood in between
		// if others, w, are in the saturation, then sat({v}+uc) is computed twice.
		// also in in extendby(... w )
		if(solution){
			break;
		}
	}
}
/*--------------------------------------------------------------------------*/
} // treedec
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#undef stcnt

#endif // guard
