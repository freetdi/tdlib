// Felix Salfelder, 2016, 2017, 2021
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
// fill in heuristic (supergraph backend)

#ifndef TREEDEC_FILL_IN_SG_HPP
#define TREEDEC_FILL_IN_SG_HPP

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#error "not meant to be used like that."
#endif

#include "../algo.hpp"
#include "greedy_base.hpp"
#include "../message.hpp"
#include "../elim_util.hpp"

#ifndef NDEBUG
#define DEBUG_MARKER
#elif defined DEBUG_FILL
#define DEBUG_MARKER
#elif defined DO_TRACE
#define DEBUG_MARKER
#endif

static const auto minlong = -std::numeric_limits<long>::max();
// static const bool me = false;
static const bool me = 1;
static const bool me3 = 0;
static const bool me2 = 1 || me3;
static const bool deref = 0;
static const bool dc2 = 1; // drop cliques in p2
static const bool dc3 = 1; // drop cliques in p3
static const bool flb = 0; // flatten bags
static const bool czo = 0; // n has cliques only

static const bool tree = 1;

namespace treedec{
namespace pending{
namespace impl{
namespace detail{

template < class G, class MarkerP, class NumberD, class Stack,
  class VertexIndexMap >
class remove_and_collect {
	typedef typename boost::graph_traits< G >::vertex_descriptor vertex_t;
	typedef typename boost::graph_traits< G >::edge_descriptor edge_t;

public:
	remove_and_collect(G& g, MarkerP& _marker, NumberD const& _numbering,
	                   Stack& n_e, VertexIndexMap id)
	  : _g(g)
	  , marker(&_marker)
	  , numbering(_numbering)
	  , neighbor_elements(&n_e)
	  , id(id) {
	}

	bool operator()(edge_t e) {
		vertex_t s = boost::source(e, *_g);
		vertex_t t = boost::target(e, *_g);
		if(!_g.supernode_size(t)){
			return true;
//		}else if (marker->is_tagged(t)){ untested();
//			trace2("collect delete", s, t);
//			unreachable();
//			return true;
		}else{
			assert(!marker->is_tagged(t));
			marker->mark_tagged(t);
		}

		if(tree && _g.supernode_size(t)<0){
			assert (numbering.is_numbered(t));
			trace2("collect redundant", s, t);
			return true;
		}else if (numbering.is_numbered(t)) {
			neighbor_elements->push(get(id, t));
			trace2("collect delete numbered", s, t);
			return true;
		}else{
			assert( t != s );
			return false;
		}
	}

private:
	G& _g;
	MarkerP* marker;
	NumberD const& numbering;
	Stack* neighbor_elements;
	VertexIndexMap id;
	size_t _cnt{0};
}; //remove_and_collect



// count N(n) \ N(c)
// runs on non-unique neighbours of n.
// multiple-tags outbound neighbours.
// N(c) is tagged... delete edges to N(c)
template<class G, class M>
class predicate_scan_neigh {
public:
	typedef typename G::vertex_descriptor vertex_descriptor;
	typedef typename G::vertices_size_type vertices_size_type;

	predicate_scan_neigh(G const& g, M& m, vertex_descriptor c /*, H& h*/)
	  : _g(g), _fill_marker(m)
	  , _c(c)
#ifdef DEREF
	  , _stack_hack(h)
#endif
	{ itested();
	}

	template<class E>
	bool operator()(E e, bool do_cliques=true) const {
		auto& marker = _g._marker;
		auto s = boost::source(e, *_g);
		auto t = boost::target(e, *_g);

		trace5("p2 oc?", s, t, do_cliques, _g.is_numbered__(t), _only_cliques);
		if(!do_cliques){
			// assert(_g.is_numbered(s)); no. mmark doesn't filter (should?)
		}else if (_g.is_numbered__(t)){
		}else{
			_only_cliques = false;
			trace5("p2 oc", s, t, do_cliques, _g.is_numbered__(t), _only_cliques);
			// assert(!_g.is_numbered(t));
		}

		int reason = 0;
		if (_c == t){
			assert(do_cliques);
			reason = 99;
		}else if(is_covered(t)) {
			if(do_cliques){
				reason = 5;
			}else{ untested();
			  unreachable();
				// before/after?
				reason = -5;
			}
//        }else if (tree & supernode_size(t) <= 0 ){
		}else if (_g.is_numbered_(t)){
			assert(do_cliques);
			{
				assert(marker.is_extra(s)); // must be in N(c)...
				marker.mark_multiple_tagged(s);

				reason = mmark_clique(e);

#ifndef NDEBUG
				if(reason<=0){
					trace2("not dropping clique", s, t);
					auto vv = boost::adjacent_vertices(t, *_g);
					for(;vv.first != vv.second; ++vv.first){
						auto n= *vv.first;
						trace4(" --- ndc", s, t, n, _g.is_numbered(n));
					}

				}else{
					trace3("dropping clique", s, t, reason);
					auto vv = boost::adjacent_vertices(t, *_g);
					for(;vv.first != vv.second; ++vv.first){
						auto n= *vv.first;
						trace4(" --- dc", s, t, n, _g.is_numbered(n));
					}
				}
#endif
				if(dc2){
				}else{ untested();
					reason = -21;
				}
			}
		}else if (marker.is_extra(t)){
			// N(c) is extra
			assert(!_g.is_numbered(t));
			assert(marker.is_tagged(t));

			// extra < multi.
			// need tag to identify missing edges.
			trace2("p2 multitag", s, t);
			marker.mark_multiple_tagged(t);
			reason = 19;
		}else if (marker.is_tagged(t)){
			assert(!marker.is_multiple_tagged(t));

			// not in N(c), already been here
			reason = 42;
		}else{
			// untagged. just count.
			marker.mark_tagged(t);
			assert(_g.supernode_size(t)>0);
			trace2("p2 found dn", s, t);

			assert(!_g.is_numbered(t));
			assert(!marker.is_multiple_tagged(t));
			assert(!marker.is_extra(t));
			_cn += _g.supernode_size(t);
			++_rcn;
			reason = -12;
		}

		trace6("p2 delete", s, t, reason, do_cliques, _g.is_numbered__(t), _only_cliques);
		return reason>0;
	}

	bool only_cliques() const{
		return _only_cliques;
	}
	size_t cnt() const{
		return _cn;
	}
	size_t rcnt() const{
		return _rcn;
	}
	size_t rcl() const{
		return _rcl;
	}

	// debugging?
	void reset() {
		_cn = 0;
		_rcn = 0;
		_rcl = 0;
		//_ncc = 0;
#ifdef DEREF
		_stack_hack.clear();
#endif
		_only_cliques = true;
	}

private:
	template<class E>
	int mmark_clique(E const&) const;

	bool is_covered(vertex_descriptor t) const{
		if(tree){
			return _g.supernode_size(t) <= 0;
		}else{
			return _g.supernode_size(t) <= 0;
		}
	}

private:
	G const& _g;
	M& _fill_marker;
	mutable vertices_size_type _cn{0};
	mutable vertices_size_type _rcn{0};
	mutable vertices_size_type _rcl{0};
	vertex_descriptor _c;
	mutable bool _only_cliques;
public:
#ifdef DEREF
	std::vector<vertex_descriptor>& _stack_hack;
#endif
}; // predicate_scan_neigh


// mark N(n) clique component rooted at v
// todo: prune expired cliques.
template<class G, class M>
template<class E>
int predicate_scan_neigh<G, M>::mmark_clique(E const& e) const
{
	int ret = 17;

	auto s = boost::source(e, *_g);
	auto v = boost::target(e, *_g);
#ifndef NDEBUG
	auto& marker = _g._marker;
#endif
	assert( marker.is_multiple_tagged(s));
	assert(_g.is_numbered(v));
	trace1("p2 scan", v);

	int count = 0;
	unsigned rcl = 0;

	auto vv = boost::adjacent_vertices(v, *_g);
	vertex_descriptor which;
	auto& marker = _g._marker;

	for(; vv.first!=vv.second; ++vv.first) {
		++rcl;
		auto t = *vv.first;
		trace1("p2 mmark", t);

		auto p = std::make_pair(s, t);
		assert(boost::source(p, *_g)==s);
		assert(boost::target(p, *_g)==t);


	 	if(s==t){
//			++count;
//			which = t;
//			ret = -18;
		}else if(_g.supernode_size(t)<=0){
				//   marker.mark_multiple_tagged(t);
//			++count;
//			which = t;
//			ret = -18;
			//??

	 	}else if(_g.is_numbered(t)){
#if 1
		}else if (marker.is_extra(t)){
			// N(c) is extra
			assert(!_g.is_numbered(t));
			assert(marker.is_tagged(t));

			// extra < multi.
			// need tag to identify missing edges.
			trace2("p2 multitag", s, t);
			marker.mark_multiple_tagged(t);
			//ret = 19;
		}else if (marker.is_tagged(t)){
			assert(!marker.is_multiple_tagged(t));

			// not in N(c), already been here
			//ret = 42;
		}else{
			// untagged. just count.
			marker.mark_tagged(t);
			assert(_g.supernode_size(t)>0);
			trace2("p2 found dn", s, t);

			assert(!_g.is_numbered(t));
			assert(!marker.is_multiple_tagged(t));
			assert(!marker.is_extra(t));
			_cn += _g.supernode_size(t);
			++_rcn;
			ret = -12;
		}

#else
		}else{
			bool o = operator()(p, false);

			if(!o){
				++count;
				which = t;
				ret = -17;
			}else{
			}

			assert(t!=_c);
		}
#endif
	}

#ifdef DEREF
	if(count==1){
		auto deg = boost::out_degree(v, *_g);
		std::cout << "c deref clique " << v << "? " << deg << " " << which << "\n";
		_stack_hack.push_back(which);
		return 50;
	}
#else
	(void) which;
#endif

	{
		if(ret>0){
			trace2("unlink clique", v, ret);
		}else{
			trace2("no unlink clique", s, v);
			_rcl += rcl;
		}
		return ret;
	}
}

// delete out edges to done nodes
// multitag marked nodes, N(c).
// find overlap: multitag fill marked nodes (DN(n)) and count.
template<class G, class M, class F>
class predicate_collect_overlap {
private:
	typedef typename G::vertex_descriptor vertex_descriptor;
public:
	explicit predicate_collect_overlap(G const& g, M& m, F& f,
	                                   vertex_descriptor c,
	                                   vertex_descriptor n )
		: _g(g), _fill_marker(m), _fill(f), _c(c), _n(n) { itested();
	}
	template<class E>
	bool operator()(E e, bool do_cliques=true) const { itested();
		auto s = boost::source(e, *_g);
		auto t = boost::target(e, *_g);
		trace4("overlap collect N(s)", s, t, _c, do_cliques);
		int reason = 0; // don't delete.

		// auto ii = _supergraph.adjacent_vertices(*t);
		{ //
			// assert(i!=c || _marker.is_done(i));
			// N(s) should be tagged by now.
			//
			//assert(_g._marker.is_done(t) || _g.supernode_size(t)>0);
			assert(_fill_marker.is_done(t) || _g.supernode_size(t)>=0);

			if (_g._marker.is_done(t)){
				if(do_cliques){
					reason = 1;
				}else{ untested();
					unreachable();
					reason = -1;
				}
			}else if (_fill_marker.is_multiple_tagged(t)){
				assert(do_cliques);
				// been there.
				reason = 2;
			}else if (tree && _g.supernode_size(t)<0){ untested();
				assert(do_cliques);
				reason = -3;
			}else if(_g.is_numbered(t)) {
				assert(do_cliques);
				{

					_fill_marker.mark_multiple_tagged(s); // avoid revisit.
					trace2("p3 down", s, t);
					assert(_g.supernode_size(t)>=0);
					reason = mcount_overlap(e);

					if(dc3){
					}else{
						reason = -30;
					}
				}
			}else if (_g._marker.is_multiple_tagged(t)){
				// overlap inside N(c)
				if(me2){
					assert(_g.supernode_size(t)>0);
				}else{
					assert(_g.supernode_size(t)==1);
				}
				assert(_g._marker.is_extra(t));
				trace1("p3 fu shift internal", t);

				assert(_g.supernode_size(t)>=0);
				_ncoverlap += _g.supernode_size(t);

				_fill_marker.mark_multiple_tagged(t);

				assert(_fill.get_value(t));
				assert(_g.supernode_size(_c)>0);
				assert(_g.supernode_size(_n)>=0);
				assert(_g.supernode_size(s)>=0);

				long w = _g.supernode_size(s) * _g.supernode_size(_n);
				assert(long(_fill.get_value(t)) >= w);
				_fill.shift(t, -w);

				reason = -8;
			}else if(_g._marker.is_extra(t)){
				// nodes on N(c) not reached from n.
				// --> not interesting. still needed for subsequent scans.
				// extra < multi.
				assert(!_fill_marker.is_multiple_tagged(t));

				// delete in next go.
				_fill_marker.mark_multiple_tagged(t);
			}else if(_g._marker.is_tagged(t)){
				// external overlap
				_fill_marker.mark_multiple_tagged(t);
				assert(_g.supernode_size(t)>0); // for now

				_overlap += _g.supernode_size(t);
				reason = -2;

				trace1("p3 fu shift", t);

				assert(_g.supernode_size(_c)>0);
				assert(_g.supernode_size(s)>=0);

				auto w = _g.supernode_size(s) * _g.supernode_size(_n);
				assert(long(_fill.get_value(t)) >= w);
				_fill.shift(t, -w);

			}else{
				assert(!_fill_marker.is_multiple_tagged(t));
				assert(_g.supernode_size(t)!=0);

				_nc += _g.supernode_size(t);
				// _g._marker.mark_tagged(t);
				_fill_marker.mark_multiple_tagged(t);
				reason = -17;
			}
		}

		trace4("p3 delete", s, t, reason, do_cliques);
		return reason>0;
	}
	size_t ncoverlap() const{ untested();
		return _ncoverlap;
	}
	size_t overlap() const{
		return _overlap;
	}
	size_t cnt() const{
		return _nc + _overlap; // 3 ??
	}
	void reset() {
		_overlap = 0;
		_ncoverlap = 0;
		_nc = 0;
	}
private:
	template<class E>
	int mcount_overlap(E const&) const;
private:
	G const& _g;
	M& _fill_marker;
	F& _fill;
	//std::vector<int>& _fu;
	mutable size_t _nc{0};
	mutable size_t _overlap{0};
	mutable size_t _ncoverlap{0};
	vertex_descriptor _c;
	vertex_descriptor _n;
}; // predicate_collect_overlap

// mark ( N(n) \union N(c) ) \intersect N(m)
// multiple tags used for deduplicated counting.
// todo: prune expired cliques.
// member function?
// mmd does not nest calls, why?
template<class G, class M, class F>
template<class E>
int predicate_collect_overlap<G,M,F>::mcount_overlap(E const& e) const
{
	auto s = boost::source(e, *_g);
	auto v = boost::target(e, *_g);
	assert(s!=v);
//	assert(s==_n); no, m
	assert(_g.is_numbered(v));
	if(!_g.supernode_size(v)){
		assert(me3);
		return 1;
	}else{
	}

	auto vv = boost::adjacent_vertices(v, *_g);
	assert(v!=_c);

	int ret = 50;
	for(; vv.first!=vv.second; ++vv.first) {
		auto t = *vv.first;
		auto p = std::make_pair(s, t);

#ifndef NDEBUG
		if(_g.supernode_size(t)<0){
			assert(_g.is_numbered(t));
		}else{
		}
#endif

		if(_g.is_numbered(t)){
		}else if(_fill_marker.is_multiple_tagged(t)){
			// in N(c)?
		}else if (_g._marker.is_multiple_tagged(t)){
			// overlap inside N(c)
			if(me2){
				assert(_g.supernode_size(t)>0);
			}else{
				assert(_g.supernode_size(t)==1);
			}
			assert(_g._marker.is_extra(t));
			trace1("p3 fu shift internal", t);

			assert(_g.supernode_size(t)>=0);
			_ncoverlap += _g.supernode_size(t);

			_fill_marker.mark_multiple_tagged(t);

			assert(_fill.get_value(t));
			assert(_g.supernode_size(_c)>0);
			assert(_g.supernode_size(_n)>=0);
			assert(_g.supernode_size(s)>=0);

			long w = _g.supernode_size(s) * _g.supernode_size(_n);
			assert(long(_fill.get_value(t)) >= w);
			_fill.shift(t, -w);

			ret = -8;
		}else if(_g._marker.is_extra(t)){
			// nodes on N(c) not reached from n.
			// --> not interesting. still needed for subsequent scans.
			// extra < multi.
			assert(!_fill_marker.is_multiple_tagged(t));

			// delete in next go.
			_fill_marker.mark_multiple_tagged(t);
			ret = -5;
		}else if(_g._marker.is_tagged(t)){
			// external overlap
			_fill_marker.mark_multiple_tagged(t);
			assert(_g.supernode_size(t)>0); // for now

			_overlap += _g.supernode_size(t);
			ret = -2;

			trace1("p3 fu shift", t);

			assert(_g.supernode_size(_c)>0);
			assert(_g.supernode_size(s)>=0);

			auto w = _g.supernode_size(s) * _g.supernode_size(_n);
			assert(long(_fill.get_value(t)) >= w);
			_fill.shift(t, -w);

		}else{
			assert(!_fill_marker.is_multiple_tagged(t));
			assert(_g.supernode_size(t)!=0);

			_nc += _g.supernode_size(t);
			// _g._marker.mark_tagged(t);
			_fill_marker.mark_multiple_tagged(t);
			ret = -17;
		}

		// }else if (!operator()(p, false)){
		// 	// indicate nodelete.
		// 	ret = -50;

	}
	return ret;
} // predicate_collect_overlap::mcount_overlap

} // impl::detail

template<class V, class G>
static void push_front_edge(V N, V c, G& g)
{
#if 1
	auto f = g.out_edges(N).front();
	g.out_edges(N).push_back(f);
	g.out_edges(N).front() = c;
#else
	auto& E = g.out_edges(N);

	for(auto& e : E){
		std::swap(e,c);
	}
	
	E.push_back(c);
#endif
}
				//boost::add_edge(N, c, *_g); // **

// the fillIn heuristic.
template<typename G,
         template<class GG, class ...> class CFGT=algo::default_config>
class fillIn : public treedec::impl::greedy_base< G,
               std::vector<typename boost::graph_traits<G>::vertex_descriptor>,
               CFGT>{ //
public:
	typedef std::vector<typename boost::graph_traits<G>::vertex_descriptor> O;
private:
	typedef CFGT<G> CFG; // got through baseclass?
	typedef treedec::impl::greedy_base<G, O, CFGT> baseclass;
	// typedef typename baseclass::VertexIndexMap VertexIndexMap;
	typedef typename boost::property_map<G, boost::vertex_index_t>::type VertexIndexMap;
public: //types
	typedef O O_t; //?
	typedef G G_t; //?
	typedef typename baseclass::vertex_descriptor vertex_descriptor;
	typedef typename directed_view_select<G_t>::type D_t;
	typedef typename boost::graph_traits<D_t>::vertices_size_type vertices_size_type;
   // typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
	typedef long diff_t;
	typedef Marker< diff_t, vertex_descriptor, VertexIndexMap > marker_type;
	//    marker_type _fill_marker;
	typedef typename baseclass::numbering_type numbering_type;
	// BUG:: use CFGT::fill or fallback to current fill
	typedef typename fill_chooser<typename baseclass::supergraph_type>::type fill_type;

	typedef typename std::make_unsigned<
                          typename boost::graph_traits<G>::edges_size_type
                          >::type fill_value_t;

#ifdef DEBUG_FILL
	std::set<vertex_descriptor> debug_fill(vertex_descriptor c) {
		(void)c;
		std::set<vertex_descriptor> c_neigh;
#ifndef NDEBUG
//        typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
//        marker_type _debug_marker(_num_vert);

        auto cr = boost::adjacent_vertices(c, _supergraph);
        for(; cr.first!=cr.second; ++cr.first){
            auto n = *cr.first;
            if(n==c){
            }else if(_supergraph._marker.is_done(n)){ untested();
				}else{
                c_neigh.insert(n);
            }
        }
		  auto me = _supergraph.count_missing_edges(c, &_debug_marker);
		  trace4("DEBUG_FILL pre-elim", c, fill_cached_(c), fill_cached_is_lb_(c), me);
		  assert( fill_cached_(c) == me);

        for(auto n: c_neigh) {
            auto me = _supergraph.count_missing_edges(n, &_debug_marker);
            trace4("DEBUG_FILL pre-elim ordered", n, fill_cached_(n), fill_cached_is_lb_(n), me);
            if(fill_cached_is_lb_(n)){ untested();
                assert(me >= fill_cached_(n));
            }else{
                assert(me == fill_cached_(n));
            }

				bool ok=false;
				auto cr = boost::adjacent_vertices(n, _supergraph);
				for(; cr.first!=cr.second; ++cr.first){
					auto nn = *cr.first;
					if(nn==c){
						ok = true;
					}else{
					}
				}
				assert(ok);
        }
        //                assert(c_neigh.size()==degc);

#endif
		return c_neigh;
	}
#endif
	void debug_fill2(unsigned degc, std::set<vertex_descriptor> const& c_neigh, vertex_descriptor c) {
//		typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
#if defined(DO_TRACE) || !defined(NDEBUG)
		bool wrong=false;
#endif
//		marker_type debug_marker(_num_vert);
		for(auto n : c_neigh){
			if(supernode_size(n)<=0){
				// --degc;
				// hmm anything left to check?
				continue;
			}else{
			}

#ifndef NDEBUG

			// n is connected to c exactly once.
			auto q2 = boost::adjacent_vertices(n, _g);
			int k = 0;
			for(; q2.first!=q2.second; ++q2.first){
				auto m = *q2.first;
				assert(n != m);
				if(m==c){
					++k;
				}else{
				}
			}
			if(k==0){
				trace5("BUG missing link to c", c, n, k, c_neigh.size(), supernode_size(n));
			}else if(k==1){
			}else{ untested();
				trace5("BUG too many c", c, n, k, c_neigh.size(), supernode_size(n));
			}
			assert(k==1);
#else
			(void)c;
#endif


			if(_supergraph._marker.is_done(n)){ untested();
				trace1("DEBUG_FILL post-elim... done", n);
			}else{
				trace3("DEBUG_FILL post-elim...", n, fill_cached_(n), fill_cached_is_lb_(n));
#ifdef DO_TRACE
				auto me = _supergraph.count_missing_edges(n, &_debug_marker);
				trace3("DEBUG_FILL post-elim", n, me, fill_cached_(n));

				if(fill_cached_is_lb_(n)){ untested();
					assert(me >= fill_cached_(n));
				}else if(me == fill_cached_(n)){
				}else{ untested();
					trace3("DEBUG_FILL post-elim WRONG", n, me, fill_cached_(n));
					wrong=true;
				}
#endif
			}
		}
		trace2("DEBUG_FILL post-elim", degc, c_neigh.size());
		assert(!wrong);
		// assert(c_neigh.size()==degc); not really
	} // debug_fill2

public: // construct
	fillIn(G_t const &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
	  : baseclass(g, ub, ignore_isolated_vertices)
	  , _fill(baseclass::_supergraph, boost::num_vertices(g))
	  , _marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _fu(boost::num_vertices(g))
	  , _pick_marker(boost::num_vertices(g), get(boost::vertex_index, g))
#ifdef DEBUG_MARKER
	  , _debug_marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _dn_debug(boost::num_vertices(g), -1)
#endif
	  , _dn(boost::num_vertices(g), 0)
	  , _me2(boost::num_vertices(g), 0)
	{ untested();

//        boost::print_graph(g);
        treedec::check(g);
	}
	fillIn(G_t &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
	  : baseclass(g, ub, ignore_isolated_vertices)
	  , _fill(baseclass::_supergraph, boost::num_vertices(g), false)
	  , _marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _fu(boost::num_vertices(g))
	  , _pick_marker(boost::num_vertices(g), get(boost::vertex_index, g))
#ifdef DEBUG_MARKER
	  , _debug_marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _dn_debug(boost::num_vertices(g), -1)
#endif
	  , _dn(boost::num_vertices(g), 0)
	  , _me2(boost::num_vertices(g), 0)
	{
//        boost::print_graph(g);
		treedec::check(g);
	}

	fillIn(G_t &g, bool ignore_isolated_vertices, unsigned ub=-1u)
	  : baseclass(g, ub, ignore_isolated_vertices)
	  , _fill(baseclass::_subgraph, boost::num_vertices(g), false)
	  , _marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _fu(boost::num_vertices(g))
	  , _pick_marker(boost::num_vertices(g), get(boost::vertex_index, g))
#ifdef DEBUG_MARKER
	  , _debug_marker(boost::num_vertices(g), get(boost::vertex_index, g))
	  , _dn_debug(boost::num_vertices(g), -1)
#endif
	  , _dn(boost::num_vertices(g), 0)
	  , _me2(boost::num_vertices(g), 0)
	{ untested();
          // _cb(fill_update_cb(&_fill, baseclass::_subgraph))
//        boost::print_graph(g);
	}
private: // debugging
    size_t fill_cached_(vertex_descriptor v) const {
        return _fill.get_value(v);
    }
    bool fill_cached_is_lb_(vertex_descriptor v) const {
        return _fill.is_lb(v);
    }

public: // implementation
	using baseclass::_min;
	using baseclass::_o;
	using baseclass::_num_edges;
	using baseclass::_num_vert;
	using baseclass::_degree; // bug?
	using baseclass::_ub_tw;
	using baseclass::vertices_left;
	using baseclass::timer_on;
	using baseclass::timer_off;
	typedef typename baseclass::supergraph_type supergraph_type;
	typedef detail::predicate_collect_overlap<supergraph_type, marker_type, fill_type> predicate_collect_overlap;
	typedef detail::predicate_scan_neigh<supergraph_type, marker_type> predicate_scan_neigh;

	void do_it(){
		_numbering.init(1); // HACK.

		if(tree){
			_numbering.set_mode_tree();
		}else{
		}

		check(_g);
		_fill.init(_g);

		timer_on();

		if(!_num_vert){ untested();
			timer_off();
			return;
		}else{
		}

		assert(_o);
//        O_t& o = *_o;

#ifndef NDEBUG
		check(_g);
#endif

		baseclass::initialize();
		{
			// used in initialise
			_supergraph._marker.increment_tag();// ?
			_supergraph._marker.assert_clear();// ?
		}
		_o->resize(_num_vert);
		assert(_o->size() == _num_vert);
		vertex_descriptor c = 0;

		while(next(c)){
				assert(supernode_size(c)>0);
				assert(!_supergraph.is_numbered(c));

				int f = fill_cached_(c);
				assert(f>=0);
				if(unsigned(f)>_max_fill){
					_max_fill = f;
				}else{
				}

#ifdef DEBUG_FILL
				_supergraph._marker.assert_clear();// no. N(c) is marked (why?)
				trace2("DEBUG_FILL pre elim", c, fill_cached_(c));
				auto me = _supergraph.count_missing_edges(c, &_debug_marker);
				trace3("DEBUG_FILL c elim", fill_cached_(c), fill_cached_is_lb_(c), me);
				assert( me == fill_cached_(c) );
				auto c_neigh = debug_fill(c);
				_supergraph._marker.increment_tag();// ?

				auto ww = boost::adjacent_vertices(c, _supergraph);
				for(; ww.first!=ww.second; ++ww.first) {
					auto n = *ww.first;
					trace3("check: prenumber neigh", c, n, _fill.get_value(n));
					assert(!_supergraph.is_numbered(n));
				}
#endif

				_numbering.put(c); // here??

				if(_numbering.all_done(_supergraph.supernode_size(c))){ untested();
					_numbering.increment(_supergraph.supernode_size(c));
					break;
				}else{
				}

				_supergraph._marker.assert_clear();// no. N(c) is marked (why?)

				eliminate(c);

				auto degc = boost::out_degree(c, _g);
				auto sns = supernode_size(c);
				assert(sns);

				// if(c==5844){
				// 	auto vv = boost::adjacent_vertices(c, *_g);
				// 	for(; vv.first!=vv.second; ++vv.first) {
				// 		std::cout <<"c " << *vv.first <<" " << _numbering.get_position(*vv.first)  <<"\n";
				// 	}
				// }else{
				// }


				_marker.increment_tag();
				_marker.assert_clear();
				assert(sns>0);
				size_t real_degc = sns-1;
            auto bb = boost::adjacent_vertices(c, _g);
			  // 	_supergraph.bag_vertices(c); not initialised yet.
				for(; bb.first!=bb.second; ++bb.first){
					auto i = *bb.first;
					if(_marker.is_done(i)){ untested();
						// unreachable(); // no. why?
					}else if(_marker.is_tagged(i)){ untested();
					}else if(tree){
						if(supernode_size(i) < 0){ untested();
							// real_degc += 1 - supernode_size(i);
						}else{
							real_degc += supernode_size(i);
						}
						_marker.mark_tagged(i);
					}else{ untested();
						//incomplete();
						real_degc += supernode_size(i);
						_marker.mark_tagged(i);
					}
				}
				trace4("gone", c, _numbering.get_position(c), degc, real_degc);

				CFG::message(bLOG, "next %d. degree %d, fill %d, sns %d\n", c, real_degc, f, sns);

				if(real_degc > _ub_tw){
					_ub_tw = real_degc;
				}else{
				}

				if (_numbering.all_done()){
					break;
				}else{
					// tidy up marker?
					_supergraph._marker.set_tag_as_multiple_tag(); // catch up.
					_supergraph._marker.increment_tag();           // drop
					_supergraph._marker.assert_clear();

#ifdef Q_SHIFT
					assert(false)
						q_fill_update(c);
#endif


#ifdef DEBUG_FILL
					_supergraph._marker.clear();
					debug_fill2(degc, c_neigh, c);
#endif
					//--cnt;
					//assert(cnt == vertices_left());

					++_i;
				}
		} // next loop

		// BUG
		trace2("done loop", _i, _num_vert);

        // add one more node maybe.
        // should list remaining clique, not implemented
		postprocessing();

		timer_off();

		CFG::message(bDEBUG, "done. max fill=%d\n", _max_fill);

//		  boost::print_graph(_g);

	} // do_it

	using baseclass::_numbering;
    // fillIn::
	bool next(typename baseclass::vertex_descriptor &c){
		trace3("next loop", c, _num_edges, vertices_left());
		if(_numbering.all_done()){
			return false;
		}else{
			// todo: what do we know about lower bound?
			auto p = _fill.pick_min(0, -1u, true);
			_supergraph._marker.increment_tag(); // may have used cme.
			_supergraph._marker.assert_clear();

			c = p.first;
			trace2("next picked", c, p.second);

			if(_min<0){ untested();
				// 
			}else{
			}

			_min = p.second; // the fill of c.
			return true;
		}
	}

	using baseclass::_g;
	using baseclass::_degreemap;
	using baseclass::_subgraph;
	using baseclass::_supergraph;
	typedef typename supergraph_type::marker_type sg_marker_type;

	// update a neighbor of c.
	void update_fill_n(typename baseclass::vertex_descriptor n ) {
		incomplete();
		auto& marker = _supergraph._marker;
		assert(!_numbering.is_numbered(n));
		//            assert(marker.is_marked(n)); not yet.
		// assert(fill_cached_is_lb_(n) || marker.is_done(n));

#if 0
        // also correct? probably not
        auto mm = boost::adjacent_vertices(n, _g);
#else
		auto mm = boost::adjacent_vertices(n, _supergraph);
#endif
		for(; mm.first!=mm.second; ++mm.first){ untested();
			auto m = *mm.first;
			trace2("qeval 2?", n, m);
			//            assert(m!=c);
			if(marker.is_done(m)){ untested();
				//mass elim?
			}else{ untested();
			}
			if(marker.is_marked(m)){ untested();
			}else if(_numbering.is_numbered(m)){ untested();
				// updated from the other side if needed?
			}else{ untested();
				assert(m!=n);
				marker.mark(m);
				_fill.shift(m, minlong );
				_fill.q_eval(m);
				trace2("qeval 2", n, m);
				assert(fill_cached_is_lb_(m));
			}
		}

	}
#ifdef Q_SHIFT
	void q_fill_update(typename baseclass::vertex_descriptor c ) { untested();
        auto& marker = _supergraph._marker;
        trace1("q_fill_update", c);
        assert(_numbering.is_numbered(c));

        auto nn = boost::adjacent_vertices(c, _g);
        marker.clear();
        marker.mark(c);
        marker.set_multiple_tag(1);

        for(; nn.first!=nn.second; ++nn.first){ untested();
            auto n = *nn.first;
            trace2("q_fill_update", c, n);
            if(marker.is_done(n)){ untested();
                //mass elim?
            }else{ untested();
                assert(!_numbering.is_numbered(n));
            }

            if(marker.is_marked(n)){ untested();
                trace1("no q_fill_update", n);
            }else{ untested();
                marker.mark(n);
                _fill.shift(n, minlong );
                _fill.q_eval(n);
                // _fill.prefer(n);
                trace2("qeval1", n, boost::out_degree(n, _g));
            }

            if (marker.is_multiple_tagged(n)){ untested();
            }else{ untested();
                update_fill_n(n);
                marker.mark_multiple_tagged(n);
            }
        }
#if 0
        nn = boost::adjacent_vertices(c, _g);
        for(; nn.first!=nn.second; ++nn.first){ untested();
            auto n = *nn.first;
            if(marker.is_done(n)){ untested();
            }else{ untested();
                // update_fill_n(n);
                // _fill.shift(n, minlong );
                // _fill.q_eval(n);
                // _fill.prefer(n);
                trace1("qeval2", n);
            }
        }
#endif
        marker.set_tag_as_multiple_tag();
	} // q_fill_update
#endif

	void cleanup_mark_fill(vertex_descriptor c) {
		auto element_neighbor = _supergraph._work_space.make_stack();

		// Create two function objects for edge removal
		typedef typename supergraph_type::Workspace::stack WorkStack;
		detail::remove_and_collect<supergraph_type, sg_marker_type,
		    typename supergraph_type::numbering_type, WorkStack,
		    typename supergraph_type::VertexIndexMap>
		    p(_supergraph, _supergraph._marker, _supergraph.numbering(),
		      element_neighbor, get(boost::vertex_index, _g));

		// Reconstruct the adjacent node list, push element neighbor in a
		// List.
		trace1("rm and collect", c);
		remove_out_edge_if(c, p, _g);
		// during removal element neighbors are collected.
#ifndef NDEBUG
		{
			auto ii = boost::adjacent_vertices(c, _g);
			for(; ii.first!=ii.second; ++ii.first) {
				vertex_descriptor i = *ii.first;
				trace2("cleanup_mark first pass", c, i);
				assert(supernode_size(i) > 0);
			}
		}
#endif

			while (!element_neighbor.empty()) {
				// element absorb
				size_t e_id = element_neighbor.top();
				//            vertex_descriptor element = get(_index_vertex_map, e_id);
				auto ii = boost::adjacent_vertices(e_id, _g);
				trace2("cleanup_mark", c, e_id);
				for(; ii.first!=ii.second; ++ii.first) {
					vertex_descriptor n = *ii.first;
					if (_supergraph._marker.is_tagged(n)){
					}else if( tree && supernode_size(n)<0 ){
						// child node
						trace4("element edge0", c, n, supernode_size(c), supernode_size(n));
						_supergraph._marker.mark_tagged(n);
#ifndef NDEBUG
						auto dd = boost::out_degree(n, *_g);
#endif
//						boost::add_edge(c, n, _g); // ???
						assert(dd == boost::out_degree(n, *_g));
					}else if( !_numbering.is_numbered(n)) {
						_supergraph._marker.mark_tagged(n);

						// add node to newly created clique
						trace4("element edge1", c, n, supernode_size(c), supernode_size(n));
#ifndef NDEBUG
						auto dd = boost::out_degree(n, *_g);
#endif
						boost::add_edge(c, n, _g);
						assert(dd == boost::out_degree(n, *_g));
					}else{
					}
				}
				element_neighbor.pop();
			}
	} // cleanup_mark_fill

	 // fillcollect missing forward edges. return number of existing forward edges.
	template<class VI>
	void fillcollect(VI m, VI end, vertex_descriptor c, vertex_descriptor n,
	                         vertices_size_type DN, fill_value_t& missing_edges_left,
	                         bool only_cliques) {
		int k=0;
		
		predicate_collect_overlap p3(_supergraph, _marker, _fill, c, n);

		int num_edges_filled = 0;
		for(; (missing_edges_left||me3) && m!=end; ++m){
			p3.reset();

#ifndef NDEBUG
			{
				auto vv = boost::adjacent_vertices(n, *_g);
				for(; vv.first!=vv.second; ++vv.first) {
					assert( n != *vv.first);
				}
			}
			{
				auto vv = boost::adjacent_vertices(c, *_g);
				for(; vv.first!=vv.second; ++vv.first) {
					auto n = *vv.first;
					assert(_supergraph._marker.is_extra(n));
				}
			}
#endif

			if (_supergraph._marker.is_multiple_tagged(*m)){
				trace3("no fillcollect", c, n, *m);
				continue; // not a missing edge to n.
			}else if(czo && only_cliques){
				CFG::message(bLOG, "czo: missing, DN %d\n", DN);

				// _fu[*m] += (DNm - overlap) * _supergraph.supernode_size(n);
				_dn[*m] -= _supergraph.supernode_size(n);

				// _fu[n] += (DN - overlap) * _supergraph.supernode_size(*m);
				_fu[n] += DN * _supergraph.supernode_size(*m);

				auto wt = supernode_size(*m) * supernode_size(n);
				num_edges_filled += wt;
				missing_edges_left -= wt;
				continue;
			}else{
				// will fill edg n--m
				auto wt = supernode_size(*m) * supernode_size(n);
				num_edges_filled += wt;
				trace3("fillcollect overlap nonedge", c, n, *m);
				missing_edges_left -= wt;
			}
			assert(*m!=c);

			//multiple tag is relative to tag.
			// just tag, forget about multiple?
			_marker.set_multiple_tag(++k); // drop count marks from multiple-tagged N(c)
			
			trace3("fillcollect", c, n, *m);
			trace1("fill marker set mt", k);
			assert(n!=*m);
#ifndef NDEBUG
			std::set<unsigned> K;
			auto posc = _numbering.get_position(c);
#endif

//			_marker.clear(); // fill marker?

			assert(_marker.is_done(c));

			// not yet.
			_numbering.unput(c); // BUG
			assert(!_supergraph.is_numbered(c));

			boost::remove_out_edge_if(*m, p3, *_g);

			trace4("done p3", c, n, *m, p3.cnt());
			vertices_size_type overlap = p3.overlap();
			if(only_cliques){
//				assert(!overlap);
			}else{
			}

			_numbering.put(c);
			assert(_supergraph.is_numbered(c));
			assert( posc == _numbering.get_position(c));

#ifndef Q_SHIFT
			auto DNm = p3.cnt();
#ifndef NDEBUG
			// assert(_dn_debug[*m] == -1);
			_dn_debug[*m] = DNm;
#endif

			trace2("check: p3 says", *m, DNm);
//			trace4("fillcollect sum overlap", c, n, *m, ncoverlap);
			trace6("fillcollect sum overlap", c, n, *m, overlap, DNm, DN);

			assert(overlap <= DN);
			assert(overlap <= DNm);
			_fu[n] += (DN - overlap) * _supergraph.supernode_size(*m);

			if(me3 && ( overlap==DNm || overlap==DN )){
				CFG::message(bLOG, "outmatch? %d -> %d dn=overlap %d\n", *m, n, DN);

				_me2[*m] += *m - n;
				_me2[n] += n - *m;
				std::swap(_me2[n], _me2[*m]);
			}else if(me2 && overlap==DNm && DNm==DN) {

				CFG::message(bTRACE, "found me2 m=%d -> n=%d. was: %d\n", *m, n, _me2[*m]);
				trace3("me2: found", c, *m, n);

				_me2[*m] += *m - n;
				_me2[n] += n - *m;
				std::swap(_me2[n], _me2[*m]);
			}else{
				_fu[*m] += (DNm - overlap) * _supergraph.supernode_size(n);
				if(me2){
					assert(_supergraph.supernode_size(n) >= 0);
				}else{
					assert(_supergraph.supernode_size(n) == 1);
				}
			}

#else
			incomplete();
			assert(false);
#endif
		} // m loop, n<m

		trace4("fillcollect shift no disconnect", DN, c, n, _fu[n]);
//		_fu[n] -= DN * _supergraph.supernode_size(c);
		_marker.set_tag_as_multiple_tag();
	} // fillcollect

	void sort_out_me2 (vertex_descriptor c) {
		assert(me2);
		auto vv = boost::adjacent_vertices(c, *_g);
		for(; vv.first!=vv.second; ++vv.first) {
			auto N = *vv.first;
			if(_marker.is_done(N)){
				// assert(!boost::out_degree(N, *_g)); // no.
			}else if(_me2[N]){
				assert(boost::out_degree(N, *_g));
				auto n = N;
//				boost::add_edge(c, n, *_g); // **
				n += _me2[N];
				_me2[N] = 0;
				push_front_edge(N, c, *_g);
				while(_me2[n]){
					if(_marker.is_done(n)){ untested();
						// BUG: still need it in the bag?
						boost::add_edge(c, n, *_g);
					}else{
						_marker.mark_done(n);
						CFG::message(bLOG, "elim %d: merge %d into %d\n", c, n, N);

						assert(_supergraph._supernode_size[n]>=0);
						assert(_supergraph._supernode_size[N]);

						merge_vertices(n, N);
//						boost::add_edge(N, n, *_g); // **

						assert(_numbering.is_numbered(n));

						_marker.mark_done(n);
						_fill.remove(n);

					}

					auto nn = _me2[n];
					_me2[n] = 0;
					n += nn;
				}
			}else{
				assert(boost::out_degree(N, *_g));
				push_front_edge(N, c, *_g);
			}
		}

#ifndef NDEBUG
		for(auto i:_me2){
			assert(!i);
		}
#endif
	}

	// baseclass?
	void merge_vertices(vertex_descriptor n, vertex_descriptor N){
		trace2("merge", n, N);
		trace2("merge", _supergraph.supernode_size(N), _supergraph.supernode_size(n));

		_numbering.indistinguishable(n, N); // _data[n] = N
		assert(_numbering.is_numbered(n));

		// BUG _supernode_size is public.
		// TODO: supergraph.merge_vertices(n, N);
		if(!tree){
			// old allocator. does not work with me2
			_supergraph._supernode_size[N] += _supergraph._supernode_size[n];
			_supergraph._supernode_size[n] = 0;
		}else if(_supergraph._supernode_size[n]>0){
			_supergraph._supernode_size[N] += _supergraph._supernode_size[n];
			_supergraph._supernode_size[n] = 1 - _supergraph._supernode_size[n];
		}else{
			assert(_supergraph._supernode_size[n] < 0);
			// it has children, collect them
			_supergraph._supernode_size[N] += 1 - _supergraph._supernode_size[n];
		}
		trace3("merged", tree, _supergraph._supernode_size[N], _supergraph._supernode_size[n]);

		// TODO: mark_done(n)?
	}


	void eliminate(vertex_descriptor c) {
		//_marker.assert_clear();
#ifndef NDEBUG
		for(auto i:_me2){
			assert(!i);
		}
		{ itested();
			auto vv = boost::adjacent_vertices(c, *_g);
			for(; vv.first!=vv.second; ++vv.first) { itested();
				auto n = *vv.first;
				trace3("check: eliminate (clique) neigh", c, n, _fill.get_value(n));
			}
		}
		std::set<vertex_descriptor> K;
#endif

		fill_value_t missing_edges_left = _fill.get_value(c);
		_supergraph._marker.assert_clear();
		_supergraph._marker.mark(c);

		// mark all neighs of c.
		// add edges to them.
		// build NC==adjacent_vertices(c)
		// todo: tag N(c) with _marker?
		cleanup_mark_fill(c);

		vertices_size_type degc = boost::out_degree(c, *_g);

		if(1) { // raise.
			_supergraph._marker.set_extra_tag(1+degc); // sth N(c) tags survive increments.
			auto vv = boost::adjacent_vertices(c, *_g);
			for(; vv.first!=vv.second; ++vv.first) {
				auto n = *vv.first;

				trace3("mark extra", c, n, supernode_size(n));
				assert(!_supergraph._marker.is_done(n));
				_supergraph._marker.mark_extra(n);
				assert(_supergraph._marker.is_tagged(n));
#ifndef NDEBUG
				K.insert(n);
				_dn_debug[n] = -1;
#endif
			}
		}
#ifndef NDEBUG
		assert(K.size() == degc);
#endif

		size_t DN = 0;
		predicate_scan_neigh p2(_supergraph, _marker, c);

		int lvl = 0;
		size_t me_check = 0;
//        std::vector<int> fu(boost::out_degree(c, *_g));

		trace2("update", c, boost::out_degree(c, *_g));
		auto vv = boost::adjacent_vertices(c, *_g);

		// N(c) loop
		//
		for(; vv.first!=vv.second; ) {
			_supergraph._marker.increment_tag();          // lvl < degc
			_supergraph._marker.assert_extra();
			_supergraph._marker.set_multiple_tag(1+degc);   // lvl + degc
			++lvl;

			_marker.increment_tag();
			_marker.assert_clear();
			_marker.set_multiple_tag(1);

#ifndef NDEBUG
			auto dbg_first = vv.first;
			auto dbg_first_pre = vv.first;
#endif

			vertex_descriptor n = *vv.first;
			++vv.first;

			// _marker.clear();

			// update out edges of n
			trace4("update call p2", c, n, boost::out_degree(n, *_g), lvl);
			_marker.mark(n); // for clearing connections from n // extra?
			_marker.assert_mclear(); // prepare p2 call.

			// if (!degree_lists_marker.need_update(v_node) 
			//         && !degree_lists_marker.outmatched_or_done(v_node)) { untested();
			//     degreelists.remove(v_node);
			// }
			// c is also marked (reinsert below **?)
			// assert(_supergraph._marker.is_marked(c));
			assert(_supergraph.is_numbered(c));

			/// untagged -> tag (N(n) \ N(c))
			/// done -> delete?
			/// extra -> mt  (N(n) intersect N(c)) & delete.
			p2.reset();
#ifndef NDEBUG
			for(; dbg_first_pre!=vv.second; ++dbg_first_pre) { itested();
				vertex_descriptor mp = *dbg_first_pre;
				trace3("check: mp>n pre p2", c, n, mp);
			}
#endif


			auto twin=false;
			if(_fill.get_value(n)){
				// has relevant neighbours
				trace3("p2 call", n, boost::out_degree(n,_g), boost::out_degree(n,*_g));
				boost::remove_out_edge_if(n, p2, *_g);
				DN = p2.cnt();
				trace3("check: p2 says", n, DN, p2.only_cliques());

//				if(p2.rcl() > 4* p2.rcnt()){
//					CFG::message(bLOG, "scan %d: deg %d cnt %d rcnt %d same %d, clique sum %d\n", n, boost::out_degree(n,_g),
//							p2.cnt(), p2.rcnt(), p2.cnt()==p2.rcnt(), p2.rcl());
//				}

#ifndef NDEBUG
				assert(_dn_debug[n] == int(DN) || _dn_debug[n] == -1);
				_dn_debug[n] = DN;
#endif
				trace3("p2 done", c, n, DN);

#ifdef DEREF
				if(deref){
					for(auto v: p2._stack_hack){
						trace3("stackhack", c, n, v);
						//boost::add_edge(n, v, *_g);
						assert(!_supergraph.is_numbered(v));
						assert(_supergraph.supernode_size(v)>0);
						assert(v!=c);
						assert(!(boost::edge(n, v, _g).second));
#ifndef NDEBUG
						auto dn=boost::out_degree(n, _g);
						auto dv=boost::out_degree(v, _g);
#endif
#ifdef DEBUG__
						auto vv = boost::adjacent_vertices(n, *_g);
						for(; vv.first!=vv.second; ++vv.first) {
							// ...
						}
#endif

						incomplete(); // needs gala, and does not help a lot.
						// _g->out_edges(n).push_back(v);
						assert((1+dn)==boost::out_degree(n, _g));
						assert(dv==boost::out_degree(v, _g));
						assert((boost::edge(n, v, _g).second));
					}
				}else
				
#endif // DEREF
				{ itested();
				}
			}else{
				// fill_value(n) == 0
				// like p2, but faster
				// (todo: omit)
//				CFG::message(bTRACE, "nc TWIN again? %d -> %d\n", n, c);
				boost::clear_out_edges(n, *_g);
				twin = true;
				DN = 0;
				_dn[n] = 0;

				auto vv = boost::adjacent_vertices(c, *_g);
				for(; vv.first!=vv.second; ++vv.first) {
					auto n = *vv.first;
					trace2("check: N(c)", c, n);
					if(_supergraph._marker.is_done(n)){
					}else{
					  	assert(_supergraph._marker.is_extra(n));
						_supergraph._marker.mark_multiple_tagged(n);
					}
				}
			}

			_dn[n] = DN;

#ifdef DO_TRACE
			{
				auto vv = boost::adjacent_vertices(c, *_g);
				for(; vv.first!=vv.second; ++vv.first) {
					auto n = *vv.first;
					trace3("check: list N(c)", c, n, supernode_size(n));
				}
			}
#endif

#ifndef NDEBUG
			++dbg_first;
			// neighbours > n
			for(; dbg_first!=vv.second; ++dbg_first) { itested();
				vertex_descriptor mev = *dbg_first;
				if (_supergraph._marker.is_done(mev)){ itested();
				}else if (twin || _supergraph._marker.is_multiple_tagged(mev)){ itested();
					trace3("check: existing edg", c, n, mev);
				}else if(me2){
					assert(supernode_size(n)>0);
					assert(supernode_size(mev)>0);
					me_check += supernode_size(n) * supernode_size(mev);
					assert(!_supergraph._marker.is_done(mev));
					assert(!_supergraph._marker.is_done(n));
				}else{
					trace3("check: new edg", c, n, mev);
					assert(!_supergraph._marker.is_done(mev));
					assert(!_supergraph._marker.is_done(n));
					// stack.push(mev) // needed?
					++me_check;
				}
			}
#else
			(void) twin;
#endif

			trace3("done p2", c, n, boost::out_degree(n, *_g));

#if 0 // not yet?
			if(_supergraph.supernode_size(n) <= 0){
			}else if (boost::out_degree(n, *_g) == 0) {
				// n was only connected to elements in N(c)?
			}else if (_dn[n] == 1 && boost::out_degree(n, *_g) == 1) {
				CFG::message(bLOG, "singleton left %d %d\n", n, c);
				assert(_dn[n] == 1);
#ifndef NDEBUG
				auto vv = boost::adjacent_vertices(n, *_g);
				assert(!deref || !_numbering.is_numbered(*vv.first));
#endif

			}else if (boost::out_degree(n, *_g) == 1) {
				auto vv = boost::adjacent_vertices(n, *_g);
				auto b = *vv.first;
				assert(_dn[n] > 1);
				auto degb = boost::out_degree(b, *_g);
				auto sns = supernode_size(b);
				trace5("one clique", n, b, degb, _dn[n], sns);
				assert(n!=b);
				if(sns == _dn[n]){
				}else{
					assert(_numbering.is_numbered(b));
					assert(c!=b);
					CFG::message(bLOG, "one clique %d (%d) size: %d %d\n", n, c, degb, _dn[n]);
					if(!flb){
					}else if( _dn[n]*4<degb){
						cherry_pick(n, _dn[n], b);
					}else{
					}
				}
			}else{
			}
#endif

			assert(_supergraph.supernode_size(n) > 0);

			if (::me && boost::out_degree(n, _g) == 0) {
				CFG::message(bLOG, "merge %d (%d) into %d (%d)\n", n, supernode_size(n), c, supernode_size(c));
#ifndef NDEBUG
				auto dd = boost::out_degree(c, *_g);
#endif
//				boost::add_edge(n, c, *_g);
				assert( dd == boost::out_degree(c, *_g));
				
				merge_vertices(n, c);

				_supergraph._marker.mark_done(n);
				_marker.mark_done(n); // avoid in sort_out_me2
				_fill.remove(n);
#ifndef NDEBUG
				_debug_marker.mark_done(n);
#endif
				assert(_numbering.is_numbered(n));
				// needed in sg2tree..?
				trace4("merged", n, c, _supergraph.supernode_size(c), _supergraph.supernode_size(n));

#ifndef NDEBUG
				{
					int k=0;
					auto vv = boost::adjacent_vertices(c, *_g);
					for(; vv.first!=vv.second; ++vv.first) {
						auto x = *vv.first;
						k += (x==n);
					}
					assert(k);
				}
#endif
			}else{
				// not indistinguishable nodes
				/// visit all m > n
				//// NC(n) is marked by marker with mt=degc
				_marker.mark_done(c); // BUG?

				bool only_cliques = p2.only_cliques();
				trace3("call fillcollect", c, n, only_cliques);


				fillcollect(vv.first, vv.second, c, n, DN, missing_edges_left, only_cliques);

				_marker.mark(c); // not done, because numbered.
				// todo? mark once all neighbours expire

				// just store in _fill?
//				unsigned missing_edg_n = degc - 1 - ne - nef;
				assert(me3 || (missing_edges_left >= 0));

				// done with m<=n. add edge to center
				if(!me2){
					// connect them all.
					boost::add_edge(n, c, *_g); // **
				}else{
					// connect in sort_out_me2
				}

				// incomplete();
				// _degree_lists_marker.mark_need_update(v_node);
			}

			// _marker.increment_tag(); // forget about N(n)
			// _supergraph._marker.set_tag_as_multiple_tag(); // probably not?
			//
			trace0("set tag as multiple");
			_marker.set_tag_as_multiple_tag();
			trace1("fill mark done 2", n);
		} // first-neigh loop.

		if(me2){
			sort_out_me2(c);
		}else{
		}

		_marker.increment_tag(); // forget about N(n)

		{ // fillflush
			auto vv = boost::adjacent_vertices(c, *_g);
			for(; vv.first!=vv.second; ++vv.first) {
				auto n = *vv.first;
				trace4("fillflush", n, _fu[n], _fill.get_value(n), _supergraph.supernode_size(n));
//				if (boost::out_degree(n, _g)) { untested(); }
				if (_supergraph._marker.is_done(n)){
					_fu[n] = -1;
				}else if(supernode_size(n)<=0){
					// ?
					_fu[n] = -1;
				}else{
					assert(_me2[n] == 0);
					assert(boost::out_degree(n, _g)); // ==1?
					//					  boost::add_edge(n, c, *_g); // **
					_marker.mark(n); // untag done.

					assert(-long(_fu[n]) <= long(_fill.get_value(n)));
					trace3("fillflush disconnect", n, _fu[n], _fill.get_value(n));
					assert(_supergraph.supernode_size(c)>0);
					_fu[n] -= _dn[n] * _supergraph.supernode_size(c);

					_fill.shift(n, _fu[n]);
					_fu[n] = 0;
				}

#ifndef NDEBUG
#endif
				_dn[n] = 0;
			}
		} // fillflush

		auto fillc = _fill.get_value(c);
		_numbering.increment(_supergraph.supernode_size(c));

		trace2("increment", c, _supergraph.supernode_size(c));
		trace3("increment", c, me_check, fillc);
		assert(_numbering.is_numbered(c));

		if ( me_check == fillc ){
		}else{
//			std::cerr << "mismatch in " << c << ": " << me_check << " " << fillc << "\n";
		}
		assert(me3 || ( me_check == fillc )); // really?

		_supergraph._marker.set_tag_as_multiple_tag();

		{ // is this needed?
			_supergraph._marker.increment_tag();
			_supergraph._marker.assert_clear();
		}

	} // eliminate(c)

	using baseclass::_i;
	void postprocessing(){
		trace2("post", _i, baseclass::_num_vert);
		auto const& sns = _supergraph.supernode_sizes();
		for(unsigned i=0; i<sns.size(); ++i){
			trace2("post", i, sns[i]);
		}
		if(!baseclass::_iiv){
			assert(boost::num_vertices(_g) == baseclass::_num_vert);
		}else{ untested();
		}
		if(_i == baseclass::_num_vert){
			// no nodes or so.
			unreachable(); //?
		}else if(0){ untested();

			auto w = _fill.pick_min(0, 0, true).first;
			trace2("post", w, baseclass::_degreemap[w]);
			(*_o)[baseclass::_i++] = w;
			baseclass::_numbering.put(w);
			baseclass::_numbering.increment();
			auto x = _subgraph.adjacent_vertices(w);
			for(;x.first!=x.second;++x.first){ untested();
				trace1("last node neigh", *x.first);
			}

            for(; _i < baseclass::_num_vert; ++_i){
                auto v = _fill.pick_min(0, 0, true).first;
//                treedec::add_edge(w, v, baseclass::_g); already there.
                assert(_i < _o->size());
                (*_o)[baseclass::_i] = v;
            }
		}else{
			assert(_o);
			assert(_o->size() == boost::num_vertices(_g));
			std::vector<long> oo(_o->size()); // must be signed.
			_numbering.get_ordering(oo, _supergraph._supernode_size);
			*_o = O_t(oo.begin(), oo.end());
			baseclass::_i = _o->size(); // length of the ordering used in baseclass...
		}

		if(!baseclass::_iiv){
			//assert(baseclass::_i == _o->size());
			//assert(baseclass::_i == baseclass::_num_vert);
			//assert(baseclass::_i == boost::num_vertices(_g));

			// assert(baseclass::_i == _numbering.total()); for some reason, numbering numbers the bags only.
		}else{ untested();
		}
		trace2("post done", _i, baseclass::_num_vert);
		for(unsigned i=0; i<sns.size(); ++i){
			trace2("post", i, sns[i]);
		}
	}

	template<class T>
	void get_tree_decomposition(T& t){
		assert(_o);
		unsigned i=0;
		for(auto oi : *_o){
			trace3("dbg", i, oi, _numbering.get_position(i));
			++i;
		}
		assert(treedec::is_vertex_permutation(*_o, _g));
		auto const& sns = _supergraph._supernode_size;
		for(unsigned i=0; i<sns.size(); ++i){
			trace2("gtd", i, sns[i]);
		}

		auto bs = 1+_ub_tw;
#if 0
        tree_from_sg(_supergraph, *_o, t, bs, _supergraph._supernode_size);
#else
       // treedec::draft::inplace_bmdo_tree(_g, *_o, t, bs, _supergraph.numbering());
		 //
		 // calls inplace_bmdo for on suitable graphs..
		treedec::draft::tree_from_sg(_supergraph, *_o, t, bs, _supergraph.numbering());
#endif
    }

private: // debugging
    void check_vertex(vertex_descriptor c) { untested();
        size_t k=0;
        for(auto p = boost::adjacent_vertices(c, _subgraph); p.first!=p.second; ++p.first){ untested();
            ++k;
        }
        assert(k==baseclass::_degreemap[c]);

        auto p=boost::adjacent_vertices(c, _subgraph);
        for(; p.first!=p.second; ++p.first){ untested();
            auto n=*p.first;
            long degn = baseclass::_degreemap[n];
            trace4("adj", c, *p.first, degn, baseclass::_degreemap[c]);
            trace4("checking neigh", n, _fill.get_value(n), degn, _fill.is_lb(n));
            assert(2*_fill.get_value(n)<=size_t(degn*(degn-1)));

            auto q = boost::adjacent_vertices(n, _subgraph);
            for(; q.first!=q.second; ++q.first){ untested();
                auto n2=*q.first;
//                trace1("neigh", n2);
                --degn;
            }
            assert(!degn);
        }
        trace1("checked", c);
    }
	long /*vertices_size_type?*/ supernode_size(vertex_descriptor n) const{
		return _supergraph.supernode_size(n);
	}
	void cherry_pick(vertex_descriptor n, size_t howmany, vertex_descriptor clique){
		trace3("flb", howmany, n, clique);
		assert(clique!=n);
		assert(_supergraph.is_numbered(clique));
		_pick_marker.increment_tag();
		_pick_marker.mark(n);
		assert(boost::out_degree(n, _g) == 1);
		boost::clear_out_edges(n, _g);

		auto vv = boost::adjacent_vertices(clique, *_g);
		for(;vv.first != vv.second; ++vv.first){
			auto x = *vv.first;
			trace2("clique?", _supergraph._marker.is_multiple_tagged(x), _supergraph._marker.is_extra(x));
			trace2("clique?", x, _supergraph._marker.is_tagged(x));
			if(_pick_marker.is_marked(x)){
				trace1("pick marked", x);
			}else if(_supergraph.is_numbered(x)){
				trace1("numbered", x);
			}else if(!_supergraph._marker.is_multiple_tagged(x)){
				trace1("PICK", x);
				assert(howmany);
				boost::add_edge(n, x, *_g);
				_pick_marker.mark(x);
				 --howmany;
			}
		}

		trace1("flb", howmany);
		assert(howmany==0);
	}
private:
	unsigned _max_fill{0};
	fill_type _fill;
	marker_type _marker;
	std::vector<int> _fu;
	marker_type _pick_marker; // use stack instead?
#ifdef DEBUG_MARKER
	marker_type _debug_marker;
	std::vector<int> _dn_debug;
#endif
	std::vector<int> _dn;
	std::vector<int> _me2;
	std::vector<int> _tdata;
}; // fillIn

} // impl
} // pending
} // treedec

#endif // guard
