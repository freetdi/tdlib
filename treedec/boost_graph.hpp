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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

// Wrap something with a boost graph interface into a class
// This is unneeded, could wrap boost graphs directly.

#ifndef TREEDEC_BOOST_GRAPH_H
#define TREEDEC_BOOST_GRAPH_H

#include "boost_compat.h"
#include <boost/python.hpp>
#include "bits/any_iterator.hpp"
#include "exception.hpp"
#include <utility> // pair
#include <boost/variant/static_visitor.hpp>
#include "graph_iter.hpp"
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "graph.hpp"

// duplicate in frontend.hpp
namespace bits{

template<class I, class X=void>
struct fill_help{

    template<class vit>
    static void fill_idmap(vit i, vit e, I& idm, unsigned & /*max*/, unsigned&)
    { untested();
        for(;i!=e; ++i){ untested();
            idm.push_back(*i);
        }
    }

template <typename G_t, class IDM, class it, class it2>
static void fill(G_t &G, IDM const& idmap, it, it, it2 i2, it2 e2, size_t max)
{ untested();
    trace3("fill", max, idmap.size(), boost::num_vertices(G));
    typedef std::map<typename it::value_type, unsigned> idxmap_type;
    // hmm only unsigned?
    idxmap_type idxMap;
    auto bi=boost::vertices(G).first;
    // for(; ii!=ee; ++ii){ untested(); }
    for(auto ii=idmap.begin(); ii!=idmap.end(); ++ ii){ untested();
        idxMap[*ii] = *bi;
        ++bi;
        trace1("map", max);
    }

    for(; i2!=e2; ++i2) { untested();
        auto p=*i2;

        if(idxMap.find(p.first)==idxMap.end()){ untested();
            throw exception_invalid("source invalid");
            // + std::to_string(p.first) + ":" + std::to_string(p.first));
        }else if(idxMap.find(p.second)==idxMap.end()){ untested();
            throw exception_invalid("target invalid");
            // + std::to_string(p.second) + ":" + std::to_string(p.second));
        }else{ untested();
            treedec::add_edge(idxMap[p.first], idxMap[p.second], G);
        }
    }
}
};

// for "graph"
template<>
struct fill_help< boost::typed_identity_property_map<long unsigned int> > {
template<class vit, class I>
static void fill_idmap(vit i, vit e, I&, unsigned & max, unsigned& cnt)
{ untested();
    for(;i!=e; ++i){ untested();
        ++cnt;
//        _idmap.push_back(*i);
        if(unsigned(*i) > max){ untested();
            max = *i;
        }else{ untested();
        }
    }
}

template <typename G_t, class I, class it>
static void fill(G_t &G, I const&, it i2, it e2, it, it, size_t max)
{
    trace1("fill", max);

    for(; i2!=e2; ++i2) {
        auto p=*i2;

        if(unsigned(p.first)>=max){ untested();
            throw exception_invalid("source invalid "
             + std::to_string(p.first) + ":" + std::to_string(p.first));
        }else if(unsigned(p.second)>=max){ untested();
            throw exception_invalid("target invalid"
             + std::to_string(p.second) + ":" + std::to_string(p.second));
        }else{
            // takes care of reverse edges in directed graphs.
            treedec::add_edge(p.first, p.second, G);
        }
    }
}
};

//template<class idmap_t>
// "subgraph". essentially the same as Graph...
template<>
struct fill_help<std::vector<unsigned> > {
    //, std::enable_if< std::numeric_limits<idmap_t::value_type >::is_integer >::type > 
    template<class vit, class I>
    static void fill_idmap(vit i, vit e, I& _idmap, unsigned & max, unsigned& cnt) { untested();
        for(;i!=e; ++i){ untested();
            ++cnt;
            _idmap.push_back(*i);
            if(unsigned(*i) > max){ untested();
                max = *i;
            }else{ untested();
            }
        }
    }
template <typename G_t, class I, class it, class eit>
static void fill(G_t &G, I const& idmap, it, it, eit i2, eit e2, size_t max)
{ untested();
    trace3("fill", max, idmap.size(), boost::num_vertices(G));
    // hmm only unsigned?
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxmap_type;
    idxmap_type idxMap(max+1, -1u);
    auto bi=boost::vertices(G).first;
    for(auto iii=idmap.begin(); iii!=idmap.end(); ++iii){ untested();
        trace2("map", *iii, *bi);
        idxMap[*iii] = *bi;
        ++bi;
    }

    for(; i2!=e2; ++i2) { untested();
        auto p=*i2;

        if(idxMap.size() <= p.first || idxMap[p.first]==-1u){ untested();
            trace3("wrong?", p.first, p.second, max);
            throw exception_invalid("source invalid " + std::to_string(p.first)
                    + ":" + std::to_string(p.second));
        }else if(idxMap.size() <= p.first || idxMap[p.second]==-1u){ untested();
            throw exception_invalid("target invalid " + std::to_string(p.first)
                    + ":" + std::to_string(p.second));
        }else{ untested();
            treedec::add_edge(idxMap[p.first], idxMap[p.second], G);
        }
    }
}
};


} // bits;

namespace wrap {

using bits::fill_help;

template<class T>
struct vd_sel{
    typedef typename T::value_type type;
};

using IteratorTypeErasure::any_iterator;

// this is probably unneeded.
template<class backend_t, class idmap_type=std::vector<unsigned long> >
class graph{
public:
	typedef typename vd_sel<idmap_type>::type vertex_descriptor;
	typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;

	// backend iteration. private?!
	typedef any_iterator<const unsigned long,
			  std::input_iterator_tag,
			  const unsigned long> raw_vertex_iterator;

	typedef any_iterator<const edge_descriptor,
			  std::input_iterator_tag,
			  const edge_descriptor> edge_iterator;

private: // visitors
	struct get_verts;

	struct get_verts
		: boost::static_visitor<
		  std::pair<  boost::transform_iterator<get_verts, raw_vertex_iterator> ,
		  boost::transform_iterator<get_verts, raw_vertex_iterator> > > { //
		typedef boost::transform_iterator<get_verts, raw_vertex_iterator> vertex_iterator;

		get_verts(idmap_type const& m) : _m(m) {
		}
		get_verts(get_verts const& o) : _m(o._m) {
		}
		template<typename T>
		std::pair<vertex_iterator, vertex_iterator>
		operator()(T const& t) const {
			auto x=boost::vertices(t);
			auto b=raw_vertex_iterator(x.first);
			auto e=raw_vertex_iterator(x.second);

			vertex_iterator bt=make_transform_iterator(b, *this);
			vertex_iterator et=make_transform_iterator(e, *this);

			return std::make_pair(bt, et);
		}

		// hack. need proper idmaps, then pass
		// idmap to make_transform_iterator directly.
		typename idmap_type::value_type operator()(unsigned x) const {
//			trace2("map", x, _m.size());
//			if(x<_m.size()){ untested();
				return _m[x];
//			}else{
				// hack
//				return x;
//			}
		}

		idmap_type const& _m;
	};

	struct get_edgs {
		template<class G>
		struct edgpairmap{
			edgpairmap(idmap_type const& m, G const& g)
				: _m(m), _g(g)
			{
			}
			template<class E>
			std::pair<typename idmap_type::value_type,
				typename idmap_type::value_type>
					operator()(E edg) const {
						auto s=boost::source(edg, _g);
						auto t=boost::target(edg, _g);
						return std::make_pair(_m[s], _m[t]);
					}
			idmap_type const& _m;
			G const& _g;
		};
		//        typedef boost::transform_iterator<edgpairmap, raw_edge_iterator> edge_iterator;

		get_edgs(idmap_type const& m) : _m(m) {
		}
		get_edgs(get_edgs const& o) : _m(o._m) { untested();
		}

		template<typename G>
		std::pair<edge_iterator, edge_iterator>
		operator()(G const& g) const {
			typedef typename boost::graph_traits<G>::edge_iterator geit_t;
			auto x=boost::edges(g);
			geit_t b(x.first);
			geit_t e(x.second);

			// translation
			auto bt=make_transform_iterator(b, edgpairmap<G>(_m, g));
			auto et=make_transform_iterator(e, edgpairmap<G>(_m, g));

			// "cast" to any_iterator
			return std::make_pair(edge_iterator(bt),
					edge_iterator(et));
		}

		idmap_type const& _m;
		// G const& _g;
  };
public:
	// typedef boost::transform_iterator<idmap_type, raw_vertex_iterator> vertex_iterator;
	typedef typename get_verts::vertex_iterator vertex_iterator;
	//    typedef typename get_edgs::edge_iterator edge_iterator;
private:

	struct num_vert : boost::static_visitor<size_t> {
		template<typename T>
		size_t operator()(T const& t) const { untested();
			return boost::num_vertices(t);
		}
	};
	struct num_edg : boost::static_visitor<size_t> {
		template<typename T>
		size_t operator()(T const& t) const { untested();
			return treedec::num_edges(t);
		}
	};
public: // construct
	graph(std::string const&){ untested();
		incomplete();
	}
	graph(size_t, std::string const&){ untested();
		incomplete();
	}
	graph(size_t cnt=0, unsigned type=0) {
		(void) type;
		_g = backend_t(cnt);
	}
	template<class eit>
	graph( eit i2, eit e2, unsigned x, unsigned type=0) {
		trace2("creating graph from edgiter ", x, type);
		typedef fill_help<idmap_type> fh;
		_g = backend_t(x);
		// no. does not work on directed graphs
		// fill(boost::get<TD_graph_t>(_g), i, e, E, max);
		fh::fill(_g, _idmap, i2, e2, i2, e2, x);
	}

	// legacy stuff. only works with list<int>, list<tuple<int,int>>
	// no type checking...
	// this should create a graph with vertices mapped to what's in the first array
	// passed to it. perhaps the identifiers need to be unique...
	//graph(const boost::python::list& v,
	//         const boost::python::list& edg)
	// turnaround edges first.
	template<class vit, class eit>
	graph(eit i2, eit e2, vit i, vit e, unsigned type=0)
	{ untested();
		unsigned int max = 0;
		unsigned int cnt = 0;
		typedef fill_help<idmap_type> fh;
		fh::fill_idmap(i, e, _idmap, max, cnt);

		_g = backend_t(cnt);
		fh::fill(_g, _idmap, i, e, i2, e2, max);
	}

	graph(std::list<int>const& v, std::list<unsigned> const& e) { untested();
		(void)v;
		(void)e;
		incomplete();
	}
	graph(std::vector<unsigned> v, std::vector<unsigned> e) { untested();
		(void)v;
		(void)e;
		_g = backend_t();
	}
public:
	graph(const graph& o) : _g(o._g) { untested();
	}

public: //ops
	graph& operator=(const graph& o){ untested();
		// not used from within python.
		// not even accessible?!
		_g = o._g;
		return *this;
	}
public: // const access
	std::pair<vertex_iterator, vertex_iterator> vertices() const{
		return get_verts(_idmap)(_g);
	}
	std::pair<edge_iterator, edge_iterator> edges() const{
		return get_edgs(_idmap)(_g);
	}
	size_t num_vertices() const{ untested();
		return boost::num_vertices(_g);
	}
	size_t num_edges() const{ untested();
		return boost::num_edges(_g);
	}
	std::string backend_typename() const{ untested();
	   return typeid(backend_t).name();
	}
public: // modify
	void add_edge(vertex_descriptor, vertex_descriptor){
		incomplete();
	}
//	void add_edge(std::pair<vertex_descriptor, vertex_descriptor>){
//		incomplete();
//	}
public: // graph algorithms. maybe wrong place.
	template<class tdt>
	size_t preprocessing(tdt&){ incomplete(); return 0;}
	void min_degree();
	void exact_cutset();
	void exact_ta();
public: // hack
	template<class U>
	vertex_descriptor maphack(U x) const{ untested();
		assert(x<num_vertices());
		return _idmap[x];
	}
private:
	backend_t _g;
	idmap_type _idmap;
	// numbering?
	// order?
	// bags?? (use for undirected?)
};

} // wrap


/* ---- other stuff ----- */

using IteratorTypeErasure::any_iterator;

template<class X>
	using pair_ = std::pair<X, X>;

template<class vertex_descriptor>
using vertex_iterator_ = any_iterator< const vertex_descriptor,
			  std::input_iterator_tag, const vertex_descriptor>;

template<class G>
struct mapvertex{
	mapvertex(G const& g) : _g(g)
	{
	}
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	template<class V>
	unsigned long operator()(V v) const {
		return v;
	}

private:
	//idmap_type const& _m; // later
	G const& _g;
};

template<class G>
pair_< vertex_iterator_< typename boost::graph_traits<G>::vertex_descriptor> >
boost_vertices(G const& g)
{
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
//	typedef my_vd vertex_descriptor;
	typedef typename boost::graph_traits<G>::vertex_iterator vit_t;

	typedef any_iterator< const vertex_descriptor,
			  std::input_iterator_tag, const vertex_descriptor > vertex_iterator;

	auto x = boost::vertices(g);
	vit_t b(x.first);
	vit_t e(x.second);

	// translation
	auto bt = make_transform_iterator(b, mapvertex<G>(g));
	auto et = make_transform_iterator(e, mapvertex<G>(g));
	return std::make_pair(vertex_iterator(bt), vertex_iterator(et));
}


template<class G>
typename boost::graph_traits<G>::vertex_descriptor
boost_add_vertex(G& g)
{
	auto r = boost::add_vertex(g);
	return r;
}
#if 1
template<class G>
std::pair< pair_<typename boost::graph_traits<G>::vertex_descriptor>  , bool>
boost_add_edge(G& g, typename boost::graph_traits<G>::vertex_descriptor s,
                     typename boost::graph_traits<G>::vertex_descriptor t)
{
	auto r = boost::add_edge(s, t, g);
	auto s_ = boost::source(r.first, g);
	auto t_ = boost::target(r.first, g);
	return std::make_pair( std::make_pair(s_, t_), r.second);
}

#if 0
template<class G>
std::pair< pair_<typename boost::graph_traits<G>::vertex_descriptor>  , bool>
boost_add_edge(G& g,
		std::pair< typename boost::graph_traits<G>::vertex_descriptor,
	              typename boost::graph_traits<G>::vertex_descriptor> p)
{ untested();
	auto r = boost::add_edge(p.first, p.second, g);
	auto s = boost::source(r.first, g);
	auto t = boost::target(r.first, g);
	return std::make_pair( std::make_pair(s,t), r.second);
}
#endif

#else
template<class G, class V>
std::pair< std::pair<V, V>, bool> boost_add_edge(G& g, V s, V t)
{
	auto r = boost::add_edge(s, t, g);
	return r; // std::make_pair(0, r.second);
}
#endif
#endif
