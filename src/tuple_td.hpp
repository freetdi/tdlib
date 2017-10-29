#pragma once

#include <boost/graph/graph_traits.hpp>
#include "treedec_traits.hpp"

// a list of tuples vertex_descriptor, set<vertex_descriptor>
// somehow encodes a tree decomposition.
// if (v, b) is such a tuple, then v \notin b.

namespace boost{

template<class U, class V>
unsigned add_vertex(
		std::vector<boost::tuples::tuple<U, std::set<V> > >& x
		)
{
	x.emplace_back();
	return x.size()-1;
}
template<class U, class V>
unsigned add_vertex(
		std::vector<boost::tuples::tuple<U, std::vector<V> > >& x
		)
{untested();
	x.emplace_back();
	return x.size()-1;
}

template<typename U, typename V>
struct graph_traits<
	std::vector<boost::tuples::tuple<U, std::set<V> > > >{
    //typedef typename T::vertex_property_type vertex_property_type;
		typedef unsigned vertex_descriptor;
};
template<typename U, typename V>
struct graph_traits<
	std::vector<boost::tuples::tuple<U, std::vector<V> > > >{
    //typedef typename T::vertex_property_type vertex_property_type;
		typedef unsigned vertex_descriptor;
};

template<class U, class V>
U& get(vertex_underlying_t,
	unsigned n,
	std::vector<boost::tuples::tuple<U, std::set<V> > >& bags)
{
	return get<0>(bags[n]);
}
template<class U, class V>
U& get(vertex_underlying_t,
	unsigned n,
	std::vector<boost::tuples::tuple<U, std::vector<V> > >& bags)
{
	return get<0>(bags[n]);
}

} //boost

namespace treedec{

template<class U, class V>
std::set<V>& bag( unsigned n,
		std::vector<boost::tuples::tuple<U, std::set<V> > >& x)
{
	return boost::get<1>(x[n]);
}
template<class U, class V>
std::vector<V>& bag(
		std::vector<boost::tuples::tuple<U, std::vector<V> > >& x,
		unsigned n)
{ untested();
	return boost::get<1>(x[n]);
}

template<typename U, typename V>
struct treedec_traits<
	std::vector<boost::tuples::tuple<U, std::set<V> > > >{
    //typedef typename T::vertex_property_type vertex_property_type;
		typedef unsigned vd_type;
		typedef std::set<V> bag_type;
};
template<typename U, typename V>
struct treedec_traits<
	std::vector<boost::tuples::tuple<U, std::vector<V> > > >{
    //typedef typename T::vertex_property_type vertex_property_type;
		typedef unsigned vd_type;
		typedef std::vector<V> bag_type;
};

} //treedec
