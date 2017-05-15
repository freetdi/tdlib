#pragma once
#include "marker.hpp"

// mark neighbours of v up to v.
template<class M, typename V, class G>
void mark_smaller_neighbours(M& marker, V v, G const& g)
{ untested();
    auto pp=boost::adjacent_vertices(v, g);
    for(; pp.first!=pp.second; ++pp.first){
		 if(*pp.first>=v){ untested();
			 // break; // need sorting...
			 continue;
		 }else{ untested();
			 marker.mark(*pp.first);
		 }
    }
}

// FIXME: obsolete, use without mask
// mark neighbours of v up to v.
template<class M, typename V, class G, class MASK>
size_t mark_smaller_neighbours(M& marker, V v, G const& g, MASK const& m)
{
	size_t cnt=0;
	//std::cerr << "marking for " << v << "\n";
    auto pp=boost::adjacent_vertices(v, g);
    for(; pp.first!=pp.second; ++pp.first){
		 assert(*pp.first!=v);
		 if(!m[*pp.first]){
			 // masked...
		 }else if(*pp.first>=v){
			 // break; // need sorting...
		 }else{
			 // std::cerr << "marking " << *pp.first << "\n";
			 marker.mark(*pp.first);
			 ++cnt;
		 }
    }
	 return cnt;
}
