// Felix Salfelder, 2017
//
// (c) 2017 Felix Salfelder
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
#ifndef TREEDEC_MARKER_UTIL_HPP
#define TREEDEC_MARKER_UTIL_HPP

#include "marker.hpp"
//#include "graph_util.hpp" // no. used by graph_util

namespace treedec{

template<class I, class M>
void mark_range(I i, I e, M& marker)
{
    for(; i!=e; ++i){
		 marker.mark(*i);
    }
}

// mark neighbours of v up to v.
// count P, hack.
template<class M, typename V, class G, class P>
size_t mark_neighbours_c(M& marker, V v, G const& g, P const& p /*bug*/)
{
    size_t count=0;
    auto pp=boost::adjacent_vertices(v, g);
    for(; pp.first!=pp.second; ++pp.first){
        marker.mark(*pp.first);
        if(p(*pp.first)){
            ++count;
        }else{
        }
    }
    return count;
}

template<class M, typename V, class G>
void mark_neighbours(M& marker, V v, G const& g)
{
    auto pp=boost::adjacent_vertices(v, g);
    for(; pp.first!=pp.second; ++pp.first){
        marker.mark(*pp.first);
    }
}

} // treedec

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

#endif
