// Felix Salfelder, 2017
//
// (c) 2014-2016 Goethe-Universität Frankfurt
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
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
//
//
//
// almost-tree-decompositions
// a skeleton is a list of vertices with vertices attached to it.

#ifndef TD_SKELETON_HPP
#define TD_SKELETON_HPP

#include <boost/graph/graph_traits.hpp>

namespace treedec{ //

namespace draft{

template<class G, class N, class O>
class SKELETON{
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
    SKELETON(G const& g, N const& n, O const& o)
      : _g(g), _n(n), _o(o)
    { untested();
    }
public:
    size_t size() const{ untested();
        return _o.size();
    }
//    bag_range operator[](vertices_size_type x){
//
//    }
private:
    G const& _g;
    N const& _n;
    O const& _o;
};

}


} // treedec

namespace boost{

template<class VD, class B>
std::pair<typename std::vector< std::pair<VD, B> >::const_iterator,
          typename std::vector< std::pair<VD, B> >::const_iterator >
vertices( std::vector< std::pair<VD, B> > const& skel )
{ untested();
    return std::make_pair(skel.begin(), skel.end());
}

template<class VD, class B>
VD get(vertex_owner_t, std::pair<VD, B> const& p,
         std::vector< std::pair<VD, B> > const& )
{
    return p.first;
}

template<class VD, class B>
B get(::treedec::bag_t, std::pair<VD, B> const& p,
         std::vector< std::pair<VD, B> > const& )
{
    return p.second;
}

template<class VD, class B>
B get(::treedec::bag_t, size_t u,
        std::vector< std::pair<VD, B> > const& X )
{
    return X[u].second;
}

} // boost

#endif // guard

// vim:ts=8:sw=4:et
