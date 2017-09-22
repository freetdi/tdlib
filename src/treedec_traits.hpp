// Felix Salfelder, 2016
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
//   traits for tdlib treedecs.
//
//

#ifndef TD_TREEDEC_TRAITS_HPP
#define TD_TREEDEC_TRAITS_HPP

#include "graph_traits.hpp"

namespace treedec{ //

namespace detail{ //

template<class B, class T, class V>
struct tmpbaghack{ //
    static typename treedec_traits<T>::bag_type& get_bag(T& t, V& v)
    {
        return t[v];
    }
    static typename treedec_traits<T>::bag_type const& get_bag(T const& t, V const& v)
    {
        return t[v];
    }
};

template<class T_t, class V>
struct tmpbaghack<bag_t, T_t, V>{ //
    static typename treedec_traits<T_t>::bag_type& get_bag(T_t& t, V& v)
    {
        return t[v].bag;
    }
    static typename treedec_traits<T_t>::bag_type const& get_bag(T_t const& t, V const& v)
    {
        return t[v].bag;
    }
};
} //detail

template<typename T_t>
inline typename treedec_traits<T_t>::bag_type& bag(
	const typename boost::graph_traits<T_t>::vertex_descriptor& v,
        T_t& T)
{
    typedef typename T_t::vertex_property_type b; //>::bag_type b;
    return detail::tmpbaghack<b,T_t,const typename boost::graph_traits<T_t>::vertex_descriptor&>::get_bag(T, v);
}

template<typename T_t>
inline typename treedec_traits<T_t>::bag_type const& bag(
        const typename boost::graph_traits<T_t>::vertex_descriptor& v,
        T_t const& T)
{
    typedef typename T_t::vertex_property_type b; //>::bag_type b;
    return detail::tmpbaghack<b,T_t,const typename boost::graph_traits<T_t>::vertex_descriptor&>::get_bag(T, v);
}

template<class V, class G>
size_t bag_size(V const & v, G const& g)
{
    return bag(g, v).size();
}

} //treedec

#endif
// vim:ts=8:sw=4:et
