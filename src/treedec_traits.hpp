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
//   traits for treedec treedecs.
//
//

#ifndef TREEDEC_TREEDEC_TRAITS_HPP
#define TREEDEC_TREEDEC_TRAITS_HPP

#include "graph_traits.hpp"

namespace treedec{ //

namespace detail{ //

template<class B, class T, class V>
struct tmpbaghack{ //
    static typename treedec_traits<T>::bag_type& get_bag(T& t, V& v)
    { incomplete();
        return t[v];
    }
    static typename treedec_traits<T>::bag_type const& get_bag(T const& t, V const& v) {
        incomplete();
        return t[v];
    }
};

template<class T_t, class V>
struct tmpbaghack<bag_t, T_t, V>{
    static typename treedec_traits<T_t>::bag_type& get_bag(T_t& t, V& v) {
        incomplete();
        return t[v].bag;
    }
    static typename treedec_traits<T_t>::bag_type const& get_bag(T_t const& t, V const& v) {
        incomplete();
        return t[v].bag;
    }
};
} //detail

#if 1 // BUG. one call left in graph.hpp
template<typename T_t>
inline typename treedec_traits<T_t>::bag_type const& bag(
        const typename boost::graph_traits<T_t>::vertex_descriptor& v,
        T_t const& T)
{
    incomplete();
    typedef typename T_t::vertex_property_type b; //>::bag_type b;
    return detail::tmpbaghack<b,T_t,const typename boost::graph_traits<T_t>::vertex_descriptor&>::get_bag(T, v);
}
#endif

template<class V, class G>
size_t bag_size(V const & v, G const& g)
{
    return bag(g, v).size();
}

} //treedec

//namespace boost{
//template<class T, class V>
//void get(treedec::bag_t, T const& t, V v)
//{ untested();
//}
//template<class T, class V>
//void get(treedec::bag_t, T& t, V v)
//{ untested();
//}
//}
//
// uses treedec::push, include container...
// BUG: only works for "set" and "vector" of unsigned
// BUG: BAG must be bag, still used in other traits :|
#define TREEDEC_TREEDEC_BAG_TRAITS(T, BAG)\
namespace boost{\
\
    inline \
    typename property_map< T, vertex_all_t>::type\
    get(vertex_all_t, T& g) \
    {\
        typedef typename property_map< T, vertex_all_t>::type\
          pmap_type;\
        return pmap_type(g);\
    }\
\
    inline \
    bagstuff::const_treebagpmap<T>\
    get(vertex_all_t, T const& g) \
    {\
        typedef typename property_map< T, vertex_all_t>::const_type\
            pmap_type;\
        return pmap_type(g);\
    }\
\
    template<class U> \
    inline void\
    put(put_get_helper<bagstuff::gtob<T>::type,\
            bagstuff::treebagpmap<T> >& pa, U k, treedec::bag_t const& v)\
    { untested1("....");\
        auto& PA=static_cast<bagstuff::treebagpmap<T>  const&>(pa);\
        auto& b=const_cast<bagstuff::treebagpmap<T> &>(PA)[k];\
        b.clear();\
        for(auto const& i : v.BAG){ untested();\
            treedec::push(b, i);\
        }\
    }\
\
    template<class U> \
    inline void\
    put(const put_get_helper<bagstuff::gtob<T>::type,\
          bagstuff::treebagpmap<T> >& pa, U k,\
          const property<treedec::bag_t, std::set<unsigned> >& v)\
    {\
        auto& PA=static_cast<bagstuff::treebagpmap<T>  const&>(pa);\
        auto& b=const_cast<bagstuff::treebagpmap<T> &>(PA)[k];\
        b.clear();\
        for(auto const& i : v.m_value){ untested();\
            treedec::push(b, i);\
        }\
    }\
\
    template<class U> \
    inline void\
    put(const put_get_helper<bagstuff::gtob<T>::type,\
        bagstuff::treebagpmap<T> >& pa, U k,\
        const property<treedec::bag_t, std::vector<unsigned> >& v)\
    {\
        auto& PA=static_cast<bagstuff::treebagpmap<T>  const&>(pa);\
        auto& b=const_cast<bagstuff::treebagpmap<T> &>(PA)[k];\
        b.clear();\
        for(auto const& i : v.m_value){ untested();\
            treedec::push(b, i);\
        }\
    }\
\
    template<class U> \
    inline bagstuff::gtob<T>::type const& \
    get(treedec::bag_t, T const&t, U k)\
    {\
        return t[k].bag;\
    }\
\
    template<class U> \
    inline bagstuff::gtob<T>::type& \
    get(treedec::bag_t, T &t, U k)\
    {\
        return t[k].bag;\
    }\
\
    inline bagstuff::const_treebagpmap<T> \
    get(treedec::bag_t, T const& t)\
    {\
        return bagstuff::const_treebagpmap<T>(t);\
    }\
\
    inline bagstuff::treebagpmap<T> \
    get(treedec::bag_t, T & t)\
    { untested(); \
        return bagstuff::treebagpmap<T>(t);\
    }\
\
    template <>\
    struct property_map<T, treedec::bag_t>{ \
    };\
} /* boost */ \
namespace treedec{ \
    template<> \
    struct treedec_traits<T>{ \
        typedef typename T::vertex_property_type vertex_property_type; \
        typedef typename boost::bagstuff::gtob<T>::type bag_type; \
        typedef typename boost::bagstuff::gtob<T>::type::value_type vd_type; \
    }; \
} /* treedec */ \
void TREEDEC_DUMMY_FUNCTION_DECLARATION(void)

#endif
// vim:ts=8:sw=4:et
