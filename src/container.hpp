// Felix Salfelder 2016
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

// container access

#ifndef TREEDEC_CONTAINER_HPP
#define TREEDEC_CONTAINER_HPP

// #include "config.h" // not yet
#include "container_traits.hpp"

#ifdef HAVE_TLX_CONTAINER_BTREE_SET_HPP
#include <tlx/container/btree_set.hpp>
#endif

namespace treedec{

template<class C>
typename C::value_type& front(C& c)
{
  return detail::container_iter<C>::front(c);
}
template<class C>
inline typename C::value_type back(C& c)
{
  return detail::container_iter<C>::back(c);
}

template<class C>
inline typename C::iterator begin(C& c)
{
  return detail::container_iter<C>::begin(c);
}

template<class C>
inline typename C::const_iterator begin(C const& c)
{
  return detail::container_iter<C>::begin(c);
}

template<class C>
inline typename C::iterator end(C& c)
{
  return detail::container_iter<C>::end(c);
}

template<class C>
inline typename C::const_iterator end(C const& c)
{
  return detail::container_iter<C>::end(c);
}

// size increases by 1
// i.e. implicitly, if it's not enforced by the container
template<class C, class E>
void push(C& c, E e)
{
    // n=size
    detail::container_modify<C>::push(c, e);
    // assert (size==n+1)
}

// push a range. otherwise like push
template<class C, class B, class E>
void push(C& c, B b, E e)
{
    detail::container_modify<C>::push(c, b, e);
}

// must not create duplicates
// i.e. implicitly, if it's not enforced by the container
template<class C, class E>
void insert(C& c, E e)
{
    // assert(!contains...)
    detail::container_modify<C>::insert(c, e);
}

template<class C, class B, class E>
void merge(C& c, B b, E e)
{
    detail::container_modify<C>::merge(c, b, e);
}

template<class C>
void sort(C& c)
{
    detail::container_modify<C>::sort(c);
}

template<class C, class E>
bool contains(C const& c, E e)
{
    return detail::container_inspect<C>::contains(c, e);
}

} // treedec

#endif
// vim:ts=8:sw=4:et
