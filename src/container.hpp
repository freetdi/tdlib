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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

// container access

#ifndef TD_CONTAINER_HPP
#define TD_CONTAINER_HPP

// #include "config.h" // not yet
#include "container_traits.hpp"

#ifdef HAVE_STX_BTREE_SET_H
#include <stx/btree_set>
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

// must not create duplicates
// i.e. implicitly, if it's not enforced by the container
template<class C, class E>
void insert(C& c, E e)
{
    // assert(!contains...)
    detail::container_modify<C>::insert(c, e);
}

template<class C, class E>
bool contains(C const& c, E e)
{ itested();
    return detail::container_inspect<C>::contains(c, e);
}

} // treedec

#endif
// vim:ts=8:sw=4:et
