// Felix Salfelder 2016, 2021
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

// container access. FIXME: these are not traits.

#ifndef TREEDEC_CONTAINER_TRAITS_HPP
#define TREEDEC_CONTAINER_TRAITS_HPP

#include <vector>
#include <assert.h>
#include <set>
#include <algorithm> // sort
#include "trace.hpp"

namespace treedec{//

namespace detail{//
    template<class C, class X=void>
    struct container_iter{//
        // front, back: create copies. use (r)begin if you need references.
        static typename C::value_type front(C const& c)
        { untested();
          return *c.begin();
        }
        static typename C::value_type back(C const& c)
        { untested();
          return *c.rbegin();
        }
        // begin/end
        static typename C::iterator begin(C& c)
        { untested();
          return c.begin();
        }
        static typename C::const_iterator begin(C const& c)
        { untested();
          return c.begin();
        }
        static typename C::iterator end(C& c)
        { untested();
          return c.end();
        }
        static typename C::const_iterator end(C const& c)
        { untested();
          return c.end();
        }
        // rbegin/rend
        static typename C::iterator rbegin(C& c)
        { untested();
          return c.rbegin();
        }
        static typename C::const_iterator rbegin(C const& c)
        { untested();
          return c.rbegin();
        }
        static typename C::iterator rend(C& c)
        { untested();
          return c.rend();
        }
        static typename C::const_iterator rend(C const& c)
        { untested();
          return c.rend();
        }
    }; // access
    template<class C, class X=void>
    struct container_inspect{//
        template<class E>
        static bool contains(C const& c, E e)
        {
          return c.find(e)!=c.end();
        }
        // size?
        // is_ordered?
    };
    template<class C>
    struct container_inspect<C,
      typename std::enable_if<
        std::is_same<C, typename std::vector<
        typename C::value_type,
        typename C::allocator_type> >::value
        >::type >{//
        static bool contains(C const& c, typename C::value_type e) {
          return std::find(c.begin(), c.end(), e) != c.end();
        }
    };
    template<class C, class X=void>
    struct container_modify{
        // push, insert new item
        template<class E>
        static void push(C& c, E e) {
          // incomplete(); // not supported by all containers.
          // assert(!container_inspect<C>::contains(c, e));
          c.insert(e);
        }
        template<class E>
        static void insert(C& c, E e) { untested();
          c.insert(e);
        }
        // sort. no-op for some containers
        static void insert(C& c) { untested();
            incomplete();
        }
    };
    template<class C, class E, class X=void>
    struct container_assign{
        template<class O>
        static void assign(C& c, O const& o){ untested();
            container_modify<C>::push(c, o.begin(), o.end());
        }
        template<class O>
        static void assign(C& c, O&& o){
            container_modify<C>::push(c, o.begin(), o.end());
        }
    };
    template<class C>
    struct container_modify<C,
      typename std::enable_if<
        std::is_same<C, typename std::vector<
        typename C::value_type,
        typename C::allocator_type> >::value
        >::type >{//
        // push, insert new item
        static void push(C& c, typename C::value_type e) {
          assert(!container_inspect<C>::contains(c, e));
          c.push_back(e);
        }
        static void insert(C& c, typename C::value_type e) { untested();
          assert(!container_inspect<C>::contains(c, e));
          c.push_back(e);
        }
        template<class I>
        static void push(C& c, I b, I e) {
          for(; b!=e; ++b){
            push(c, *b);
          }
        }
        template<class I>
        static void merge(C& c, I b, I e) { untested();
          for(; b!=e; ++b){
            if(!container_inspect<C>::contains(c, *b)){
                insert(c, *b);
            }
          }
        }
        static void sort(C& c) {
            std::sort(c.begin(), c.end());
        }
    };
    template<class C, class E>
    struct container_assign<C, E,
      typename std::enable_if< std::is_same<
         typename std::remove_reference<C>::type,
         typename std::remove_reference<E>::type >::value >::type >{

        template<class O>
        static void assign(C& c, O const& o){ untested();
            for(auto x : o){ untested();
                push(c, x);
            }
        }
        template<class O>
        static void assign(C& c, O&& o){
            c = std::move(o);
        }
    };
    template<class C>
    struct container_modify<C,
      typename std::enable_if<
        std::is_same<C, typename std::set<
        typename C::value_type,
        typename C::key_compare,
        typename C::allocator_type> >::value
        >::type >{//
        // push, insert new item
        static void push(C& c, typename C::value_type e)
        {
          assert(!container_inspect<C>::contains(c, e));
          c.insert(e);
        }
        template<class I>
        static void push(C& c, I b, I e) {
          for(; b!=e; ++b){
            c.insert(*b);
          }
        }
        static void insert(C& c, typename C::value_type e)
        {
          c.insert(e);
        }
        template<class I>
        static void merge(C& c, I b, I e) {
          for(; b!=e; ++b){
            insert(c, *b);
          }
        }
        static void sort(C&) {
        }
    }; // container_mod<set>
} // detail

} // namespace treedec

#endif

// vim:ts=8:sw=4:et
