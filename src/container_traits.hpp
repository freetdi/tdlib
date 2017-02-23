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

#ifndef TD_CONTAINER_TRAITS_HPP
#define TD_CONTAINER_TRAITS_HPP

#include <vector>
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
        static bool contains(C& c, typename C::value_type e)
        {
          (void) c; (void) e;
          incomplete();
          return false;
        }
        // size?
        // is_ordered?
    };
    template<class C, class X=void>
    struct container_modify{//
        // push, insert new item
        static void push(C& c, typename C::value_type e)
        {
          assert(!container_inspect<C>::contains(c, e));
          c.insert(e);
        }
        static void insert(C& c, typename C::value_type e)
        {
          c.insert(e);
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
        static void push(C& c, typename C::value_type e)
        { untested();
          assert(!container_inspect<C>::contains(c, e));
          c.push_back(e);
        }
        static void insert(C& c, typename C::value_type e)
        { untested();
          assert(!container_inspect<C>::contains(c, e));
          c.insert(e);
        }
    };
} // detail

} // namespace treedec

#endif

// vim:ts=8:sw=2:et
