// Felix Salfelder, 2014 - 2016
//
// (c) 2014-2016 Felix Salfelder
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

// graph related iteration

#ifndef TD_ITER_H
#define TD_ITER_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp> // adjacent_vertices
#include <tr1/utility> // pair

template<class A>
void debug_count(A i, A e, unsigned s)
{
    unsigned c=0;
    for(; i!=e; ++i){
        ++c;
    }
    assert(c==s);
}


template<class T>
class subsets_iter{ //
public: // types
   typedef typename std::vector<T> scratch_type;
   typedef typename scratch_type::const_iterator base;
   class subset_iter : public std::vector<T>::const_iterator{ //
   public:
      subset_iter(typename std::vector<T>::const_iterator const & i)
          : std::vector<T>::const_iterator(i)
      {
      }
      subset_iter(typename std::vector<T> const & t)
          : std::vector<T>::const_iterator(t.begin())
      { untested();
      }
      subset_iter(T const & i)
      { untested();
          _t.resize(0);
          _t.push_back(i);
      }
   public: // ops
      typename T::value_type const& operator*()
      {
         return *std::vector<T>::const_iterator::operator*();
      }
      bool operator!=(const T& other)
      { untested();
         return !std::vector<T>::const_iterator::operator==(other);
      }
      bool operator==(const T& other)
      { untested();
         return std::vector<T>::const_iterator::operator==(other);
      }
   }; // subset_iter
public: // construct
   subsets_iter(const subsets_iter& p)
       : _tt(NULL),
         _t(p._t),
         _i(p._i), _e(p._e), _l(p._l), _u(p._u)
   {
   }
#if __cplusplus >= 201103L
   subsets_iter(const subsets_iter&& p)
       : _tt(p._tt),
         _t(p._t),
         _i(p._i), _e(p._e), _l(p._l), _u(p._u)
   {
       p._tt = NULL;
   }
#endif
   subsets_iter(T i, T e, size_t min=0, size_t max=-1, scratch_type* t=NULL)
       : _tt(NULL), _t(t?(*t):(*(new scratch_type))),
         _i(i), _e(e), _l(min), _u(max)
   {
#if __cplusplus >= 201103L
      if(t){
          t->resize(0);
          // use external scratch
          _tt = NULL;
      }else{
          // own scratch. record it.
          _tt = &_t;
      }
#endif
      assert(_l<=_u);
      fill();
   }
#if 0
   subsets_iter(T e)
       : _e(e), _l(0), _u(-1)
   { untested();
      _t.resize(1);
      _t[0] = e;
   }
#endif
   ~subsets_iter(){
      if(_tt){
          // own scratch. delete it.
          delete(&_t);
      }else{
      }
   }
public: // assign
   subsets_iter& operator=(const subsets_iter& other)
   { untested();
      _i = other._i;
      _e = other._e;
      _l = other._l;
      _u = other._u;
      _t = other._t;
//      _tt = other._tt;

//      other._tt=NULL;
      return *this;
   }
   subsets_iter& operator=(const T& other)
   { untested();
      base::operator=(other);
      _i = other._i;
      _e = other._e;
      _l = other._l;
      _u = other._u;
      _t = other._t;
      return *this;
   }
public: // ops
   bool operator==(const T& other)
   { itested();
      if(!_t.size()){ untested();
         // BUG?
         return false;
      }else if(_t.size()){
         return *_t.begin() == other;
      }else{ untested();
         return false;
      }
   }
   bool operator!=(const T& other)
   {
      return !operator==(other);
   }
   bool operator==(const subsets_iter& other)
   { untested();
      if(_t.size()){ untested();
          return _t.size() == other._t.size()
              && *_t[0]==*other._t[0];
      }else{ untested();
          return _i==other._i && _e==other._e;
      }
   }
   bool operator!=(const subsets_iter& other)
   { untested();
      return !operator==(other._i);
   }
   subsets_iter operator++()
   {
      if(_t.size()==0){ untested();
         _t.push_back(_i);
         assert(_i!=_e); // probably not.
                         // HACK: always push back something...
         if(_u==0){ untested();
            _t.back()=_e;
         }else{untested();
         }
      }else if(_t.size()<=_u){
         BOOST_AUTO(back,_t.back());
         ++back;
         if(back!=_e){
            if(_t.size()==_u){
               ++_t.back();
            }else{
               _t.push_back(back);
            }
         }else if(_t.back()==_e){ untested();
            unreachable();
         }else if(_t.size()==0){ untested();
            unreachable();
         }else if(_t.size()==1){
            ++_t.back();
         }else if(_t.size()==_l){
             BOOST_AUTO(back2, _t.back());
             carry(back2);
             if(_t.size()<_l){
                 _t[0]=_e;
             }else{
             }
         }else{
            _t.pop_back();
            BOOST_AUTO(back2, _t.back());
            ++back2;
            if(back2!=_e){
               ++_t.back();

            }else{ unreachable();
               // carry();
            }
         }
      }else if(_t.back() != _e){ untested();
         // assert(_t.size()==_u);
         incomplete();
      }else{ untested();
      }

//      fill();
      assert(_t.size()<=_u || _u==0);

      return *this;
   }
   void carry(T end)
   {
       BOOST_AUTO(b, _t.back());
       ++b;
       if(_t.size() == 1){
           ++_t.back();
       }else if(_t.back() == end){
           _t.pop_back();
           BOOST_AUTO(newend, _t.back());
           ++newend;
           if(newend==end){
               newend=_t.back();
           }else{
           }
           carry(newend);

           b = _t.back();
           ++b;
           if(_t.back()!=end){
               _t.push_back(b);
           }
       }else{
           assert(b==end);
           ++_t.back();
       }

   }
   void fill()
   {
      while(_t.size() < _l){
         if(!_t.size()){
            _t.push_back(_i);
         }else{
            BOOST_AUTO(back, _t.back());
            if(back==_e){ itested();
               break;
            }else{
               ++back;
               if(back==_e){
                       _t[0] = _e;

                   break;
               }
               _t.push_back(back);
            }
         }
      }

      assert(_t.size()>=_l || _t[0] == _e);
   }
   std::pair<subset_iter, subset_iter> operator*()
   {
      return std::make_pair(subset_iter(_t.begin()),
                            subset_iter(_t.end()));
   }
private: // detail
   void carry()
   { untested();
   }
private: // state
   mutable scratch_type* _tt;
   scratch_type& _t;
   T _i;
   /*const*/ T _e;
   /*const*/ size_t _l;
   /*const*/ size_t _u;
}; // subsets_iter

template<class A>
std::pair<subsets_iter<A>, subsets_iter<A> >
    make_subsets_iter(A a, A b, unsigned l, unsigned u,
          typename subsets_iter<A>::scratch_type* s=NULL)
{
   return std::make_pair(
       subsets_iter<A>(a,b,l,u,s),
       subsets_iter<A>(b,b));
}

#if 0
template<class G>
class nvsinconwsut_iter{
public: // types
   typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
public: // ops
private: // status
   std::vector<vertex_descriptor> _base;
};

template<class G>
std::pair<nvsinconwsut_iter<G>, nvsinconwsut_iter<G> >
    nonempty_vertex_sets_in_complement_of_neighborhood_with_size_up_to(*I, 2*k, g)
{
   return std::make_pair(nvsinconwsut_iter(*I, 2*k, g), nvsinconwsut_iter());
}
#endif

namespace detail{
// iterate neighbourhood of C, including C
// only works on ordered containers...
template<class A, class G>
class neighbourhood01_iter { //
public: // types
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
//     neighbourhood01_iter(const neighbourhood01_iter& p)
//     { untested();
//      incomplete. reuse _a?
//     }
    neighbourhood01_iter() : _v(NULL)
    { untested();
    }
    neighbourhood01_iter(A b, A e, unsigned size, const G& g, bool incb)
       : _b(b), _i(b), _e(e),
         _a(size), _g(g), _include_base(incb)
    {
        if(b==e){
            assert(size==0);
            return;
        }else{
        }

        bool found=false;

        if(_include_base){ untested();
            _v = *_i;
        }else{ untested();
            A ii(_b);
            for(; ii!=_e; ++ii){
                // just take one.... (inefficient)
                if(boost::adjacent_vertices(*ii, g).first
                 !=boost::adjacent_vertices(*ii, g).second){
                    _v = *boost::adjacent_vertices(*ii, g).first;
                    found = true;
                    break;
                }
            }
        }

        A ii(_b);

        unsigned n=0;
        for(; ii!=_e; ++ii){
            assert(!size || n<size);
            if(!size){ untested();
                _a.push_back(boost::adjacent_vertices(*ii, g).first);
            }else{ untested();
                _a[n] = boost::adjacent_vertices(*ii, g).first;
            }

            if(_a[n] == boost::adjacent_vertices(*ii, g).second){ untested();
            }else if(*(_a[n])<_v){ untested();
                _v = *_a[n];
                found = true;
            }else{ untested();
            }

            ++n;
        }
        if(_include_base){ untested();
        }else if(!found){ untested();
            _b = _e;
        }
    }

    ~neighbourhood01_iter()
    { itested();
    }
public: // ops
    bool operator!=(const neighbourhood01_iter& p) const
    {
        return _b != p._b;
    }
    bool operator==(const neighbourhood01_iter& p) const
    {
        return(_b == p._b);
    }
    neighbourhood01_iter& operator++()
    {
        if(_b==_e){ untested();
            return *this;
        }else{
        }

        vertex_descriptor previous = _v;
        bool found = false;
        if(_include_base){ untested();
            found = update(_i, _e, previous, _v);
        }else{ untested();
        }
        A ii(_b);

        unsigned n=0;
        for(; ii!=_e; ++ii){ untested();
            BOOST_AUTO(aend, boost::adjacent_vertices(*ii, _g).second);
            found |= update(_a[n], aend, previous, _v);
            ++n;
        }

        if(!found){ untested();
            _b = _e;
        }else{ untested();
            assert(_v>previous);
        }

        return *this;
    }
    const vertex_descriptor& operator*()
    {
        //std::cerr << "iter deref ++ " << _v << "\n";
        return _v;
    }
private: // impl
    template<class iter>
    bool update(iter& i, iter e, vertex_descriptor previous, vertex_descriptor& v)
    {
        if(i==e){
            return false;
        }else if(*i==previous){
            ++i;
            if(i==e){
                return false;
            }
        }else{
            //std::cerr << *i << " " << previous << "\n";
            assert(*i>previous);
        }

        if(v==previous){
            v = *i;
        }else if(*i<v){
            assert(*i>previous);
            v = *i;
        }else{
        }
        return true;
    }
private: // data
    A _b;
    A _i;
    A _e;
    std::vector<adjacency_iterator> _a;
    vertex_descriptor _v;
    G const& _g;
    bool _include_base;
};
} // detail



template<class A, class G>
std::pair<detail::neighbourhood01_iter<A, G>,
          detail::neighbourhood01_iter<A, G> >
inline make_neighbourhood01_iter(A b, A e, G const& g, unsigned size=0,
        bool include_base=true)
{
    typedef detail::neighbourhood01_iter<A, G> nIter;
    return std::make_pair(
            nIter(b,e,size,g,include_base),
            nIter(e,e,0,g,include_base) );
}

template<class A, class G>
std::pair<detail::neighbourhood01_iter<A, G>,
          detail::neighbourhood01_iter<A, G> >
inline make_neighbourhood_iter(A b, A e, G const& g, unsigned size=0)
{
    return make_neighbourhood01_iter(b,e,g,size,false);
}

#endif

// vim:ts=8:sw=4:et
