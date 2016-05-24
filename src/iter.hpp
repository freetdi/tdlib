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
#include <tr1/utility> // pair

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
   { untested();
       p._tt = NULL;
   }
#endif
   subsets_iter(T i, T e, size_t min=0, size_t max=-1, scratch_type* t=NULL)
       : _tt(NULL), _t(t?(*t):(*(new scratch_type))),
         _i(i), _e(e), _l(min), _u(max)
   {
#if __cplusplus >= 201103L
      if(t){untested();
          t->resize(0);
          // use external scratch
          _tt = NULL;
      }else{untested();
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
      if(_tt){ untested();
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
   {
      if(!_t.size()){ untested();
         // BUG?
         return false;
      }else if(_t.size()){
         return *_t.begin() == other;
      }else{ untested();
         untested();
         return false;
      }
   }
   bool operator!=(const T& other)
   {
      return !operator==(other);
   }
   bool operator==(const subsets_iter& other)
   {
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
         }else if(_t.back()==_e){
            unreachable();
         }else if(_t.size()==0){
            unreachable();
         }else if(_t.size()==1){
            ++_t.back();
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

      fill();
      assert(_t.size()<=_u || _u==0);

      return *this;
   }
   void fill()
   {
      while(_t.size() < _l){
         if(!_t.size()){
            _t.push_back(_i);
         }else{ untested();
            BOOST_AUTO(back, _t.back());
            if(back==_e){ untested();
               break;
            }else{ untested();
               ++back;
               if(back==_e){ untested();
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
   {
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

#endif

// vim:ts=8:sw=4:et
