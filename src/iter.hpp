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
   typedef typename std::vector<T>::const_iterator base;
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
   };
public: // construct
   subsets_iter(const subsets_iter& p)
       : _t(p._t),_i(p._i), _e(p._e), _l(p._l), _u(p._u)
   { untested();
   }
   subsets_iter(T i, T e, size_t min=0, size_t max=-1)
       : _i(i), _e(e), _l(min), _u(max)
   { untested();
      assert(_l<=_u);

      fill();
   }
   subsets_iter(T e)
       : _e(e), _l(0), _u(-1)
   { untested();
      _t.resize(1);
      _t[0] = e;
   }
public: // ops
   void operator=(const T& other)
   { untested();
      return base::operator=(other);
		_i = other._i;
		_e = other._e;
		_l = other._l;
		_u = other._u;
		_t = other._t;
   }
   bool operator==(const T& other)
   {
      if(!_t.size()){ untested();
         // BUG?
         return false;
      }else if(_t.size()){ untested();
         return *_t.begin() == other;
      }else{ untested();
         untested();
         return false;
      }
   }
   bool operator!=(const T& other)
   { untested();
      return !operator==(other);
   }
   bool operator==(const subsets_iter& other)
   { untested();
      return _i==_e && other._i==other._e;
   }
   bool operator!=(const subsets_iter& other)
   { incomplete();
      return _i!=other._i || _e!=other._e;
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
      }else if(_t.size()<=_u){ untested();
         BOOST_AUTO(back,_t.back());
         ++back;
         if(back!=_e){ untested();
            if(_t.size()==_u){ untested();
               ++_t.back();
            }else{ untested();
               _t.push_back(back);
            }
         }else if(_t.back()==_e){
            unreachable();
         }else if(_t.size()==0){
            unreachable();
         }else if(_t.size()==1){
            ++_t.back();
         }else{ untested();
            _t.pop_back();
            BOOST_AUTO(back2, _t.back());
            ++back2;
            if(back2!=_e){ untested();
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
      while(_t.size() < _l){ untested();
         if(!_t.size()){ untested();
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
   std::vector<T> _t;
   T _i;
   /*const*/ T _e;
   /*const*/ size_t _l;
   /*const*/ size_t _u;
}; // subsets_iter

template<class A>
std::pair<subsets_iter<A>, subsets_iter<A> >
    make_subsets_iter(A const& a, A const& b, unsigned l, unsigned u)
{
   return std::make_pair(
		 subsets_iter<A>(a,b,l,u),
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
