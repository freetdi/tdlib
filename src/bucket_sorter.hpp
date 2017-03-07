//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//
//
// Revision History:
//   13 June 2001: Changed some names for clarity. (Jeremy Siek)
//   01 April 2001: Modified to use new <boost/limits.hpp> header. (JMaddock)
//
//   2016: iterable bucket stacks (Felix Salfelder)
//
#ifndef BOOST_GRAPH_DETAIL_BUCKET_SORTER_HPP
#define BOOST_GRAPH_DETAIL_BUCKET_SORTER_HPP

#include <vector>
#include <cassert>
#include <boost/limits.hpp>
#include "trace.hpp"

namespace boost {

  template <class BucketType, class ValueType, class Bucket, 
            class ValueIndexMap>
  class bucket_sorter {
  public:
    typedef BucketType bucket_type;
    typedef ValueType value_type;
    typedef Bucket Bucket_type;
    typedef ValueIndexMap value_index_map;
    typedef typename std::vector<value_type>::size_type size_type;
    
    bucket_sorter(size_type _length, bucket_type _max_bucket, 
                  const Bucket& _bucket = Bucket(), 
                  const ValueIndexMap& _id = ValueIndexMap()) 
      : next(_length+_max_bucket, invalid_value()), 
        prev(_length, invalid_value()),
        head(next.size()?(&next[_length]):NULL),
        id_to_value(_length),
        bucket(_bucket), id(_id) {
          trace2("created bs", _length, _max_bucket);
          assert(head==&next[_length]);
          assert(!next.size() || head[0]==invalid_value());
        }
    bucket_sorter(){untested();
    }

    void remove(const value_type& x) { untested();
      const size_type i = get(id, x);
      assert(x<size());
#if 0
      trace3("rm", x, i, bucket[x]);
      trace1("rm", head[bucket[x] ]);
      trace2("rm", prev[i], next[i]);
#endif
      const size_type& next_node = next[i];
      const size_type& prev_node = prev[i];
      assert(prev_node!=i);
    
      //check if i is the end of the bucket list 
      if ( next_node != invalid_value() ){ untested();
        assert(next_node != prev_node);
        prev[next_node] = prev_node; 
      }
      //check if i is the begin of the bucket list
      if( prev_node == invalid_value() ){ untested();
        unreachable();
        // double remove?
      }else if(prev_node>=size()){ untested();
        next[prev_node] = next_node;
      }else{ untested();
        next[prev_node] = next_node;
      }
      if(next_node == prev_node){ untested();
        // double remove?
      }

      assert(head[bucket[x]]!=x);
#if 0
      assert(head[bucket[x]]!=x);
      untested();
      trace3("rdm", x, i, bucket[x]);
      trace1("rmd", head[bucket[x] ]);

      auto const& T=*this;
      const const_stack& b = T[bucket[x]];
      trace1("rmd", b.top());
      trace2("rmd", prev[i], next[i]);
#endif
      next[i]=invalid_value(); // BUG, why needed?
#ifndef NDEBUG
      prev[i]=invalid_value();
#endif
    }

    void push_back(const value_type& x) { untested();
      assert(x<next.size());
      id_to_value[get(id, x)] = x;
      (*this)[bucket[x]].push_back(x);
    }

    void push(const value_type& x) { untested();
      assert(x<next.size());
      id_to_value[get(id, x)] = x;
      (*this)[bucket[x]].push(x);
    }
#ifndef NDEBUG
    bool is_known(const value_type& v) const { untested();

      //trace4("is_known", v, k, id_to_value[v], bucket[v]);
      //trace4("is_known", v, k, id_to_value[v], head[bucket[v]]);
      // if it exists, it is at the top, or
      // it has a predecessor.
      const const_stack& b = (*this)[bucket[v]];
      if( prev[id_to_value[v]] < next.size()){ untested();
        return true;
      }else if( b.empty()){ itested();
      //  trace3("empty", v, bucket[v], head[bucket[v]]);
        return false;
      } else if( head[bucket[v]] == v){ itested();
        return true;
      }else{ itested();
        //trace3("buckettop is", v, k, b.top());
        return false;
      }
    }
#endif

    void update(const value_type& x) { untested();
      remove(x); // can kill head.
      (*this)[bucket[x]].push(x);
    }

#if 0
    // does not help
    void update(const value_type& x, const size_type& old_val ) { untested();
      assert(bucket[x] == old_val);
      remove(x, old_val);
      (*this)[bucket[x]].push(x);
    }

    // this does not work, since bucket[x] is not necessarily
    // writable
    void update(const value_type& x, const size_type& new_val ) { untested();
      remove(x);
      bucket[x] = new_val;
      (*this)[bucket[x]].push(x);
    }
#endif
    //  private: 
    //    with KCC, the nested stack class is having access problems
    //    despite the friend decl.
    static size_type invalid_value() { untested();
      return (std::numeric_limits<size_type>::max)();
    }
    
    typedef typename std::vector<size_type>::iterator Iter;
    typedef typename std::vector<size_type>::const_iterator ConstIter;
    typedef typename std::vector<value_type>::iterator IndexValueMap;
    typedef typename std::vector<value_type>::const_iterator ConstIndexValueMap;
    
  public:

    template<class Iter_, class IndexValueMap_>
    class stack_ {
    public:
      typedef bucket_sorter base;
      typedef bucket_sorter::value_type value_type;
    public:
      class const_iterator{
        // bug? feature?
        const_iterator(){unreachable();}
      public:
        const_iterator(size_type t, stack_ const& s_)
           : s(s_), b(t) {}
        const_iterator(const const_iterator& p)
           : s(p.s), b(p.b) {}
        const_iterator(const const_iterator&& p)
           : s(p.s), b(p.b) {untested();}
        ~const_iterator(){}
      public:
        const_iterator& operator=(const const_iterator& o){ untested();
          assert(&s==&o.s); // how to compile time check?!
                            // or just fallback to pointer?
          b = o.b;
          return *this;
        }
        const_iterator& operator=(const const_iterator&& o){ untested();
          assert(&s==&o.s); // how to compile time check?!
                            // or just fallback to pointer?
          b = o.b;
          return *this;
        }
        value_type operator*() const{ untested();
          assert(b<s.size());
          return s.value[b];
        }
        const_iterator& operator++(){ untested();
          assert(b!=invalid_value());
          assert(b!=s.next[b]);
          b = s.next[b];
          return *this;
        }
        bool operator!=(const_iterator const& o){ untested();
          return o.b!=b;
        }
        bool operator==(const_iterator const& o)
        { return o.b==b; }
      private:
        stack_ const& s;
        size_type b;
      };
    public:
      stack_(const stack_& p)
        : bucket_id(p.bucket_id),
          head(p.head),
          next(p.next),
          prev(p.prev),
          value(p.value),
          id(p.id)
      {untested();
      }
      // stack_(const stack_&&){untested();}
    public:
      stack_(bucket_type _bucket_id, typename Iter_::value_type* h, Iter_ n, Iter_ p, IndexValueMap_ v,
            const ValueIndexMap& _id)
      : bucket_id(_bucket_id), head(h), next(n), prev(p), value(v), id(_id) {}

      // Avoid using default arg for ValueIndexMap so that the default
      // constructor of the ValueIndexMap is not required if not used.
      stack_(bucket_type _bucket_id, typename Iter_::value_type* h, Iter_ n, Iter_ p, IndexValueMap_ v)
        : bucket_id(_bucket_id), head(h), next(n), prev(p), value(v) { untested(); }


      void push_back(const value_type& x) { untested();
        incomplete();
      }
      void push(const value_type& x) { untested();
        const size_type new_head = get(id, x);
        assert(new_head < size());
        const size_type current = head[bucket_id];
        if(new_head == current){ untested();
//          assert(false);
        }
        if ( current != invalid_value() ){ untested();
          assert(current!=new_head);
          prev[current] = new_head;
        }

        prev[new_head] = bucket_id + (head - next);
        next[new_head] = current;
        head[bucket_id] = new_head;
      }
      void pop() { untested();
        assert(bucket_id<size());
        size_type current = head[bucket_id];
        size_type next_node = next[current];
        head[bucket_id] = next_node;
        // prev[bucket_id] = bucket_id + (head - prev);
        if ( next_node != invalid_value() ){ untested();
          assert(next_node != prev[current]);
          prev[next_node] = bucket_id + (head - next);
        }
        assert(next_node != prev[current]);
      }
      value_type const& top() const { return value[ head[bucket_id] ]; }
      value_type& top() { return value[ head[bucket_id] ]; }
      bool empty() const { return head[bucket_id] == invalid_value(); }
    public: // iterator access
      const_iterator begin() const{ untested();
        return const_iterator(head[bucket_id], *this);
      }
      // BUG: template override in degree.hpp does not match (why?)
      const_iterator rbegin() const{ untested();
        return const_iterator(head[bucket_id], *this);
      }
      const_iterator end() const{ return const_iterator(invalid_value(), *this); }
    public: // debug
      size_t size()const{return head-next;}
    private:
      bucket_type bucket_id;
      Iter_ head;
      Iter_ next;
      Iter_ prev;
      IndexValueMap_ value;
      ValueIndexMap id;
    }; // stack

    typedef stack_<Iter, IndexValueMap> stack;
    typedef stack_<ConstIter, ConstIndexValueMap> const_stack;
    
    const_stack operator[](const bucket_type& i) const{ untested();
      assert(i < next.size());

      return const_stack(i, head, next.begin(), prev.begin(),
                   id_to_value.begin(), id);
    }
    stack operator[](const bucket_type& i) { untested();
      assert(i < next.size());
      return stack(i, head, next.begin(), prev.begin(),
                   id_to_value.begin(), id);
    }
    unsigned size() const{ return prev.size(); }
  protected:
    std::vector<size_type>   next;
    std::vector<size_type>   prev;
    typename std::vector<size_type>::value_type* head;
    std::vector<value_type>  id_to_value;
    Bucket bucket;
    ValueIndexMap id;
  };
  
}

#endif

// vim:ts=8:sw=2:et
