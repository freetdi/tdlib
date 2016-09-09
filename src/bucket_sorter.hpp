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
      : head(_max_bucket, invalid_value()),
        next(_length, invalid_value()), 
        prev(_length, invalid_value()),
        id_to_value(_length),
        bucket(_bucket), id(_id) { }
    
    void remove(const value_type& x) {
      const size_type i = get(id, x);
      const size_type& next_node = next[i];
      const size_type& prev_node = prev[i];
    
      //check if i is the end of the bucket list 
      if ( next_node != invalid_value() )
        prev[next_node] = prev_node; 
      //check if i is the begin of the bucket list
      if ( prev_node != invalid_value() )
        next[prev_node] = next_node;
      else //need update head of current bucket list
        head[ bucket[x] ] = next_node;
    }

    void push(const value_type& x) {
      id_to_value[get(id, x)] = x;
      (*this)[bucket[x]].push(x);
    }
    
    void update(const value_type& x) {
      remove(x);
      (*this)[bucket[x]].push(x);
    }
    //  private: 
    //    with KCC, the nested stack class is having access problems
    //    despite the friend decl.
    static size_type invalid_value() {
      return (std::numeric_limits<size_type>::max)();
    }
    
    typedef typename std::vector<size_type>::iterator Iter;
    typedef typename std::vector<value_type>::iterator IndexValueMap;
    
  public:

    friend class stack;

    class stack {
    public:
      typedef bucket_sorter base;
      typedef bucket_sorter::value_type value_type;
    public:
      class const_iterator{
      public:
        const_iterator(size_type t, stack const& s_)
           : s(s_), b(t) {}
        const_iterator(const const_iterator& p)
           : s(p.s), b(p.b) {}
      public:
        value_type operator*() const{
          trace2("*", b, s.value[b]);
          assert(b!=invalid_value());
          return s.value[b];
        }
        const_iterator& operator++(){
          assert(b!=invalid_value());
          b = s.next[b];
          return *this;
        }
        bool operator!=(const_iterator const& o)
		  {
			  trace2("!=", b, o.b);
			  untested(); return o.b!=b; }
        bool operator==(const_iterator const& o)
		  { untested(); return o.b==b; }
      private:
        stack const& s;
        size_type b;
      };
    public:
      stack(bucket_type _bucket_id, Iter h, Iter n, Iter p, IndexValueMap v,
            const ValueIndexMap& _id)
      : bucket_id(_bucket_id), head(h), next(n), prev(p), value(v), id(_id) {}

      // Avoid using default arg for ValueIndexMap so that the default
      // constructor of the ValueIndexMap is not required if not used.
      stack(bucket_type _bucket_id, Iter h, Iter n, Iter p, IndexValueMap v)
        : bucket_id(_bucket_id), head(h), next(n), prev(p), value(v) {}

      void push(const value_type& x) {
        const size_type new_head = get(id, x);
        const size_type current = head[bucket_id];
        if ( current != invalid_value() )
          prev[current] = new_head;
        prev[new_head] = invalid_value();
        next[new_head] = current;
        head[bucket_id] = new_head;
      }
      void pop() {
        size_type current = head[bucket_id];
        size_type next_node = next[current];
        head[bucket_id] = next_node;
        if ( next_node != invalid_value() )
          prev[next_node] = invalid_value();
      }
      value_type& top() { return value[ head[bucket_id] ]; }
      const value_type& top() const { return value[ head[bucket_id] ]; }
      bool empty() const { return head[bucket_id] == invalid_value(); }
    public: // iterator access
      const_iterator begin() const{
			trace2("", bucket_id, head[bucket_id]);
			return const_iterator(head[bucket_id], *this);
		}
		// BUG: template override in degree.hpp does not match (why?)
      const_iterator rbegin() const{
			return const_iterator(head[bucket_id], *this);
		}
      const_iterator end() const{ return const_iterator(invalid_value(), *this); }
	 private: // BUG
      bucket_type bucket_id;
      Iter head;
      Iter next;
      Iter prev;
      IndexValueMap value;
      ValueIndexMap id;
    };
    
    stack operator[](const bucket_type& i) { untested();
      assert(i < head.size());
      return stack(i, head.begin(), next.begin(), prev.begin(),
                   id_to_value.begin(), id);
    }
    stack operator[](const bucket_type& i) const{
      assert(i < head.size());
		typedef std::vector<size_type> vec;
		auto& h=const_cast<vec&>(head);
		auto& n=const_cast<vec&>(next);
		auto& p=const_cast<vec&>(prev);
		auto&i2v=const_cast<std::vector<value_type>&>(id_to_value);
//		auto& vim=const_cast<ValueIndexMap&>(id);

      return stack(i, h.begin(), n.begin(), p.begin(),
                   i2v.begin(), id);
    }
  protected:
    std::vector<size_type>   head;
    std::vector<size_type>   next;
    std::vector<size_type>   prev;
    std::vector<value_type>  id_to_value;
    Bucket bucket;
    ValueIndexMap id;
  };
  
}

#endif

// vim:ts=8:sw=2:et
