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
#ifndef TD_BUCKET_SORTER
#define TD_BUCKET_SORTER

#include <vector>
#include <cassert>
#include <boost/limits.hpp>

#include "random_generators.hpp"

namespace treedec {

  template <class BucketType, class ValueType, class Bucket, 
            class ValueIndexMap>
  class bucket_sorter {
  public:
    typedef BucketType bucket_type;
    typedef ValueType value_type;
    typedef typename std::vector<value_type>::size_type size_type;
    
    bucket_sorter(size_type _length, bucket_type _max_bucket, 
                  const Bucket& _bucket = Bucket(), 
                  const ValueIndexMap& _id = ValueIndexMap()) 
      : head(_max_bucket, invalid_value()),
        next(_length, invalid_value()), 
        prev(_length, invalid_value()),
        last(_max_bucket, invalid_value()),
        id_to_value(_length),
        bucket(_bucket), id(_id) { }
    
    void remove(const value_type& x) {
      const size_type i = boost::get(id, x);
      const size_type& next_node = next[i];
      const size_type& prev_node = prev[i];

      //check if i is the end of the bucket list 
      if ( next_node != invalid_value() )
        prev[next_node] = prev_node; 
      else{
        if(prev_node == invalid_value()){
        }
        else if(head[ bucket[x] ] != invalid_value()){
          last[ bucket[x] ] = prev_node;
        }
      }
      //check if i is the begin of the bucket list
      if ( prev_node != invalid_value() )
        next[prev_node] = next_node;
      else { //need update head of current bucket list
        head[ bucket[x] ] = next_node;
        last[ bucket[x] ] = next_node;
      }
    }

    void push(const value_type& x) {
      id_to_value[boost::get(id, x)] = x;
      assert(bucket[x] < head.size());
      (*this)[bucket[x]].push(x);
    }
    
    void update(const value_type& x) {
      remove(x);
      assert(bucket[x] < head.size());
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
      stack(bucket_type _bucket_id, Iter h, Iter n, Iter p, Iter l, IndexValueMap v,
            const ValueIndexMap& _id)
      : bucket_id(_bucket_id), head(h), next(n), prev(p), last(l), value(v), id(_id) {}

      // Avoid using default arg for ValueIndexMap so that the default
      // constructor of the ValueIndexMap is not required if not used.
      stack(bucket_type _bucket_id, Iter h, Iter n, Iter p, Iter l, IndexValueMap v)
        : bucket_id(_bucket_id), head(h), next(n), prev(p), value(v) {}
      
      void push(const value_type& x) {
        const size_type new_head = boost::get(id, x);
        const size_type current = head[bucket_id];

        if ( current != invalid_value() )
          prev[current] = new_head;
        else{
          last[bucket_id] = new_head;
        }
        prev[new_head] = invalid_value();
        next[new_head] = current;
        head[bucket_id] = new_head;
      }
      void pop_front() {
        size_type current = head[bucket_id];
        size_type next_node = next[current];
        head[bucket_id] = next_node;
        if ( next_node != invalid_value() )
          prev[next_node] = invalid_value();
        else{
          last[bucket_id] = invalid_value();
        }
      }
      void pop_back() {
        size_type current = last[bucket_id];
        if ( current != invalid_value() ){
          size_type prev_node = prev[current];
          if ( prev_node != invalid_value() ){
            next[prev_node] = invalid_value();
            last[bucket_id] = prev_node;
          }
          else{
              head[bucket_id] = invalid_value();
              last[bucket_id] = invalid_value();
          }
        }
        else{
          last[bucket_id] = invalid_value();
        }
      }
      std::pair<value_type, bool> next_value() {
        bool c = treedec::random::coin();
        //if(c){ std::cout << "coin: true" << std::endl; }
        //else { std::cout << "coin: false" << std::endl; }
        //std::cout << top() << std::endl;
        //std::cout << back() << std::endl;
        if(c) {
          assert(top() != invalid_value());
          return std::make_pair(top(), true);
        }
        else {
          assert(back() != invalid_value());
          return std::make_pair(back(), false);
        }
      }

      value_type& back() { assert(last[bucket_id] != invalid_value()); return value[ last[bucket_id] ]; }
      value_type& top() {  assert(head[bucket_id] != invalid_value()); return value[ head[bucket_id] ]; }
      const value_type& top() const { return value[ head[bucket_id] ]; }
      bool empty() const { return head[bucket_id] == invalid_value(); }
    private:
      bucket_type bucket_id;
      Iter head;
      Iter next;
      Iter prev;
      Iter last;
      IndexValueMap value;
      ValueIndexMap id;
    };
    
    stack operator[](const bucket_type& i) {
      //std::cout << "i: " << i << std::endl;
      //std::cout << "head.size(): " << head.size() << std::endl;
      assert(i < head.size());
      return stack(i, head.begin(), next.begin(), prev.begin(), last.begin(),
                   id_to_value.begin(), id);
    }
  protected:
    std::vector<size_type>   head;
    std::vector<size_type>   next;
    std::vector<size_type>   prev;
    std::vector<size_type>   last;
    std::vector<value_type>  id_to_value;
    Bucket bucket;
    ValueIndexMap id;
  };
  
}

#endif
