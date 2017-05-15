// Felix Salfelder 2016, 2017
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
//
// Tree decompositions
//
//
// a tree decomposition is a forest with bags at its vertices

#ifndef TREEDEC_HPP
#define TREEDEC_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <assert.h>
#include "graph.hpp"

namespace treedec{

// VECTOR_TD.
//
// conceptually, a boost directed graph with no more than one out-edge per
// vertex and bag" vertex property. does not allow remove_vertex on vertices
// other than the last.

template<class G>
class VECTOR_TD{
public: //types
    typedef typename boost::graph_traits<G>::vertex_descriptor VD;
    typedef std::vector<VD> bag_type;
    typedef struct value_type{
        value_type* first;
        bag_type second;
    } value_type;
    typedef value_type* vertex_descriptor;
    typedef value_type const* const_vertex_descriptor;
    class adjacency_iterator{
    public:
        adjacency_iterator(vertex_descriptor v):_v(v)
        {
        }
    public:
        vertex_descriptor operator*()
        {
            return _v;
        }
    private:
        vertex_descriptor _v;
    };
    typedef std::vector<value_type> container_t;
    typedef typename container_t::iterator iterator;
// is this really needed?
    class const_vertex_iterator : public container_t::const_iterator{
    public:
        const_vertex_iterator(typename container_t::const_iterator i)
            : container_t::const_iterator(i)
        {
        }
        const_vertex_iterator(typename container_t::iterator i)
            : container_t::const_iterator(i)
        {untested();
        }
    public:
        const_vertex_descriptor operator*()
        {
            return &container_t::const_iterator::operator*();
        }
    };
    class vertex_iterator : public container_t::iterator{
    public:
        vertex_iterator() : container_t::iterator()
        {
        }
        vertex_iterator(typename container_t::iterator i) : container_t::iterator(i)
        {
        }
        vertex_descriptor operator*()
        {
            return &container_t::iterator::operator*();
        }
    };
    typedef typename container_t::const_iterator const_iterator;
    typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;

public: // not implemented
    typedef void directed_category;
    typedef void edge_parallel_category;
    typedef void traversal_category;
public: //bug. should not have to be here.
    typedef boost::no_property vertex_property_type;
public: //construct
    VECTOR_TD(G const& g, unsigned bag_size=0) : _size(0), _bs(bag_size), _g(g)
    {
        trace1("results", boost::num_vertices(_g));
    }

    void reserve(size_t x){
        _v.reserve(x);
    }

public: //access
    size_t size() const
    {
        return _size;
    }
    vertex_iterator begin()
    {
        return _v.begin();
    }
    vertex_iterator end()
    {
        return _v.begin() + _size;
    }
    const_vertex_iterator begin() const
    {
        return const_vertex_iterator(_v.begin());
    }
    const_vertex_iterator end() const
    {
        return _v.begin() + _size;
    }
    void pop_back()
    {
        --_size;
    }
    void push_back(const value_type& x)
    { incomplete();
        assert(_size<=_v.size());
        if(_size==_v.size()){
            _v.push_back(x);
        }else{
            assert(_v.size()>_size);
            _v[_size] = x;
        }
        ++_size;
    }
    void erase(iterator b, iterator)
    {
       assert(_size<=_v.size());
       _size = b - _v.begin();
       assert(_size<=_v.size());
    }
    value_type& back()
    {
        assert(_size<=_v.size());
        assert(_size>0);
        return _v[_size-1];
    }
    // new self-loop-node
    vertex_descriptor new_one()
    {
        assert(_size<=_v.size());
        if(_size==_v.size()){
            _v.push_back(value_type());
            _v.back().second.reserve(_bs);
        }else{
            _v[_size].second.resize(0);
        }

        _v[_size].first = &_v[_size];
        ++_size;
        trace1("new_one", _size);
        assert(_size<=2*boost::num_vertices(_g));
        return &_v[_size-1];
    }
#if 0
    vertex_descriptor new_one(const_vertex_descriptor x)
    { untested();
        assert(_size<=_v.size());
        if(_size==_v.size()){
            _v.push_back(value_type());
        }else{
            _v[_size].second.resize(0);
        }

        _v[_size].first = const_cast<vertex_descriptor>(x);
        ++_size;
        trace1("new_one", _size);
        assert(_size<=2*boost::num_vertices(_g));
        return &_v[_size-1];
    }
#endif
    void clear()
    {
        _size=0;
    }
public: // ops
    value_type& operator[](size_t i)
    {
        assert(i<_size);
        assert(_size<=_v.size());
        return _v[i];
    }

private:
    container_t _v;
    unsigned _size;
    unsigned _bs; // bag size
    G const& _g;
};

} // treedec

namespace boost{
    template<class G>
    struct graph_traits<treedec::VECTOR_TD<G> >
    {
        typedef typename treedec::VECTOR_TD<G>::vertex_property_type
            type;
    };
    template<class G>
    std::pair<typename treedec::VECTOR_TD<G>::const_vertex_iterator,
              typename treedec::VECTOR_TD<G>::const_vertex_iterator>
    vertices(typename treedec::VECTOR_TD<G> const& t)
    {
        return std::make_pair(t.begin(), t.end());
    }
    template<class G>
    std::pair<typename treedec::VECTOR_TD<G>::vertex_iterator,
              typename treedec::VECTOR_TD<G>::vertex_iterator>
    vertices(typename treedec::VECTOR_TD<G> & t)
    {
        return std::make_pair(t.begin(), t.end());
    }
    // adjacent vertices as in "out_edges"
    template<class G>
    std::pair<typename treedec::VECTOR_TD<G>::adjacency_iterator,
              typename treedec::VECTOR_TD<G>::adjacency_iterator>
    adjacent_vertices(typename treedec::VECTOR_TD<G>::vertex_descriptor v,
                      typename treedec::VECTOR_TD<G> const&)
    {
        return std::make_pair(v->first, v->first+1);
    }
    // there's only one edge. replace it.
    template<class G>
        std::pair<
    std::pair<typename treedec::VECTOR_TD<G>::vertex_descriptor,
              typename treedec::VECTOR_TD<G>::vertex_descriptor>
                  ,bool >
    add_edge(typename treedec::VECTOR_TD<G>::vertex_descriptor v,
             typename treedec::VECTOR_TD<G>::const_vertex_descriptor b,
                      typename treedec::VECTOR_TD<G> &)
    {
        assert(v!=b);
        typedef typename treedec::VECTOR_TD<G>::vertex_descriptor vertex_descriptor;
        // should be OK, as the graph is mutable.
        vertex_descriptor B=const_cast<vertex_descriptor>(b);
        v->first = B;
        return std::make_pair( std::make_pair(v, B), true);
    }
    template<class G>
    bool degree( typename treedec::VECTOR_TD<G>::const_vertex_descriptor v,
                treedec::VECTOR_TD<G> const&)
    {
        return v!=v->first;
    }
    template<class G>
    unsigned num_vertices(treedec::VECTOR_TD<G> const& g)
    {
        return g.size();
    }
    template<class G>
    typename treedec::VECTOR_TD<G>::vertex_descriptor
        add_vertex(treedec::VECTOR_TD<G>& g)
    {
        return g.new_one();
    }
    template<class G>
    void remove_vertex(
            typename treedec::VECTOR_TD<G>::vertex_descriptor v,
            treedec::VECTOR_TD<G>& g)
    {
        // not possible to remove vertices other than top
        (void) v;
        assert(&g.back() == v);
        return g.pop_back();
    }

    template<class G>
    inline void
    get(vertex_index_t, treedec::VECTOR_TD<G> & g,
            typename treedec::VECTOR_TD<G>::vertex_descriptor v)
    //  boost::get(boost::vertex_index, G, v);
    { incomplete();
    }
} // boost

#endif
