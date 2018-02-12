// Felix Salfelder, 2014 - 2016
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

// graph related iterator classes
//
// !!!! does not fully work with c++<11 !!!!
//
// make_components_range
//     iterate through components of induced subgraph
// make_bfs_range
//     breadth-first iteration based on vertex range
// make_neighbourhood_range
//     iterate immediate neighbours with or without base range
// make_onion_range
//     iterate through "onion layers" around vertex range
//

#ifndef TREEDEC_ITER_HPP
#define TREEDEC_ITER_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp> // adjacent_vertices
#include <tr1/utility> // pair

#include "graph.hpp"

#ifndef NDEBUG
template<class S, class T>
size_t count_range(S i, T const& e)
{
    unsigned ii=0;
    for(;i!=e;++i) ++ii;
    return ii;
}
template<class P>
size_t count_range(P const& p)
{
    return count_range(p.first, p.second);
}
#else
template<class P>
size_t count_range(P const&)
{ untested();
    return 99;
}
template<class S, class T>
size_t count_range(S i, T const& e)
{ untested();
    return 88;
}
#endif

namespace util{ //

template<class I, class ref>
class skip_it : public I{//
public:
//	skip_it(const skip_it& p) : I(p)
//	{ untested();
//	}
	skip_it(I& i, I const& e, const ref& what): I(i), _skip(what), _end(e)
	{ untested();
		skip();
	}
	skip_it& operator++()
	{ untested();
		if (I(*this)==(_end)){ untested();
		}else{ untested();
			I::operator++();
			skip();
		}
		return *this;
	}
	bool operator==(const I& i) const
	{ untested();
		return I(*this)==(i);
	}
	bool operator!=(const I& i) const
	{ untested();
		return I(*this)!=(i);
	}

private:
	void skip()
	{ untested();
		if (I(*this)==(_end)){ untested();
		}else if (I::operator*() != _skip){ untested();
		}else{ untested();
			I::operator++();
		}
	}
	ref _skip;
	I _end;
};

template<class I, class S>
inline skip_it<I,S> make_skip_it(I b, I e, S s)
{ untested();
	return skip_it<I,S>(b, e, s);
}

template<class C, class B>
inline skip_it<C,B> make_skip_it(C c, B b)
{untested();
	return make_skip_it(c.begin(), c.end() ,b);
}

} // sethack

template<class A>
void debug_count(A i, A e, unsigned s)
{ untested();
    unsigned c=0;
    for(; i!=e; ++i){ untested();
        ++c;
    }
    assert(c==s);
}

typedef enum{
    bEXCLUDE=-1,
    bPASS=0,
    bINCLUDE=1
}base_t;

namespace detail{

struct localmask{
    typedef std::vector<BOOL> S;

    localmask() {}
    localmask(unsigned x) : _s(x) {}
    bool operator()(unsigned x) const{
        return _s[x];
    }
    bool operator[](unsigned x) const{
        return operator()(x);
    }
    void visit (unsigned x){
        _s[x] = true;
    }
    void resize(size_t t){ _s.resize(t); }
    size_t size() const{ return _s.size(); }
    S _s;
};

template<class M>
struct mask_help{
    template<class A, class B>
    static void init(A, B){
    }
};

template<>
struct mask_help<localmask>{
    template<class A, class B>
    static void init(A& a, B size){
        a.resize(size);
    }
};

// iterate through the connected components of a graph
// THIS IMPLEMENTS A ONE-PASS ITERATOR
template<class G, class VR, class vis_t=localmask >
class components_iter{ //
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    typedef typename VR::first_type some_iterator;
    typedef typename std::pair<adjacency_iterator, adjacency_iterator> adj_range;
    typedef typename std::vector<adj_range> stack_type;

    typedef typename std::vector<adj_range> scratch_type;
//    typedef std::vector<BOOL_> vis_t;
    class component_iter{ //
    public:
        typedef component_iter hackself;
    public:
        component_iter(typename VR::first_type v, components_iter& cs, bool visit=true)
           : _v(v) /* necessary? */, _cs(cs)
        {
            assert(cs._stack.size()==0);
            if(v==_cs._range.second){
            }else{
                auto pos=boost::get(boost::vertex_index, _cs._g, *_v);
                if(visit){
                    _cs._visited.visit(pos);
                }else{
                }
            }
        }
    private:
        bool is_end() const
        {
            return(_v==_cs._range.second);
        }
    public: // ops
        vertex_descriptor operator*() {
            assert(!is_end());
            vertex_descriptor v;
            if(_cs._stack.size()){
                v = *(_cs._stack.back().first);
            }else{
                assert(_v!=_cs._range.second);
                v = *_v;
            }
            assert(treedec::is_valid(v, _cs._g));
            auto pos=boost::get(boost::vertex_index, _cs._g, v);
            (void)pos;
            assert(_cs._visited[pos]);
            return v;
        }
        component_iter& operator++() {
            vertex_descriptor pre=**this;
            size_t pos=boost::get(boost::vertex_index, _cs._g, pre);
            (void)pos;
            // trace1("cmp++ ", pos);
            assert(_cs._visited[pos]);

            BOOST_AUTO(p, boost::adjacent_vertices(pre, _cs._g));
            _cs._stack.push_back(p);

            // DFS
            while(true){
                BOOST_AUTO(& candi, _cs._stack.back());
                if(candi.first==candi.second){
                    _cs._stack.pop_back();

                    if(_cs._stack.size()){
                        // incomplete?
                    }else{
                        assert(treedec::is_valid(*_v, _cs._g));
                        _v = _cs._range.second; // indicate eoc

                        if(_v!=_cs._range.second){
                            // auto pos=treedec::get_pos(*_v, _cs._g);
                        //    _cs._visited[pos] = true;
                        }else{
                        }
                        return *this;
                    }
                }else{
                    unsigned pos=boost::get(boost::vertex_index, _cs._g, *candi.first);
                    if (!_cs._visited[pos]){
                        // trace1("visiting", pos);
                        _cs._visited.visit(pos);
                        break;
                    }else{
                    }

                }

                next_nonvisited_or_end(_cs._stack.back());
            }

            assert(is_end() || _cs._visited[boost::get(boost::vertex_index, _cs._g, **this)]);

            return *this;
        }
#if 1
        bool operator<(const component_iter& other) const
        { incomplete();
            // uuh vector?
            if(&_v < &other._v){ untested();
            }else{ untested();
            }
            return &_v < &other._v;
        }
        bool operator!=(const component_iter& other) const
        {
            if(_v != other._v){
            }else{
            }
            return _v != other._v;
        }
#else
        bool operator!=(const component_iter& other) const
        {
            if(_cs._stack.size()){
                return _cs._stack.back() != other._cs._stack.back();
            }else{
                return true;
            }
        }
        bool operator==(const component_iter& other) const
        {
            if(_cs._stack.size()){
                return _cs._stack.back() == other._cs._stack.back();
            }else{
                return false;
            }
        }
#endif

    private:
        template<class T>
        // component_iter<T>::
        void next_nonvisited_or_end(T& p)
        {
            while(true){
                if(p.first==p.second){
                    break;
                }else{
                }

                BOOST_AUTO(f, *p.first);
                unsigned pos=boost::get(boost::vertex_index, _cs._g, f);
                assert(pos<_cs._visited.size());
                if(!_cs._visited[pos]){
             //       _cs._visited[pos] = true;
                    return;
                }else{
                    // trace1("visited already", pos);
                    ++p.first;
                }
            }
        }
    private:
        some_iterator _v; // use _cs?
        //vertex_iterator _e;
        components_iter& _cs;
    }; // component_iter
public: // construct
//    components_iter(vertex_iterator v, G const& g)
//        : _base(v),
//          _g(g),
//          _stack(0)
//    {
//        _visited.assign(boost::num_vertices(_g), false);
//    }
   components_iter(VR v, G const& g, vis_t vis, scratch_type* c=NULL)
       : _range(v),
         _visited(vis),
         _stack(c?(*c):(*(new scratch_type(0)))), _cc(&_stack),
         _g(g)
   {
       if(c){
           _cc = NULL;
           _stack.resize(0);
       }else{
       }
       // trace3("cmpsiter", *v.first, visited, v.first==v.second);
   }

   // needed for end?
   components_iter(VR v, G const& g, scratch_type* c=NULL)
       : _range(v),
         _visited(vis_t(0)),
         _stack(c?(*c):(*(new scratch_type(0)))), _cc(&_stack),
         _g(g)
   {
       if(c){
           _cc = NULL;
           _stack.resize(0);
       }else{
       }
   }
   // hack
   components_iter(VR v, G const& g, G const&, scratch_type* c=NULL )
       : _range(v),
         _visited(vis_t(boost::num_vertices(g))),
         _stack(c?(*c):(*(new scratch_type(0)))), _cc(&_stack),
         _g(g)
   {
       if(c){
           _cc = NULL;
           _stack.resize(0);
       }else{
       }
   }
public: // copy
   components_iter(components_iter const& p)
       : _range(p._range),
         _visited(vis_t(p._visited)),
         _stack(*(new scratch_type(p._stack))), _cc(&_stack),
         _g(p._g)
   {
   }
   components_iter(components_iter const&& p)
       : _range(p._range),
         _visited(p._visited),
         _stack(p._stack), _cc(p._cc),
         _g(p._g)
   {
       p._cc = NULL;
   }
   ~components_iter()
   {
       if(_cc){
           assert(_cc=&_stack);
           delete(&_stack);
       }else{
       }
   }
public: // ass
    components_iter& operator=(const components_iter&)
    { unreachable();
        incomplete();
        return *this;
    }
#if __cplusplus >= 201103L
    components_iter& operator=(const components_iter&&)
    { unreachable();
        incomplete();
        return *this;
    }
#endif
public: // ops
    bool operator==(const vertex_iterator& other)
    {
        return _range.first == other;
    }
    bool operator<(const vertex_iterator& other)
    { incomplete();
        if(_range.first < other){ untested();
        }else{
        }
        return _range.first < other;
    }
    bool operator!=(const vertex_iterator& other)
    {
        if(_range.first != other){
        }else{
        }
        return _range.first != other;
    }
    bool operator==(const components_iter& other) const
    {
        return _range.first == other._range.first;
    }
    bool operator<(const components_iter&) const
    { incomplete();
        return false;
        // uuh ooh. vector?! (does not work)
//        if( &*_range.first < &*other._range.first){ untested();
//        }else{ untested();
//        }
//        return &*_range.first < &*other._range.first;
    }
    bool operator!=(const components_iter& other) const
    {
        if( _range.first != other._range.first){
        }else{
        }
        return _range.first != other._range.first;
    }
    components_iter& operator++()
    {
        unsigned p=-1;
        while(true){
            if(_range.first==_range.second){
                break;
            }
            p = boost::get(boost::vertex_index, _g, *_range.first);
            // trace1("cmps++ ", p);
            assert(p<_visited.size());
            if(_range.first==_range.second){
                break;
            }else if(_visited[p]){
                assert(_range.first!=_range.second);
                ++_range.first;
            }else{
                break;
            }
        }

        assert(_range.first==_range.second || !_visited[p]);
        return *this;
    }
    std::pair<component_iter, component_iter> operator*()
    {
        BOOST_AUTO(p, boost::get(boost::vertex_index, _g, *_range.first));
        (void)p;
        // trace2("op*", *_range.first,  _visited[p]);
        assert(_range.first==_range.second || !_visited[p]);

        BOOST_AUTO(second, _range.first);
        if(_range.first!=_range.second){
            ++second;
        }else{
        }
        return std::make_pair(component_iter(_range.first, *this),
                              component_iter(_range.second, *this, false));
    }
private: // state
    VR _range;
    vis_t _visited;
//    mutable vis_t* _vv; // bool?
    stack_type& _stack;
    mutable stack_type* _cc; // bool?
    G const& _g;
}; // components_iter

} // detail


#define COMMA ,
#define VRP_ std::pair< \
typename boost::graph_traits<G>::vertex_iterator COMMA \
typename boost::graph_traits<G>::vertex_iterator >
template<class G, class MASK=detail::localmask >
std::pair<detail::components_iter<G, VRP_, MASK>,
          detail::components_iter<G, VRP_, MASK> >
    make_components_range(G const& g,
        typename detail::components_iter<G, VRP_, MASK>::scratch_type* s=NULL,
            BOOL b=true
        )
{
    (void) b;
    BOOST_AUTO(p, boost::vertices(g));
    return std::make_pair(
        detail::components_iter<G, VRP_, MASK>(p, g, g, s),
        detail::components_iter<G, VRP_, MASK>(std::make_pair(p.second, p.second), g)); // FIXME
}
#undef VRP_

#define VRIP_ std::pair< VRI COMMA VRI >
#define VRIEP_ std::pair< VRI COMMA VRI >

// iterate through components of graph induced by non-masked nodes that
// intersect with a vertex range *[i, e).
template<class G, class VRI, /*class VRIE,*/ class MASK=detail::localmask>
std::pair<detail::components_iter<G, VRIP_, MASK>,
          detail::components_iter<G, VRIEP_, MASK> >
    make_components_range(
            VRI i, VRI e,
            G const& g,
            MASK mask=detail::localmask(),
            typename detail::components_iter<G, VRIP_, MASK>::scratch_type* s=NULL)
{
    detail::mask_help<MASK>::init(mask, boost::num_vertices(g));
    assert(i==e || treedec::is_valid(*i, g));
    while( i!=e ){
        assert(treedec::is_valid(*i, g));
        if(mask[boost::get(boost::vertex_index, g, *i)]){
            ++i;
        }else{
            break;
        }
    }
    if(i==e){
    }else{
        assert(!mask[boost::get(boost::vertex_index, g, *i)]);
    }
    return std::make_pair(
        detail::components_iter<G, VRIP_, MASK>(std::make_pair(i, e), g, mask, s),
        detail::components_iter<G, VRIEP_, MASK>(std::make_pair(e, e), g, mask));
}
#undef VRIP_
#undef VRIEP_

namespace detail{

// sort of bfs
// THIS IMPLEMENTS A ONE-PASS ITERATOR
template<class G, class VR, class BOOL_=bool>
class bfs_iter{ //
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    typedef typename VR::first_type some_iterator;
    typedef typename std::pair<adjacency_iterator, adjacency_iterator> adj_range;
    typedef typename std::deque<adj_range> queue_type;
    typedef std::pair<bfs_iter, bfs_iter> bfs_range;

    typedef queue_type scratch_type;
    typedef std::vector<BOOL_> vis_t;
    class bfs_end{ //
    };


   bfs_iter(VR v, G const& g, vis_t* vis, scratch_type* c=NULL)
       : _visited(vis?(*vis):(*(new vis_t(0)))), _vv(&_visited),
         _q(c?(*c):(*(new scratch_type(0)))), _qq(&_q),
         _g(g)
   {
       if(!vis){ incomplete();
       }else{
           // use external
           _vv = NULL;
       }
       if(c){
           _qq = NULL;
           _q.resize(0);
       }else{
       }

       trace1("new bfs", count_range(v));
       for(; v.first!=v.second; ++v.first){
           trace2("base node", *v.first, boost::out_degree(*v.first, _g));
           maybe_push_back(boost::adjacent_vertices(*v.first, _g));
       }

       skip_and_visit();
       trace1("new bfs done", _q.size());
   }
private:
    // skip to next nonvisited
    void skip_and_visit()
    {
        while(_q.begin() != _q.end()){
            trace1("..", count_range(front_range()));
            if(front_range().first != front_range().second){
    //            ++front_range().first;
                trace2("nnm", count_range(front_range()), _q.size());
                next_nonvisited_or_end();
                if(front_range().first != front_range().second){
                    return;
                }
            }else{
            }

            if (front_range().first == front_range().second ){
                pop_q();
            }else if (front_range().first != front_range().second ){ untested();
                trace2("?!", count_range(front_range()), _q.size());
                ++_q.front().first;
            }else{ untested();
            }
        }
    }
    void pop_q(){
        _q.pop_front();
        //if(_layerend == _q.begin()){ untested();
        //    ++_layer;
        //    _layerend = _q.end();
        //}else{ untested();
        //}
    }
public: // ops
    vertex_descriptor operator*()
    {
        if(_q.front().first!=_q.front().second){
        }else{ untested();
        }
        BOOST_AUTO(vd, *front_range().first);
        trace1("op*", vd);
        return vd;
    }
    bfs_iter& operator++()
    {
        assert(_q.begin() != _q.end());
        ++front_range().first;
        if(front_range().first == front_range().second){
        }else{
        }
//            trace1("++ in layer", count_range(front_range()));
        skip_and_visit();
        // trace1("++'d in layer", count_range(front_range()));
        return *this;
    }
    bool operator==(const bfs_iter&e) const
    { untested();
        return !operator!=(e);
    }
    bool operator!=(const bfs_iter&) const
    {
        if(_q.empty()){
            trace1("layer atend", _q.size());
            return false;
        }else if(_q.front().first!=_q.front().second){
            return true;
        }else{ untested();
            return false;
        }
    }

private:
    void set_visited(unsigned pos)
    {
        _visited[pos] = true;
    }
    BOOL visited(unsigned pos) const
    {
        assert(pos<_visited.size());
        return _visited[pos];
    }
    typename queue_type::value_type const& front_range() const
    { untested();
        assert(_q.begin() != _q.end());
        return _q.front();
    }
    typename queue_type::value_type& front_range()
    {
        assert(_q.begin() != _q.end());
        return _q.front();
    }

    void next_nonvisited_or_end() {
        trace2("next_nonvisited_or_end", _q.size(), count_range(front_range()));
        while(front_range().first!=front_range().second){
            vertex_descriptor v=*front_range().first;
            auto pos=boost::get(boost::vertex_index, _g, v);
            if(visited(pos)){
                trace1("been there", pos);
            }else{
                set_visited(pos);
                BOOST_AUTO(av, boost::adjacent_vertices(v, _g));
                assert(av.first!=av.second);
                maybe_push_back(av);
                trace2("found new", pos, count_range(av));
                trace2("found new", pos, count_range(_q.back()));
                assert(count_range(av)==boost::out_degree(v, _g));
                return;
            }
            ++front_range().first;
        }

        if(_q.empty()){ untested();
        }else if(front_range().first==front_range().second){
//                _q.pop_front();
        }else{ untested();
        }
    } // next_nonvisited_or_end
private:
    // typedef std::pair<layer_iter, layer_end> onion_layer_range;
public: // construct
//    onion_iter(vertex_iterator v, G const& g)
//        : _base(v),
//          _g(g),
//          _q(0)
//    { untested();
//        _visited.assign(boost::num_vertices(_g), false);
//    }

   template<class T>
   void maybe_push_back(T x)
   {
       if(x.first==x.second){ untested();
           //reachable?!
           return;
       }else{
       }
       while(true){
           auto pos=boost::get(boost::vertex_index, _g, *x.first);
           if(visited(pos)){
               ++x.first;
           }else{
               trace1("push", count_range(x));
               _q.push_back(x);
               break;
           }
           if(x.first==x.second){
               break;
           }
       }
   }

   // needed for end?
   bfs_iter(G const& g, scratch_type* c=NULL)
       : _visited(*(new vis_t(0))), _vv(&_visited),
         _q(c?(*c):(*(new scratch_type(0)))), _qq(&_q),
         _g(g)
   {
       if(c){ untested();
           _qq = NULL;
           _q.resize(0);
       }else{
       }
   }
public: // copy
#if __cplusplus >= 201103L
   bfs_iter(G const& g, int /*dummy*/)
       : _visited(*(new vis_t(0))), _vv(&_visited),
         _q(*(new scratch_type())), _qq(&_q),
         _g(g)
   {
   }
   bfs_iter(bfs_iter const& p)
       : _visited(*(new vis_t(p._visited))), _vv(&_visited),
         _q(*(new scratch_type(p._q))), _qq(&_q),
         _g(p._g)
   { untested();
       incomplete();
   }
   bfs_iter(bfs_iter const&& p)
       : _visited(p._visited), _vv(p._vv),
         _q(p._q), _qq(p._qq),
         _g(p._g)
   {
       p._qq = NULL;
       p._vv = NULL;
   }
public: // assign
    bfs_iter& operator=(const bfs_iter& other)
    { untested();
        (void) other;
        assert(0);
        assert(&_g == &other._g);
        return *this;
    }
    bfs_iter& operator=(const bfs_iter&& p)
    {
        // if(_vv){ untested();
        //     delete(_vv);
        // }
        // if(_qq){ untested();
        //     delete(_qq);
        // }
        _q = MOVE(p._q);
        trace1("hmm1", _visited.size());
        trace1("hmm1", p._visited.size());
        _visited=MOVE(p._visited);
        trace1("_hmm", _visited.size());

        (void) p;
        assert(&_g == &p._g);
        p._qq = NULL;
        p._vv = NULL;
        return *this;
    }
public: // destroy
#endif
   ~bfs_iter()
   {
       if(_vv){
           // incomplete();
           assert(&_visited==_vv);
           delete(_vv);
       }else{
       }
       if(_qq){
           assert(_qq==&_q);
           delete(_qq);
       }else{
       }
   }
private: // state
#ifndef NDEBUG
public:
#endif
    vis_t& _visited;
    mutable vis_t* _vv; // bool?
    queue_type& _q;
    mutable queue_type* _qq; // bool?
    G const& _g;
}; // bfs_iter

} // detail

#define ORIP_ std::pair< VRI COMMA VRI >

// typedef typename onion_iter<G, ORIP_, BOOL>::onion_range onion_range;
// iterate through layers of graph induced by non-masked nodes connected to
// vertex range *[i, e).
template<class G, class VRI, /*class VRIE,*/ class BOOL_>
    typename detail::bfs_iter<G, ORIP_, BOOL_>::bfs_range
    make_bfs_range(
            VRI i, VRI e,
            G const& g,
            std::vector<BOOL_>* mask,
            typename detail::bfs_iter<G, ORIP_, BOOL_>::scratch_type* s=NULL)
{
    return std::make_pair(
        detail::bfs_iter<G, ORIP_, BOOL_>(std::make_pair(i, e), g, mask, s),
        detail::bfs_iter<G, ORIP_, BOOL_>(g));
}
#undef ORIP_

namespace detail{

// sort of bfs
// THIS IMPLEMENTS A ONE-PASS ITERATOR
template<class G, class VR, class BOOL_=bool>
class onion_iter{ //
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    typedef typename VR::first_type some_iterator;
    typedef typename std::pair<adjacency_iterator, adjacency_iterator> adj_range;
    typedef typename std::deque<adj_range> queue_type;
    typedef std::pair<onion_iter, onion_iter> onion_range;

    typedef queue_type scratch_type;
    typedef std::vector<BOOL_> vis_t;
    class layer_end{ //
    };
    class layer_iter {
    public:
        layer_iter(onion_iter& o)
           : _onion(o)
        { untested();
            assert(_onion._q.begin() != _onion._layerend);
            assert(_onion._q.begin() != _onion._q.end());

            skip_and_visit();
            trace0("done constr");
        }
    private:
        // skip to next nonvisited
        void skip_and_visit()
        { untested();
            trace3("skip_and_visit", _onion._layerend - _onion._q.begin(),  _onion._q.size(),
                   count_range(front_range())  );
            assert(_onion._q.begin() != _onion._layerend);
            while(_onion._q.begin() != _onion._layerend){ untested();
                trace0("..");
                trace1("..", count_range(front_range()));
                if(front_range().first != front_range().second){ untested();
        //            ++front_range().first;
                    trace2("nnm", count_range(front_range()), _onion._q.size());
                    next_nonvisited_or_end();
                    if(front_range().first != front_range().second){ untested();
                        return;
                    }
                }else{ untested();
                }

                if(_onion._q.begin() == _onion._layerend){ untested();
                    // don't set visited!
                    return;
                }else if (front_range().first == front_range().second ){ untested();
                    _onion._q.pop_front();
                }else if (front_range().first != front_range().second ){ untested();
                    trace2("?!", count_range(front_range()), _onion._q.size());
                    ++_onion._q.front().first;
                }else{ untested();
                }
            }
        }
    public: // ops
        vertex_descriptor operator*()
        { untested();
            if(_onion._q.front().first!=_onion._q.front().second){ untested();
            }else{ untested();
            }
            BOOST_AUTO(vd, *front_range().first);
            trace1("op*", vd);
            return vd;
        }
        layer_iter& operator++()
        { untested();
            assert(_onion._q.begin() != _onion._layerend);
            assert(_onion._q.begin() != _onion._q.end());
            ++front_range().first;
            if(front_range().first == front_range().second){ untested();
            }else{untested();
            }
//            trace1("++ in layer", count_range(front_range()));
            skip_and_visit();
            // trace1("++'d in layer", count_range(front_range()));
            return *this;
        }
        bool operator!=(const layer_iter&) const { incomplete();
            assert(false);
            return true; // dont know how to do that.
        }
        bool operator==(const layer_end&e) const { untested();
            return !operator!=(e);
        }
        bool operator!=(const layer_end&) const
        { untested();
            if(_onion._q.empty()){ untested();
                trace1("layer atend", _onion._q.size());
                return false;
            }else if(_onion._q.begin()==_onion._layerend){ untested();
                return false;
            }else if(_onion._q.front().first!=_onion._q.front().second){ untested();
                return true;
            }else{ untested();
                return false;
            }
        }

    private:
        BOOL& visited(unsigned pos) { untested();
            assert(pos<_onion._visited.size());
            return _onion._visited[pos];
        }
        typename queue_type::value_type const& front_range() const { untested();
            assert(_onion._q.begin() != _onion._q.end());
            assert(_onion._q.begin() != _onion._layerend);
            return _onion._q.front();
        }
        typename queue_type::value_type& front_range() { untested();
            assert(_onion._q.begin() != _onion._q.end());
            assert(_onion._q.begin() != _onion._layerend);
            return _onion._q.front();
        }
        // onion_iter::
        void next_nonvisited_or_end() { untested();
            while(front_range().first!=front_range().second){ untested();
                vertex_descriptor v=*front_range().first;
                auto pos=boost::get(boost::vertex_index, _onion._g, v);
                if(visited(pos)){ untested();
                }else{ untested();
                    visited(pos) = true;
                    BOOST_AUTO(av, boost::adjacent_vertices(v, _onion._g));
                    assert(av.first!=av.second);
                    _onion.maybe_push_back(av);
                    assert(count_range(av)==boost::out_degree(v, _onion._g));
                    return;
                }
                ++front_range().first;
            } untested();

            if(_onion._q.empty()){ untested();
            }else if(front_range().first==front_range().second){ untested();
                assert(_onion._q.begin() != _onion._layerend);
            }else{ untested();
            }
        }
    private:
        // some_iterator _layerend; // use _cs?
        //vertex_iterator _e;
        onion_iter& _onion;
    }; // layer_iter
    typedef std::pair<layer_iter, layer_end> onion_layer_range;
public: // construct
//    onion_iter(vertex_iterator v, G const& g)
//        : _base(v),
//          _g(g),
//          _q(0)
//    { untested();
//        _visited.assign(boost::num_vertices(_g), false);
//    }
   onion_iter(VR v, G const& g, vis_t* vis, scratch_type* c=NULL)
       : _visited(vis?(*vis):(*(new vis_t(0)))), _vv(&_visited),
         _q(c?(*c):(*(new scratch_type(0)))), _qq(&_q),
         _g(g)
   { untested();
       if(!vis){ incomplete();
       }else{ untested();
           // use external
           _vv = NULL;
       }
       if(c){ untested();
           _qq = NULL;
           _q.resize(0);
       }else{ incomplete();
       }

       trace1("new onion", count_range(v));
       for(; v.first!=v.second; ++v.first){ untested();
           maybe_push_back(boost::adjacent_vertices(*v.first, _g));
       }

       // trace3("cmpsiter", *v.first, visited, v.first==v.second);
   }

   template<class T>
   void maybe_push_back(T x)
   { untested();
       assert(x.first!=x.second);
       while(true){ untested();
           auto pos=boost::get(boost::vertex_index, _g, *x.first);
           if(_visited[pos]){ untested();
               ++x.first;
           }else{ untested();
               trace0("push");
               _q.push_back(x);
               break;
           }
           if(x.first==x.second){ untested();
               break;
           }
       }
   }

   // needed for end?
   onion_iter(G const& g, scratch_type* c=NULL)
       : _visited(*(new vis_t(0))), _vv(&_visited),
         _q(c?(*c):(*(new scratch_type(0)))), _qq(&_q),
         _g(g)
   { untested();
       if(c){ untested();
           _qq = NULL;
           _q.resize(0);
       }else{ untested();
       }
   }
public: // copy
#if __cplusplus >= 201103L
   onion_iter(G const& g)
       : _visited(*(new vis_t(0))), _vv(&_visited),
         _q(*(new scratch_type())), _qq(&_q),
         _g(g)
   { untested();
   }
   onion_iter(onion_iter const& p)
       : _visited(*(new vis_t(p._visited))), _vv(&_visited),
         _q(*(new scratch_type(p._q))), _qq(&_q),
         _g(p._g)
   { untested();
   }
   onion_iter(onion_iter const&& p)
       : _visited(p._visited), _vv(p._vv),
         _q(p._q), _qq(p._qq),
         _g(p._g)
   { untested();
       p._qq = NULL;
       p._vv = NULL;
   }
#endif
   ~onion_iter()
   { untested();
       if(_vv){ untested();
           // incomplete();
           assert(&_visited==_vv);
           delete(_vv);
       }else{ untested();
       }
       if(_qq){ untested();
           assert(_qq==&_q);
           delete(_qq);
       }else{ untested();
       }
   }
public: // ass
    onion_iter& operator=(const onion_iter& other)
    { untested();
        (void) other;
        incomplete();
        assert(&_g == &other._g);
        return *this;
    }
#if __cplusplus >= 201103L
    onion_iter& operator=(const onion_iter&& other)
    { untested();
        (void) other;
        incomplete();
        assert(&_g == &other._g);
        return *this;
    }
#endif
public: // ops
    bool operator==(const onion_iter&) const
    { untested();
        // uuh ooh. just check if at end.
        return _q.begin() == _q.end();
    }
    bool operator!=(const onion_iter&) const
    { untested();
        trace1("not end?", _q.size());
        return _q.size();
    }
    onion_iter& operator++()
    { untested();
        trace1("next layer", _q.size());
        if(_q.empty()){untested();
            // happens if the last layer is empty but had iterators.
        }else{ untested();
        }
        assert(_q.begin()==_layerend); // one pass...
        return *this;
    }
    std::pair<layer_iter, layer_end> operator*()
    { untested();
        trace1("creating layeriter", _q.size());
        _layerend=_q.end();
        return std::make_pair(layer_iter(*this),
                              layer_end());
    }
private: // state
    VR _range;
    vis_t& _visited;
    mutable vis_t* _vv; // bool?
    queue_type& _q;
    mutable queue_type* _qq; // bool?
    G const& _g;
    typename queue_type::const_iterator _layerend;
}; // onion_iter

} // detail
#define ORIP_ std::pair< VRI COMMA VRI >

// typedef typename onion_iter<G, ORIP_, BOOL>::onion_range onion_range;
// iterate through layers of graph induced by non-masked nodes connected to
// vertex range *[i, e).
template<class G, class VRI, /*class VRIE,*/ class BOOL_>
    typename detail::onion_iter<G, ORIP_, BOOL_>::onion_range
    make_onion_range(
            VRI i, VRI e,
            G const& g,
            std::vector<BOOL_>* mask,
            typename detail::onion_iter<G, ORIP_, BOOL_>::scratch_type* s=NULL)
{ untested();
    return std::make_pair(
        detail::onion_iter<G, ORIP_, BOOL_>(std::make_pair(i, e), g, mask, s),
        detail::onion_iter<G, ORIP_, BOOL_>(g));
}
#undef ORIP_
#undef COMMA

template<class G, class VR, class BOOL>
inline typename detail::onion_iter<G, VR, BOOL>::scratch_type*
new_onion_range_scratch(G const&, VR, BOOL, unsigned size=0)
{ untested();
    return new typename detail::onion_iter<G, VR, BOOL>::scratch_type(size);
}

template<class G, class VR, class BOOL>
inline typename detail::onion_iter<G, VR, BOOL>::scratch_type*
new_bfs_range_scratch(G const&, VR, BOOL, unsigned size=0)
{
    return new typename detail::onion_iter<G, VR, BOOL>::scratch_type(size);
}

// FIXME: namespace detail

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
       : _tt(NULL),
         _t(t?(*t):(*(new scratch_type))),
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
      if(!_t.size()){
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
   {
      return !operator==(other._i);
   }
   subsets_iter operator++()
   {
      if(_t.size()==0){
         _t.push_back(_i);
//         assert(_i!=_e); why not?!
         if(_u==0){
            _t.back()=_e;
         }else{
         }
      }else if(_t.size()<=_u){
         BOOST_AUTO(back, _t.back());
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
            }
         }
      }else if(_t.back() != _e){ untested();
         incomplete();
      }else{ untested();
      }

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
            if(back==_e){
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
    make_subsets_range(A a, A b, unsigned l, unsigned u,
          typename subsets_iter<A>::scratch_type* s=NULL)
{
   return std::make_pair(
       subsets_iter<A>(a,b,l,u,s),
       subsets_iter<A>(b,b));
}

namespace detail{
// iterate neighbourhood of C, including C
// only works on ordered containers...
template<class A, class G>
class neighbourhood01_iter { //
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor value_type;
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef std::vector<adjacency_iterator> scratch_type;
public: // construct
//     neighbourhood01_iter(const neighbourhood01_iter& p)
//     { untested();
//      incomplete. reuse _a?
//     }
    neighbourhood01_iter() : _v(NULL),
        _a(*(new scratch_type(0))), _aa(&_a)
    { untested();
    }
    neighbourhood01_iter(A b, A e, unsigned size, const G& g, base_t incb,
            scratch_type* a=NULL)
       : _b(b), _i(b), _e(e),
         _a(a?(*a):(*(new scratch_type(size)))),
         _aa(&_a),
         _g(g), _include_base(incb)
    {
#if __cplusplus >= 201103L
        if(a){ untested();
            a->resize(size);
            // use external scratch
            _aa = NULL;
        }else{
            // own scratch. initialized above
            assert(_aa == &_a);
        }
#else
        // ignore external scratch
        assert(_aa == &_a);
#endif
        if(b==e){
            assert(size==0);
            return;
        }else{
        }

        bool found=false;

        if(_include_base){
            _v = *_i;
        }else{ untested();
            A ii(_b);
            for(; ii!=_e; ++ii){ untested();
                // just take one.... (inefficient)
                if(boost::adjacent_vertices(*ii, g).first
                 !=boost::adjacent_vertices(*ii, g).second){ untested();
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
            }else{
                _a[n] = boost::adjacent_vertices(*ii, g).first;
            }

            if(_a[n] == boost::adjacent_vertices(*ii, g).second){
            }else if(*(_a[n])<_v){
                _v = *_a[n];
                found = true;
            }else{
            }
            ++n;
        }
        if(_include_base==bEXCLUDE){ untested();
            incomplete();
        }else if(_include_base==bINCLUDE){
        }else if(!found){ untested();
            _b = _e;
        }
    }
    neighbourhood01_iter& operator=(const neighbourhood01_iter&)
    { untested();
        incomplete();
        return *this;
    }

    neighbourhood01_iter(const neighbourhood01_iter& x)
        : _b(x._b), _i(x._i), _e(x._e),
          _a(*(new scratch_type(x._a))), _aa(&_a),
          _v(x._v),
          _g(x._g), _include_base(x._include_base)
    {
        assert(_aa == &_a);
    }
#if __cplusplus >= 201103L
    neighbourhood01_iter(const neighbourhood01_iter&& x)
        : _b(x._b),
          _i(x._i),
          _e(x._e),
          _a(x._a),
          _aa(x._aa),
          _v((x._v)),
          _g(x._g),
          _include_base(x._include_base)
    {
        x._aa = NULL;
        assert(_aa == NULL || _aa == &_a);
    }
    neighbourhood01_iter& operator=(const neighbourhood01_iter&&)
    { untested();
        incomplete();
        return *this;
    }
#endif

    ~neighbourhood01_iter()
    {
        if(_aa){
            // own scratch. delete it.
            delete(&_a);
        }else{
        }
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
        if(_include_base){
            found = update(_i, _e, previous, _v);
        }else{ untested();
        }
        A ii(_b);
        unsigned n=0;
        for(; ii!=_e; ++ii){
            BOOST_AUTO(aend, boost::adjacent_vertices(*ii, _g).second);
            found |= update(_a[n], aend, previous, _v);
            ++n;
        }
        if(!found){
            _b = _e;
        }else{
            assert(_v>previous);
        }
        return *this;
    }
    const vertex_descriptor& operator*()
    {
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
    scratch_type& _a;
    mutable scratch_type* _aa;
    vertex_descriptor _v;
    G const& _g;
    base_t _include_base;
};
} // detail

namespace std{
    template<class A, class B>
    struct iterator_traits<detail::neighbourhood01_iter<A, B> >{ //
        typedef typename detail::neighbourhood01_iter<A, B>::value_type value_type;
    };
}

// iterate depending on include_base
// - [b, e) and neighbors           if bINCLUDE
// - only neighbors                 if bPASS
// - neighbors without [b, e)       if bEXCLUDE
// neighbors(S)={v\in V | \exist s\in S : E(s,v)}

// 01 variant obsolete?
template<class A, class G>
std::pair<detail::neighbourhood01_iter<A, G>,
          detail::neighbourhood01_iter<A, G> >
inline make_neighbourhood01_range(A b, A e, G const& g, unsigned size=0,
        base_t include_base=bINCLUDE,
        typename detail::neighbourhood01_iter<A,G>::scratch_type* nrs=NULL)
{
    typedef detail::neighbourhood01_iter<A, G> nIter;
    return std::make_pair(
            nIter(b, e, size, g, include_base, nrs),
            nIter(e, e, 0, g, include_base) );
}

template<class A, class G>
inline
typename detail::neighbourhood01_iter<A,G>::scratch_type*
new_neighbourhood_range_scratch(A const, G const&, unsigned size=0)
{ untested();
    return new typename detail::neighbourhood01_iter<A,G>::scratch_type(size);
}

template<class A, class G>
std::pair<detail::neighbourhood01_iter<A, G>,
          detail::neighbourhood01_iter<A, G> >
inline make_neighbourhood_range(A b, A e, G const& g, unsigned size=0,
        base_t include_base=bINCLUDE,
        typename detail::neighbourhood01_iter<A,G>::scratch_type* nrs=NULL)
{
    return make_neighbourhood01_range(b, e, g, size, include_base, nrs);
}

#define first_adj(i) \
    boost::adjacent_vertices(i,_g).first
#define second_adj(i) \
    boost::adjacent_vertices(i,_g).second


namespace detail{
template<class A, class G, class V>
class neighbourhood_visitor { //
public: // types
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
    neighbourhood_visitor(const neighbourhood_visitor& p) :
        _i(p._i), _e(p._e),
        _g(p._g),
        _current(p._current),
        _visited(p._visited),
        _isatend(p._isatend)
    {
        // trace3("copy", &p, p._e-p._i, _e-_i);
    }
    neighbourhood_visitor(const neighbourhood_visitor&& p) :
        _i(p._i), _e(p._e),
        _g(p._g),
        _current(p._current),
        _visited(p._visited),
        _isatend(p._isatend)
    {
    }
    neighbourhood_visitor(A i, A const& e, G const& g, V const& v) :
        _i(i), _e(e),
        _g(g),
        _visited(v),
        _isatend(false)
    {
        // trace3("NV", this, _e-_i, (_e-_i)?boost::degree(*_i,_g):-1);
        if (_i==_e){
            _isatend=true;
            return;
        }else{
            assert(treedec::is_valid(*_i,_g));
            auto P=boost::adjacent_vertices(*_i,_g);
            unsigned t=0;
            unsigned v=0;
            for(;P.first!=P.second;++P.first){
                ++t;
                auto p=boost::get(boost::vertex_index, _g, *P.first);
                if(!_visited[p]){
                    ++v;
                }else{
                }
            }
        }
        _current = first_adj(*_i);

        // skip to next valid neighbour.
        while(_current==second_adj(*_i)){ untested();
            ++_i;
            assert(size_t(_e-_i)<=boost::num_vertices(_g));
            if(_i==_e){ untested();
                trace1("initial isatend", 1);
                _isatend=true;
                return;
            }else{ untested();
            }
            _current = first_adj(*_i);
        }
        assert(size_t(_e-_i)<=boost::num_vertices(_g));
        // then increment if visited.
        if(visited()){
            operator++();
        }else{
        }
        assert(size_t(_e-_i)<=boost::num_vertices(_g));
    }
    neighbourhood_visitor& operator=(const neighbourhood_visitor& x)
    {
        if(x._e==x._i){
            _i=x._i;
        }else{ untested();
            assert(0); //incomplete?;
        }
        return *this;
    }
public: // ops
    bool operator==(const neighbourhood_visitor& x) const
    {
        return !operator!=(x);
    }
    bool operator<(const neighbourhood_visitor& x) const
    { incomplete();
        if(_i < x._i){ untested();
            return true;
        }else if(_i==_e){untested();
            return false;
        }else if(x._i==x._e){untested();
            unreachable(); //?
            return true;
        }else{untested();
            // HACK HACK
            return &_current<&x._current;
        }
        return _isatend; // hack
    }
    bool operator!=(const neighbourhood_visitor& x) const
    {

//        x.update(); ?!

        if(x._i != _i){
        }else if(x._i==x._e){
            return false;
        }else{untested();
            // assert(!visited()); externally modified!
            // assert(!x.visited()); same!
            return x._current!=_current;
        }
        return !_isatend; // hack
    }
    neighbourhood_visitor& operator++()
    {
        assert(_i!=_e);
        assert(!_isatend);
        assert(size_t(_e-_i)<=boost::num_vertices(_g));
        assert(_current!=second_adj(*_i));
        // assert(visited()); // no, externally modified?
        assert(_current!=second_adj(*_i));
        assert(treedec::is_valid(*_current, _g));

        ++_current;
        if(_current==second_adj(*_i)) {
        }else if(!visited()) {
            // trace2("skip to", _e-_i, *_current);
            return *this;
        }else{
        }

        while(true){
            assert(_i!=_e);
            if(_current==second_adj(*_i)){
                ++_i;
                assert(size_t(_e-_i)<=boost::num_vertices(_g));
                if(_i==_e){
                    _isatend = true;
                    return *this;
                }else{
                    assert(treedec::is_valid(*_i, _g));
                    assert(_i!=_e);
                    trace2("hmm", *_i, boost::num_vertices(_g));
                    _current = first_adj(*_i);
                    if(_current == second_adj(*_i)){ untested();
                    }
                }
            }else{
                assert(treedec::is_valid(*_current, _g));
            }
            assert(_current!=second_adj(*_i));
            assert(treedec::is_valid(*_current, _g));
            if(!visited()) {
                return *this;
            }else{
            }
            ++_current;
        }
        untested();
    }
    const vertex_descriptor operator*()
    {
        assert(!_isatend);
        assert(_i!=_e);
        assert(size_t(_e-_i)<=boost::num_vertices(_g));
        return *_current;
    }
private: // impl
    bool visited() const
    {
        assert(treedec::is_valid(*_current, _g));
        auto p=boost::get(boost::vertex_index, _g, *_current);
        return _visited[p];
    }
private: // data
    A _i;
    A const& _e;
    G const& _g;
    adjacency_iterator _current;
    V const& _visited;
    bool _isatend;
};
} // detail

template<class A, class G, class V>
std::pair<detail::neighbourhood_visitor<A, G, V>,
          detail::neighbourhood_visitor<A, G, V> >
inline make_neighbourhood_range(A b, A const& e, G const& g, V const& visited)
{
    typedef detail::neighbourhood_visitor<A, G, V> nIter;
    return std::make_pair(
            nIter(b, e, g, visited),
            nIter(e, e, g, visited) );
}


namespace treedec{
// TODO: where?
template<class G>
void assert_connected(G const & g)
{
    (void)g;
#ifndef NDEBUG
#if __cplusplus >= 201103L
    std::vector<BOOL> mask(boost::num_vertices(g));
    mask.assign(boost::num_vertices(g), false);
    auto b=boost::vertices(g).first;
    auto e=b; ++e;
    auto BFS=make_bfs_range(b, e, g, &mask);

    unsigned nv=0;
    for(;BFS.first!=BFS.second;++(BFS.first)){
        ++nv;
    }

    trace2("ok?", nv, boost::num_vertices(g));
    assert(nv==boost::num_vertices(g));
    trace2("ok", nv, &g);
#endif
#endif
}

namespace draft{

// better use boost::range?

template<class iter1, class iter2>
class concat_iterator{
public:
    concat_iterator(iter1 begin1, iter1 end1, iter2 begin2, iter2 end2)
     : _i1(begin1), _e1(end1), _i2(begin2), _e2(end2){}

    bool operator==(const concat_iterator& o) const{
        return _i1==o._i1 && _i2==o._i2;
    }
    bool operator!=(const concat_iterator& o) const{
        return !operator==(o);
    }

    bool operator!=(const iter2& /*end*/) const{
        // warning: only works if end==end_of range2.

        return !(_i1 ==_e1 && _i2 == _e2);
    }

    void operator++(){
        if(_i1!=_e1){
            // still busy with range 1
            ++_i1;
        }
        else{
            ++_i2;
        }
    }

    unsigned operator*() const{
        if(_i1!=_e1){
            //still busy with range 1
            return *_i1;
        }

        return *_i2;
    }

    bool is_in_underlying(){
        return _i1!=_e1;
    }

private:
    iter1 _i1, _e1;
    iter2 _i2, _e2;
};

template<class i1, class i2>
concat_iterator<i1, i2> make_concat_iterator(i1 a, i1 b, i2 c, i2 d)
{ untested();
    return concat_iterator<i1, i2>(a, b, c, d);
}

} // draft

} // treedec

namespace std{
    // incomplete.
    template<class A, class B>
    struct iterator_traits<treedec::draft::concat_iterator<A, B> >{
        typedef typename A::value_type value_type;
        typedef typename A::difference_type difference_type;
        typedef typename A::reference reference;
        typedef typename std::forward_iterator_tag iterator_category;
    };
}

#endif

// vim:ts=8:sw=4:et
