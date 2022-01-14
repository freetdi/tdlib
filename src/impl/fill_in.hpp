// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
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
//
// fill in heuristic

#ifndef TREEDEC_FILL_IN_HPP
#define TREEDEC_FILL_IN_HPP

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#error "not intended to be used like that."
#endif

#include "../algo.hpp"
#include "greedy_base.hpp"
#include "obsolete_greedy_base.hpp"

namespace treedec{
namespace impl{
namespace detail{

template<class G, class O, class F>
class eliminated_before{
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map< G, boost::vertex_index_t >::const_type::value_type vertex_index_type;
public:
    explicit eliminated_before(vertex_descriptor c, O const& num, G const& g, F& f)
        : _c(c), _numbering(num), _g(g), _fill(f) {
    }
    template<class E>
    bool operator()(E const& e) const{
        auto t = boost::target(e, _g);
        trace2("b4", t, _c);
        if( _numbering.is_before(t, _c)){
            return true;
        }else{
            _fill.mark(t);
            return false;
        }
    }
private:
    vertex_descriptor _c;
    O const& _numbering;
    G const& _g;
    F& _fill;
};

}

// the fillIn heuristic.
template<typename G_t,
         template<class GG, class ...> class CFGT=algo::default_config>
class fillIn : public greedy_base< G_t,
               std::vector<typename boost::graph_traits<G_t>::vertex_descriptor>,
               CFGT>{ //
public: //types
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O;
    typedef O O_t; //?
    typedef typename directed_view_select<G_t>::type D_t;
    typedef typename boost::graph_traits<D_t>::vertices_size_type vertices_size_type;
#ifdef DEBUG
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
    marker_type _debug_marker;
#endif
    typedef greedy_base<G_t, O_t, CFGT> baseclass;
    typedef typename baseclass::numbering_type numbering_type;
    typedef typename baseclass::vertex_descriptor vertex_descriptor;
    // BUG:: use CFGT::fill or fallback to current fill
    typedef typename fill_chooser<typename baseclass::subgraph_type>::type fill_type;

    struct fill_update_cb : public graph_callback<typename baseclass::subgraph_type>{
        typedef typename baseclass::subgraph_type G;

        fill_update_cb(fill_type* d, G const& g) :
            _fill(d), _g(g)
        { untested();
        }

        void operator()(vertex_descriptor v){ untested();
            unreachable();
            _fill->q_eval(v);
        }
        // q_decrement nodes that are incident to both endpoints.
        void operator()(vertex_descriptor, vertex_descriptor) { unreachable();
        }
    private:
        fill_type* _fill;
        G const& _g;
    }; // update_cb

public: // construct
    fillIn(G_t const &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
#ifdef DEBUG
        later.
        _debug_marker(boost::num_vertices(g)),
#endif
          _fill(baseclass::_subgraph, boost::num_vertices(g))
    { untested();
//        boost::print_graph(g);
        treedec::check(g);
    }
    fillIn(G_t &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
#ifdef DEBUG
        _debug_marker(boost::num_vertices(g)),
#endif
          _fill(baseclass::_subgraph, boost::num_vertices(g))
    {
//        boost::print_graph(g);
        treedec::check(g);
    }

    fillIn(G_t &g, bool ignore_isolated_vertices, unsigned ub=-1u)
        : baseclass(g, ub, ignore_isolated_vertices),
          _fill(baseclass::_subgraph, boost::num_vertices(g))
          // _cb(fill_update_cb(&_fill, baseclass::_subgraph))
    { untested();
//        boost::print_graph(g);
    }
#ifdef DEBUG_FILL
private: // debugging
    size_t fill_cached_(vertex_descriptor v) const override {
        return _fill.get_value(v);
    }
    bool fill_cached_is_lb_(vertex_descriptor v) const override {
        return _fill.is_lb(v);
    }
#endif

public: // implementation
    using baseclass::_min;
    using baseclass::vertices_left;
    using baseclass::_num_edges;
    using baseclass::_degree; // bug?

    // fillIn::
    bool next(typename baseclass::vertex_descriptor &c){
        auto n = vertices_left();
        trace2("next loop", _num_edges, vertices_left());

        if (0 && _num_edges == (n*(n-1))/2){ untested();
            // something is wrong with the condition
            return false;
        }else if(!_num_edges){
            return false;
        }else{
            // todo: what do we know about lower bound?
            auto p = _fill.pick_min(0, -1u, true);
            c = p.first;
            trace2("next picked", c, p.second);

            if(_min<0){
                // 
            }else{
            }

            _min = p.second; // the fill of c.
            return true;
        }
    }

    using baseclass::_g;
    using baseclass::_numbering;
//    using baseclass::_fill;
    // TODO more useful specialisation?
    // ... finish minDegree, then lets see.
    // fillIn::
    void eliminate(typename baseclass::vertex_descriptor c){

        assert( _fill.get_value(c) == _min);

        long fill_c = _min;
        trace2("elim", vertices_left(), _degree[c]);
        trace2("elim", _fill.max_fill(), _fill.is_lb(c) );
        trace2("degree", boost::out_degree(c, baseclass::_g), _degree[c]);

#ifdef DEBUG
        check_vertex(c);
        auto me = treedec::count_missing_edges(c, _debug_marker, baseclass::_g);
        trace3("DEBUG", me, _min, c);
        assert(me==_min);
#endif
        // assert(_min<=baseclass::_num_vert); // no!

        /// _min is fill(v).
        // use remove_out_edge_if?
        //_fill.mark_neighbours(c, _min); // relevant in q_decrement ?
                                          // different marker
        _fill.clear_marker();
        detail::eliminated_before<D_t, numbering_type, fill_type> P(c, _numbering, _g, _fill);
        boost::remove_out_edge_if(c, P, _g);

        assert( boost::out_degree(c, baseclass::_g) ==  _degree[c]);

        // bug: idmap!
        auto degc = baseclass::_degreemap[c];

        assert(fill_c<=long(degc*(degc-1)));

        _numbering.put(c);
        _numbering.increment(); // remove c from _subgraph

        { // make clique
            // todo?: faster special cases for small numbers?!
            // if(fill==1 && degree==2){ untested();
            //    easy.
            // }
            //
            assert(baseclass::_num_edges >= degc);
            baseclass::_num_edges -= degc;
            trace3("c...", c, fill_c, degc);
            auto p = boost::adjacent_vertices(c, _g);
            for(; p.first!=p.second; ++p.first){
                auto n = *p.first;

                baseclass::_marker.clear();
                // count overlap with c neighbourhood
                // use remove_out_edge_if?
#if 0
                size_t overlap = mark_neighbours_c(
                        baseclass::_marker, n, baseclass::_subgraph,
                        _fill.marked() );
#else
                // template<class M, typename V, class G, class P>
                // size_t mark_neighbours_c(M& marker, V v, G const& g, P const& p /*bug*/)
                size_t overlap=0;
                {
                    auto& p = _fill.marked(); // neighbours of c
                    auto pp = boost::adjacent_vertices(n, _g); // ???
                    for(; pp.first!=pp.second; ++pp.first){
                        auto nn = *pp.first;
                        baseclass::_marker.mark(nn);
                        // assert(nn!=c); // was: _subgraph.
                        if(p(nn)){
                            ++overlap;
                        }else{
                        }
                    }
                }
#endif

                long degn = baseclass::_degreemap[n];
                auto const fill_n=_fill.get_value(n);
                trace5("------------> ", n, degn, fill_n, degc, overlap);
                assert(overlap<size_t(degn));
                assert(overlap<size_t(degc));

                long DC = degc - overlap - 1; // nodes connected to c but not (to) n
                long DN = degn - overlap - 1; // nodes connected to n but not (to) c

                trace5("DC?", fill_n, fill_c, degc, overlap, degn);
                trace5("DC?", c, n, DC, DN, _fill.is_lb(n));
                assert(DC>=0);
                assert(DN>=0);

                long offset = - long(fill_c);
#if 0
              if(degc==2 && overlap==1){ untested();
                  no need to queue.
                  assert(fillc==0)
                  offset = - long(DN);
                }else
#endif
                    
                if(DC){
                    offset = - long(fill_c);
//                    offset -= overlap*DC + (DC-1)*DC/2
                    offset -= DN; // no more need to connect edges behind c
//                    offset = - long(fill_n); // TODO
                    trace4("DC", n, fill_n, offset, _fill.is_lb(n));

                    _fill.shift(n, offset);
                    _fill.q_eval(n); // treat as lb
                }else{
                    // this is exact, unless an edge is missing.
                    // (then it is flagged lb, later);
                    offset = - long(fill_c) - DN;
                    trace6("noDC, shift", n, fill_c, fill_n, offset, DN, DC);
                    _fill.shift(n, offset);
                }
                trace3("--- done -->", n, _fill.get_value(n), _fill.is_lb(n));

                ///q.first=next;

                // iterate {n2,n} \subset 1-neighborhood
                // n2 < n... add edges to previously visited n2 only
                auto q = boost::adjacent_vertices(c, _g);
                for(; q.first!=p.first; ++q.first){
//                for(; q.first!=q.second; ++q.first) // careful!
                    auto n2=*q.first;
                    if(baseclass::_marker.is_marked(n2)){
                        // done, neighbour of n
                    }else{
                        if(degc==2){
                            // --- n2 --- c ---- n ---
                            // no need to reevaluate them
                            // neighbourhood does not change
                        }else{
                            // fill could decrease here.
                            _fill.q_eval(n);
                            _fill.q_eval(n2);
                        }
                        // assert(n2>n); // visit only once
                        assert(n2 != n);

                        trace2("addedge decfill", n, n2);
                        // mark 2-neighbors that are incident to both endpoints.
                        auto r=adjacent_vertices(n2, baseclass::_subgraph);
                        for(; r.first!=r.second; ++r.first){
                            if(!baseclass::_marker.is_marked(*r.first)){
                                // not neighbour of n
                            }else{
                                // trace3("--common neigh", n, n2, *r.first);
                                assert(*r.first!=n);
                                assert(*r.first!=n2);
                                _fill.decrement_fill(*r.first);
                            }
                        }

                        assert(!boost::edge(n, n2, baseclass::_g).second);
                        assert(!boost::edge(n2, n, baseclass::_g).second);
                        // why does it default to boost::??!
                        treedec::add_edge(n, n2, baseclass::_g);
                        assert(boost::edge(n, n2, baseclass::_g).second);
                        assert(boost::edge(n2, n, baseclass::_g).second);
                        // use _subgraph?!
                        ++baseclass::_degreemap[n2];
                        ++baseclass::_degreemap[n];
                        ++baseclass::_num_edges;
                        trace4("addedge", n, n2,
                                baseclass::_degreemap[n], baseclass::_degreemap[n2]);

                    }
                } // n2
                // disconnect center.
                --baseclass::_degreemap[n];
                trace2("incomplete n (missing edges)", n, baseclass::_degreemap[n]);
                degn = baseclass::_degreemap[n];
            } // neighbor loop
        }
        trace2("elimd", c, baseclass::_num_edges);
        treedec::check(_subgraph);

#ifndef NDEBUG
        check_vertex(c);
#endif
    } // eliminate(c)

    using baseclass::_i;
    using baseclass::_o; // BUG
    using baseclass::_subgraph;
    void postprocessing(){
        trace2("post", _i, baseclass::_num_vert);
        if(!baseclass::_iiv){
            assert(boost::num_vertices(_g) == baseclass::_num_vert);
        }else{
        }
        if(_i == baseclass::_num_vert){ untested();
            unreachable(); //?
            // no nodes at all?!
        }else{
            // the last node is missing, but why?
      //       auto v = _fill.pick_min(0, 0, true).first;
      //       (*_o)[_i++] = v;
      //       baseclass::_numbering.put(v);

      //       auto x = _subgraph.adjacent_vertices(v);
      //       for(;x.first!=x.second;++x.first){ untested();
      //           baseclass::_o->push_back(*x.first);
      //           (*_o)[_i++] = *x.first;
      //       }

            auto w = _fill.pick_min(0, 0, true).first;
            trace2("post", w, baseclass::_degreemap[w]);
            (*_o)[baseclass::_i++] = w;
            baseclass::_numbering.put(w);
            baseclass::_numbering.increment();
            auto x = _subgraph.adjacent_vertices(w);
            for(;x.first!=x.second;++x.first){
                trace1("last node neigh", *x.first);
            }

            for(; _i < baseclass::_num_vert; ++_i){
                auto v = _fill.pick_min(0, 0, true).first;
//                treedec::add_edge(w, v, baseclass::_g); already there.
                assert(_i < _o->size());
                (*_o)[baseclass::_i] = v;
            }
        }

        if(!baseclass::_iiv){
            assert(baseclass::_i == _o->size());
            assert(baseclass::_i == baseclass::_num_vert);
            assert(baseclass::_i == boost::num_vertices(_g));

            // assert(baseclass::_i == _numbering.total()); for some reason, numbering numbers the bags only.
        }else{
        }
    }

private: // debugging
    void check_vertex(vertex_descriptor c) {
        size_t k=0;
        for(auto p = boost::adjacent_vertices(c, _subgraph); p.first!=p.second; ++p.first){
            ++k;
        }
        assert(k==baseclass::_degreemap[c]);

        auto p=boost::adjacent_vertices(c, _subgraph);
        for(; p.first!=p.second; ++p.first){
            auto n=*p.first;
            long degn = baseclass::_degreemap[n];
            trace4("adj", c, *p.first, degn, baseclass::_degreemap[c]);
            trace4("checking neigh", n, _fill.get_value(n), degn, _fill.is_lb(n));
            assert(2*_fill.get_value(n)<=size_t(degn*(degn-1)));

            auto q=boost::adjacent_vertices(n, _subgraph);
            for(; q.first!=q.second; ++q.first){
//                auto n2=*q.first;
//                trace1("neigh", n2);
                --degn;
            }
            assert(!degn);
        }
        trace1("checked", c);
    }
private:
    fill_type _fill;
//    fill_update_cb _cb;
}; // fillIn

} // impl
} // treedec

#endif // guard

// vim:ts=8:sw=4:et
