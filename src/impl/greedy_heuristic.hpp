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
//
// greedy heuristics

#ifndef TREEDEC_GREEDY_HEURISTIC_HPP
#define TREEDEC_GREEDY_HEURISTIC_HPP

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#error "not intended to be used like that."
#endif

#include "../algo.hpp"
#include "greedy_base.hpp"
#include "obsolete_greedy_base.hpp"

namespace treedec{

namespace impl{

template <typename G_t, template<class G, class...> class CFG=algo::default_config>
class minDegree : public greedy_heuristic_base<G_t, CFG>{
public:
    typedef greedy_heuristic_base<G_t, CFG> baseclass;

    typedef typename deg_chooser<G_t>::type degs_type;

    minDegree(G_t &g,
                    unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
         _degs(baseclass::_g)
    {
    }

    minDegree(G_t &G, bool ignore_isolated_vertices)
        : baseclass(G, -1u, ignore_isolated_vertices),
          _degs(baseclass::_g)
    {
    }

#if 0 // base
    void get_elimination_ordering(){ untested();
        // incomplete()
    }
#endif

    void initialize(){
        auto zerodegbag1=MOVE(_degs.detach_bag(0));
        BOOST_AUTO(it, zerodegbag1.begin());

        if(!baseclass::_iiv){
            for(; it!=zerodegbag1.end(); ++it){
                (*baseclass::_o)[baseclass::_i++] = *it;
            }
        }else{
            baseclass::_num_vert -= zerodegbag1.size();
        }

        baseclass::_min = 1;
    }

    void next(typename baseclass::vertex_descriptor &c){
        if(baseclass::_min>1){
            --baseclass::_min;
        }

        boost::tie(c, baseclass::_min) = _degs.pick_min(baseclass::_min, baseclass::_num_vert);
    }

    // md::
    void eliminate(typename baseclass::vertex_descriptor v){
        typename baseclass::adjacency_iterator I, E;
        for(boost::tie(I, E) = boost::adjacent_vertices(v, baseclass::_g); I!=E; ++I){
            assert(*I!=v); // no self loops...
            typename baseclass::vertex_descriptor w=*I;
            _degs.unlink(w);
        }

        baseclass::_current_N->resize(boost::out_degree(v, baseclass::_g));

        make_clique_and_detach(v, baseclass::_g, *baseclass::_current_N);

        redegree(NULL, baseclass::_g, *baseclass::_current_N, _degs);
        _degs.unlink(v, baseclass::_min);
        _degs.flush();
    }

    void postprocessing(){
        auto zerodegbag=MOVE(_degs.detach_bag(0));
        BOOST_AUTO(it, zerodegbag.begin());

        for(; it!=zerodegbag.end(); ++it){
            (*baseclass::_o)[baseclass::_i++] = *it;
        }
    }

private:
    degs_type _degs;

}; // minDegree

// the fillIn heuristic.
template<typename G_t,
         template<class GG, class ...> class CFGT=algo::default_config>
class fillIn : public greedy_base< G_t,
               std::vector<typename boost::graph_traits<G_t>::vertex_descriptor>,
               CFGT>{ //
public: //types
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_t;
    typedef typename directed_view_select<G_t>::type D_t;
    typedef typename boost::graph_traits<D_t>::vertices_size_type vertices_size_type;
    typedef greedy_base<G_t, O_t, CFGT> baseclass;
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
          _fill(baseclass::_subgraph, boost::num_vertices(g))
    { untested();
//        boost::print_graph(g);
        treedec::check(g);
    }
    fillIn(G_t &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
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
            trace2("next picked", c, _min);

            if(_min<0){
                // 
            }else{
            }

            _min = p.second; // the fill of c.
            return true;
        }
    }

    // TODO more useful specialisation?
    // ... finish minDegree, then lets see.
    // fillIn::
    void eliminate(typename baseclass::vertex_descriptor c){
        long fill_c = _min;
        trace2("elim", vertices_left(), _degree[c]);
        trace2("elim", _fill.max_fill(), _fill.is_lb(c) );

#ifndef NDEBUG
        // check_vertex(c);
#endif

        /// _min is fill(v).
        // use remove_out_edge_if?
        _fill.mark_neighbours(c, _min); // relevant in q_decrement
                                       // different marker
                                       //
        // bug: idmap!
        auto degc = baseclass::_degreemap[c];

        assert(fill_c<=long(degc*(degc-1)));

        baseclass::_numbering.put(c);
        baseclass::_numbering.increment();

        { // make clique
            // todo?: faster special cases for small numbers?!
            // if(fill==1 && degree==2){ untested();
            //    easy.
            // }
            //
            assert(baseclass::_num_edges >= degc);
            baseclass::_num_edges -= degc;
            trace3("c...", c, fill_c, degc);
            auto p=boost::adjacent_vertices(c, baseclass::_subgraph);
            for(; p.first!=p.second; ++p.first){
                auto n=*p.first;

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
    auto& g = baseclass::_subgraph;
    auto& p = _fill.marked();
    auto pp=boost::adjacent_vertices(n, g); // ???
    for(; pp.first!=pp.second; ++pp.first){
        baseclass::_marker.mark(*pp.first);
        if(p(*pp.first)){
            ++overlap;
        }else{
        }
    }
}
#endif

                long degn=baseclass::_degreemap[n];
                auto const fill_n=_fill.get_value(n);
                trace5("------------> ", n, degn, fill_n, degc, overlap);
                assert(overlap<size_t(degn));
                assert(overlap<size_t(degc));

                long DC = degc - overlap - 1; // nodes connected to c but not (to) n
                long DN = degn - overlap - 1; // nodes connected to n but not (to) c

                trace5("DC?", fill_n, fill_c, degc, overlap, degn);
                trace5("DC?", c, n, DC, DN, _fill.is_lb(n));

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
                    offset -= DN; // no more need to connect edges behind
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

                ///q.first=next;

                // iterate {n2,n} \subset 1-neighborhood
                // n2 < n... add edges to previously visited n2 only
                auto q = adjacent_vertices(c, _subgraph);
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
                                trace3("--common neigh", n, n2, *r.first);
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
            } // n
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

            for(; _i < baseclass::_num_vert; ++_i){ untested();
                auto v = _fill.pick_min(0, 0, true).first;
//                treedec::add_edge(w, v, baseclass::_g); already there.
                assert(_i < _o->size());
                (*_o)[baseclass::_i] = v;
            }
        }
        assert(baseclass::_i == baseclass::_num_vert);
    }

private: // debugging
    void check_vertex(vertex_descriptor c) {
        size_t k=0;
        for(auto p=boost::adjacent_vertices(c, _subgraph); p.first!=p.second; ++p.first){
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
                auto n2=*q.first;
                trace1("neigh", n2);
                --degn;
            }
            trace2("check", n, degn);
            assert(!degn);
        }
    }
private:
    fill_type _fill;
//    fill_update_cb _cb;
}; // fillIn

} // impl

namespace obsolete {
// the fillIn heuristic.
template <typename G_t, template<class G, class...> class CFGT_t=algo::default_config>
class fillIn : public treedec::impl::greedy_heuristic_base<G_t, CFGT_t>{
public: //types
    typedef treedec::impl::greedy_heuristic_base<G_t, CFGT_t> baseclass;
    typedef typename treedec::obsolete::FILL<G_t> fill_type;

    struct fill_update_cb : public graph_callback<G_t>{
        typedef typename baseclass::vertex_descriptor vertex_descriptor;

        fill_update_cb(fill_type* d, G_t const& g) :
            _fill(d), G(g){}

        void operator()(vertex_descriptor v){
            _fill->q_eval(v);
        }
        void operator()(vertex_descriptor s, vertex_descriptor t) {
            assert(s < t); // likely not. is this necessary below?
            // e has just been inserted.
            BOOST_AUTO(cni, common_out_edges(s, t, G));
            BOOST_AUTO(i, cni.first);
            BOOST_AUTO(e, cni.second);
            for(; i!=e; ++i){
                assert(*i != s);
                assert(*i != t);
    //            no. maybe theres only half an edge.
    //            assert(boost::edge(boost::source(edg, G), *i, G).second);
    //            assert(boost::edge(boost::target(edg, G), *i, G).second);

                // BUG: *i might be within 1-neighborhood.
                _fill->q_decrement(*i);
            }
        }
    private:
        fill_type* _fill;
        G_t const& G;
    }; // update_cb

public: // construct
    fillIn(G_t &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
         _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    {
    }

    fillIn(G_t &G, bool ignore_isolated_vertices, unsigned ub=-1u)
        : baseclass(G, ub, ignore_isolated_vertices),
          _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    {
    }

public: // implementation
    void initialize(){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(baseclass::_g); vIt != vEnd; ++vIt){
            if(boost::out_degree(*vIt, baseclass::_g) == 0){
                if(!baseclass::_iiv){
                    (*baseclass::_o)[baseclass::_i++] = *vIt;
                }
                else{
                    --baseclass::_num_vert;
                }
            }
        }
    }

    // obs::fillIn::
    void next(typename baseclass::vertex_descriptor &c){
        _fill.check();
        boost::tie(c, baseclass::_min) = _fill.pick_min(0, -1, true);
        _fill.check();
    }

    // obs::fillIn::
    void eliminate(typename baseclass::vertex_descriptor v){
        _fill.mark_neighbors(v, baseclass::_min);

        baseclass::_current_N->resize(boost::out_degree(v, baseclass::_g));

        make_clique_and_detach(v, baseclass::_g, *baseclass::_current_N, &_cb);

        _fill.unmark_neighbours(*baseclass::_current_N);
    }

    void postprocessing(){
        for(; baseclass::_i < baseclass::_num_vert; ++baseclass::_i){
            auto v = _fill.pick_min(0, 0, true).first;
            (*baseclass::_o)[baseclass::_i] = v;
        }
    }

private:
    fill_type _fill;
    fill_update_cb _cb;

}; // fillIn

} // obsolete

} // treedec

#endif // guard

// vim:ts=8:sw=4:et
