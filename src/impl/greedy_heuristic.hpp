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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
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

namespace treedec{

namespace impl{

// TODO: how does this relate to greedy_base?
template<class G_t, template<class G, class...> class CFGT_t=algo::default_config>
class greedy_heuristic_base : public ::treedec::algo::draft::algo1{
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename treedec::graph_traits<G_t>::treedec_type T_t;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename std::vector<vertex_descriptor> bag_t;
    typedef typename std::vector<vertex_descriptor> O_t;

    greedy_heuristic_base(G_t &G, unsigned ub, bool ignore_isolated_vertices=false)
      : algo1("."), _g(G), _t(NULL), _own_o(),
        _ub_in(ub), _iiv(ignore_isolated_vertices), _i(0),
        _min(0), _ub(0), _current_N(&bag_i), _num_vert(boost::num_vertices(_g))
    {
        _o = new O_t;
        _o->resize(_num_vert);
        _do_tree_decomposition=true; // for now.
    }


    ~greedy_heuristic_base(){
        if(_own_o){
            delete _o;
        }
    }

    O_t& get_elimination_ordering() {

        return *_o;
    }

#if 0
    template<class O>
    void get_elimination_ordering(O& o) const { untested();
        O = *_o; // doesnt work like this.
    }
#endif

    template<class T>
    void get_tree_decomposition(T& t)const{
        // std::cerr << "hmm " << _o->size() << " " << _bags.size() << "\n";
        assert(_o->size()<=_bags.size()); // this is obsolete anyway
        assert(_o->size()==_num_vert); // this is relevant.

        // yuck... will be obsolete with FI rework
        typename std::vector<
            std::pair<vertex_descriptor, bag_t>
                > bags(_num_vert);
        typename std::vector<unsigned> io(_num_vert);

        // stuff center and friends into "skeleton"
        // _num_vert can be less than order/bags size
        for(unsigned i = 0; i < _num_vert; i++){
            bags[i].first = (*_o)[i];
            bags[i].second = _bags[i];
            // io[ (*_o)[i] ] = i;
        }

        treedec::detail::skeleton_to_treedec(_g, t, bags, *_o, _i);
    }

    vertices_size_type get_bagsize(){
        return _ub+1;
    }

    virtual void initialize() = 0;
    virtual void next(vertex_descriptor &c) = 0;
    virtual void eliminate(vertex_descriptor v) = 0;
    virtual void postprocessing() = 0;

    void disable_td(){
        // stupid hack
        _do_tree_decomposition = true;
    }

    void do_it() {
        if(_do_tree_decomposition){
            _t=new T_t;
            // bags seem to be unnecessary
            _bags.resize(_num_vert);
        }else{untested();
        }


        timer_on();

        if(!_num_vert){
            timer_off();
            return;
        }else{
        }

        assert(_o);
        O_t& elim_vertices = *_o;

#ifndef NDEBUG
        check(_g);
#endif

        initialize();

        _o->resize(_num_vert);

        assert(elim_vertices.size() == _num_vert);

        while(boost::num_edges(_g) > 0){
            vertex_descriptor c;

            next(c);

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){ untested();
                assert(_t); // ouch?
                _t->clear(); //could be also not the case
                throw exception_unsuccessful();
            }

            elim_vertices[_i] = c;

            if(_t){
                _current_N = &_bags[_i];
            }else{
            }

            _ub = (boost::out_degree(c, _g)>_ub)?boost::out_degree(c, _g):_ub;

            // assert(bags_i);?!

            eliminate(c);

            if(!_t){
                _current_N->clear();
            }else{
            }

            ++_i;
        }

        postprocessing();

//        assert(_i == num_vert); //or _i-1==num_vert?!

        timer_off();

    } // do_it


protected:
    G_t &_g;
    T_t* _t;
    O_t* _o;
    bool _own_o;

    vertices_size_type _ub_in;
    bool _iiv;
    size_t _i;
    unsigned _min;

    std::vector<bag_t> _bags;

    vertices_size_type _ub;

    bag_t bag_i;
    bag_t* _current_N;

    unsigned _num_vert;
private:
    bool _do_tree_decomposition;
};


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
               CFGT>{
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
        {
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
    using baseclass::_num_edges;
    using baseclass::_degree; // bug?

    // fillIn::
    bool next(typename baseclass::vertex_descriptor &c){
        trace1("next", _num_edges);
        if(_num_edges){
            // todo: what do we know about lower bound?
            auto p=_fill.pick_min(0, -1u, true);
            c = p.first;
            _min = p.second; // the fill of c.
            return true;
        }else{
            return false;
        }
    }

    // TODO more useful specialisation?
    // ... finish minDegree, then lets see.
    void eliminate(typename baseclass::vertex_descriptor c){
        long fill_c=_min;
        trace4("elim", c, _num_edges, _degree[c], fill_c);

        /// _min is fill(v).
        // use remove_out_edge_if?
        _fill.mark_neighbours(c, _min); // relevant in q_decrement
                                       // different marker
                                       //
        // bug: idmap!
        auto degc=baseclass::_degreemap[c];

        assert(fill_c<=degc*(degc-1));

        baseclass::_numbering.put(c);
        baseclass::_numbering.increment();

        { // make clique
            // todo?: faster special cases for small numbers?!
            // if(fill==1 && degree==2){
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
                size_t overlap = mark_neighbours_c(
                        baseclass::_marker, n, baseclass::_subgraph,
                        _fill.marked() );
                long degn=baseclass::_degreemap[n];
                auto const fill_n=_fill.get_value(n);
                trace5("------------> ", n, degn, fill_n, degc, overlap);
                assert(overlap<size_t(degn));
                assert(overlap<size_t(degc));

                long DC = degc - overlap - 1; // nodes connected to c but not (to) n
                long DN = degn - overlap - 1; // nodes connected to n but not (to) c

                trace5("DC?", fill_n, fill_c, degc, overlap, degn);
                trace2("DC?", DC, DN);
                long offset = - long(fill_c);
#if 0
              if(degc==2 && overlap==1){
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
                    trace3("DC", n, fill_n, offset);

                    _fill.shift(n, offset);
                    _fill.q_eval(n); // treat as lb
                }else{ // !DC
                    // this is exact, unless an edge is missing.
                    // (then it is flagged lb, later);
                    offset = - long(fill_c) - DN;
                    trace3("noDC", n, fill_n, offset);
                    _fill.shift(n, offset);
                }

                auto q=adjacent_vertices(c, _subgraph);
                ///q.first=next;

                // iterate {n2,n} \subset 1-neighborhood
                // n2 < n... add edges to previously visited n2 only
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
        auto p=boost::adjacent_vertices(c, _subgraph);
        for(; p.first!=p.second; ++p.first){
            auto n=*p.first;
            long degn=baseclass::_degreemap[n];
            trace3("check", n, _fill.get_value(n), degn);
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
#endif
    } // eliminate(c)

    using baseclass::_i;
    using baseclass::_subgraph;
    void postprocessing(){
        if(_i == baseclass::_num_vert){
            // no nodes at all?!
        }else{
            // the last node is missing, but why?
            ++_i;
            auto v = _fill.pick_min(0, 0, true).first;
            assert(_i == baseclass::_o->size());
            baseclass::_o->back() = v;
            baseclass::_numbering.put(v);
        }
        assert(baseclass::_i == baseclass::_num_vert);
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

    void next(typename baseclass::vertex_descriptor &c){
        _fill.check();
        boost::tie(c, baseclass::_min) = _fill.pick_min(0, -1, true);
        _fill.check();
    }

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
