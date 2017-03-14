// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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

#ifndef TD_IMPL_ELIM_ORD
#define TD_IMPL_ELIM_ORD

#ifndef TD_ELIMINATION_ORDERINGS
#error "not intended to be used like that."
#endif

#include "../algo.hpp"

namespace treedec{ //

namespace impl{ //

// is there a common denominator?
// class greedy_heuristic{
//    endless()
//    ordering()
//    decomposition()
// };

// the minDegree heuristic.

template <typename G_t, typename T_t, typename O_t>
class greedy_heuristic_base : public ::treedec::algo::draft::algo1{
public:
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    greedy_heuristic_base(G_t &G, T_t *T, O_t *O, unsigned ub, bool ignore_isolated_vertices=false)
      : algo1("."), _g(G), _t(T), _o(O), _own_o(!O), _ub_in(ub), _iiv(ignore_isolated_vertices), _i(0), _min(0), _current_N(&bag_i), _num_vert(boost::num_vertices(_g))
    {
        if(_own_o){
            _o = new O_t;
        }

        //the following seems to be unnecessary
        if(_t){
            _bags.resize(_num_vert);
        }
        _o->resize(_num_vert);
    }


    ~greedy_heuristic_base()
    {
        if(_own_o){
            delete _o;
        }
    }

    void tree_decomposition()
    {
        treedec::detail::skeleton_to_treedec(_g, *_t, _bags, *_o, _i);
    }

    vertices_size_type get_bagsize()
    {
        return _ub+1;
    }

    O_t& elimination_ordering()
    {
        return *_o;
    }

    virtual void initialize() = 0;

    virtual void next(vertex_descriptor &c) = 0;

    virtual void eliminate(vertex_descriptor v) = 0;

    virtual void postprocessing() = 0;

    void do_it()
    {
        timer_on();

        if(!_num_vert){
            timer_off();
            return;
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
            INTERRUPTION_POINT;

            // recompute ntd can only increase from here
            vertex_descriptor c;

            next(c);

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){
                _t->clear(); //could be also not the case
                throw exception_unsuccessful();
            }

            elim_vertices[_i] = get_vd(_g, c);

            if(_t){
                _current_N = &_bags[_i];
            }
            else if(_min > _ub){
                _ub = _min;
            }
            assert(bags_i);

            eliminate(c);

            if(!_t){
                _current_N->clear();
            }

            ++_i;
        }

        postprocessing();

        assert(_i == num_vert); //or _i-1==num_vert?!

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

    std::vector<bag_type> _bags; // BUG. use _t;

    vertices_size_type _ub;

    bag_type bag_i;
    bag_type* _current_N;

    unsigned _num_vert;

};

template <typename G_t, typename T_t, typename O_t>
class greedy_heuristic_base_vec : public ::treedec::algo::draft::algo1{
public:
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename std::vector<vertex_descriptor> bag_type;

    greedy_heuristic_base_vec(G_t &G, T_t *T, O_t *O, unsigned ub, bool ignore_isolated_vertices=false)
      : algo1("."), _g(G), _t(T), _o(O), _own_o(!O), _ub_in(ub), _iiv(ignore_isolated_vertices), _i(0), _min(0), _current_N(&bag_i), _num_vert(boost::num_vertices(_g))
    {
        if(_own_o){
            _o = new O_t;
        }

        //the following seems to be unnecessary
        if(_t){
            _bags.resize(_num_vert);
        }
        _o->resize(_num_vert);
    }


    ~greedy_heuristic_base_vec()
    {
        if(_own_o){
            delete _o;
        }
    }

    void tree_decomposition()
    {
        treedec::detail::skeleton_to_treedec_vec(_g, *_t, _bags, *_o, _i);
    }

    vertices_size_type get_bagsize()
    {
        return _ub+1;
    }

    O_t& elimination_ordering()
    {
        return *_o;
    }

    virtual void initialize() = 0;

    virtual void next(vertex_descriptor &c) = 0;

    virtual void eliminate(vertex_descriptor v) = 0;

    virtual void postprocessing() = 0;

    void do_it()
    {
        timer_on();

        if(!_num_vert){
            timer_off();
            return;
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
            INTERRUPTION_POINT;

            // recompute ntd can only increase from here
            vertex_descriptor c;

            next(c);

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){
                _t->clear(); //could be also not the case
                throw exception_unsuccessful();
            }

            elim_vertices[_i] = get_vd(_g, c);

            if(_t){
                _current_N = &_bags[_i];
            }
            else if(_min > _ub){
                _ub = _min;
            }
            assert(bags_i);

            eliminate(c);

            if(!_t){
                _current_N->clear();
            }

            ++_i;
        }

        postprocessing();

        assert(_i == num_vert); //or _i-1==num_vert?!

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

    std::vector<bag_type> _bags; // BUG. use _t;

    vertices_size_type _ub;

    bag_type bag_i;
    bag_type* _current_N;

    unsigned _num_vert;

};


template <typename G_t, typename T_t, typename O_t>
class minDegree : public greedy_heuristic_base_vec<G_t, T_t, O_t>{
public:
    typedef greedy_heuristic_base_vec<G_t, T_t, O_t> baseclass;

    typedef typename deg_chooser<G_t>::type degs_type;

    minDegree(G_t &g, T_t *t, O_t *o,
                    unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, t, o, ub, ignore_isolated_vertices),
         _degs(baseclass::_g)
    {
    }

    minDegree(G_t &G, O_t& o, bool ignore_isolated_vertices)
        : baseclass(G, NULL, &o, -1u, ignore_isolated_vertices),
          _degs(baseclass::_g)
    {
    }

    void initialize(){
        auto zerodegbag1=MOVE(_degs.detach_bag(0));
        BOOST_AUTO(it, zerodegbag1.begin());

        if(!baseclass::_iiv){
            for(; it!=zerodegbag1.end(); ++it){
                (*baseclass::_o)[baseclass::_i++] = get_vd(baseclass::_g, *it);
            }
        }
        else{
            baseclass::_num_vert -= zerodegbag1.size();
            //we assume that the vertices in zerodegbag are the prefix of o, so we dont overwrite them
            //idea behind this: just use one orderding for all algos
            //baseclass::_i+=zerodegbag1.size(); //TODO: if needed...
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
            assert(*I!=c); // no self loops...
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

        for(; it!=zerodegbag.end(); ++it){ untested();
            (*baseclass::_o)[baseclass::_i++] = get_vd(baseclass::_g, *it);
        }
    }

private:
    degs_type _degs;

}; // minDegree



// the fillIn heuristic.
template <typename G_t, typename T_t, typename O_t>
class fillIn : public greedy_heuristic_base<G_t, T_t, O_t>{
public:
    typedef greedy_heuristic_base<G_t, T_t, O_t> baseclass;

    typedef typename fill_chooser<G_t>::type fill_type;


    struct fill_update_cb : public graph_callback<G_t>{ //
        typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;

        fill_update_cb(fill_type* d, G_t const& g) :
            _fill(d), G(g){}

        void operator()(vertex_descriptor v)
        {
            _fill->q_eval(v);
        }
        void operator()(vertex_descriptor s, vertex_descriptor t)
        {
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


    fillIn(G_t &g, T_t *t, O_t *o,
                    unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, t, o, ub, ignore_isolated_vertices),
         _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    {
    }

    fillIn(G_t &G, O_t& o, bool ignore_isolated_vertices, unsigned ub=-1u)
        : baseclass(G, NULL, &o, ub, ignore_isolated_vertices),
          _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    {
    }



public: // implementation

    void initialize(){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(baseclass::_g); vIt != vEnd; ++vIt){
            if(boost::out_degree(*vIt, baseclass::_g) == 0){
                if(!baseclass::_iiv){
                    (*baseclass::_o)[baseclass::_i++] = get_vd(baseclass::_g, *vIt);
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

} // namespace impl

} // namespace treedec

#endif // guard

// vim:ts=8:sw=4:et
