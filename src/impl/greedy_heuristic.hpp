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
class minDegree{//
public:
//    typedef typename std::vector<typename boost::graph_traits<G>::vertex_descriptor> O_t;
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename deg_chooser<G_t>::type degs_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    minDegree(G_t &g, T_t *t, O_t *o,
                    unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : _g(g), _t(t), _o(o), _own_o(!o), _ub_in(ub), _iiv(ignore_isolated_vertices),
         _degs(_g)
    {
        _i=0;

        vertices_size_type num_vert=boost::num_vertices(_g);
        if(_own_o){
            _o = new O_t;
        }else{ untested();
        }

        if(_t){
            _bags.resize(num_vert);
        }
        _o->resize(num_vert);
    }

    minDegree(G_t &G, O_t& o, bool ignore_isolated_vertices):
        _g(G), _t(NULL), _o(&o), _own_o(false), _ub_in(-1u), _degs(_g)
    {
        if(!boost::num_vertices(_g)){
            incomplete();
        }

        _i = 0;

        unsigned int n_ = 0;

        //Mark isolated vertices as already visited.
        // DUPLICATE in fillIn
        if(ignore_isolated_vertices){ untested();
            _visited.resize(boost::num_vertices(G));
            typename boost::graph_traits<G_t>::vertex_iterator vit, vend;
            for(boost::tie(vit, vend) = boost::vertices(_g); vit != vend; vit++){
                if(boost::degree(*vit, _g) == 0){
                    unsigned int pos = get_pos(*vit, _g);
                    _visited[pos] = true;
                }else{
                    n_++;
                }
            }
            _o->resize(n_);
            trace3("fillIn with ignore", n_, boost::num_vertices(G), _ub_in);
        }else{
            _o->resize(boost::num_vertices(G));
            auto r=_o->begin();
            typename boost::graph_traits<G_t>::vertex_iterator vit, vend;
            for(boost::tie(vit, vend) = boost::vertices(_g); vit!=vend; ++vit){
                if(boost::degree(*vit, _g) == 0){
                    *r = *vit;
                    ++r;
                    ++_i;
                }else{
                }
            }
        }

    }
    ~minDegree()
    {
        if(_own_o){
            delete _o;
        }else{
        }
    }

// private: // not yet.
    void do_it()
    {
        trace2("MD", _iiv, _i);
        if(!boost::num_vertices(_g)){
            unreachable(); // caller cannot know yet.
            incomplete();
            return;
        }else{
        }
        bag_type bag_i;
        bag_type* bags_i=&bag_i;
        assert(_o);
        O_t& elim_vertices = *_o;

        vertices_size_type num_vert=boost::num_vertices(_g);

#ifndef NDEBUG
        check(_g);
#endif

        unsigned int i = _i;
        unsigned min_ntd = 1; // minimum nontrivial vertex degree
        unsigned upper_bound = 0; // computed, if T

        // constructor?
        if(_t){ untested();
            assert(elim_vertices.size() == num_vert);
            auto zerodegbag=MOVE(_degs.detach_bag(0));
            BOOST_AUTO(it, zerodegbag.begin());

            if(_iiv){ untested();
            }else{ untested();
                for(; it!=zerodegbag.end(); ++it){
                    elim_vertices[i++] = get_vd(_g, *it);
                }
            }
        }

        //trace1("entering MD", boost::num_edges(_g));
        while(boost::num_edges(_g) > 0){ untested();
            INTERRUPTION_POINT;
            assert(min_ntd != num_vert);

            // recompute ntd can only increase from here
            vertex_descriptor c;
            if(min_ntd>1){ untested();
                --min_ntd;
            }
            boost::tie(c, min_ntd) = _degs.pick_min(min_ntd, num_vert);
            assert(min_ntd == boost::degree(c, _g));

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(min_ntd >= _ub_in){ untested();
                _t->clear();
                throw exception_unsuccessful();
            }

            if(_o){
                elim_vertices[i] = get_vd(_g, c);
            }
            if(_t){ untested();
                assert(i<_bags.size());
                bags_i = &_bags[i];
            }else if(min_ntd > upper_bound){ untested();
                upper_bound = min_ntd;
            }
            assert(bags_i);

            ++i; // number of nodes in tree decomposition tree

            adjacency_iterator I, E;
            for(boost::tie(I, E) = boost::adjacent_vertices(c, _g); I!=E; ++I){ untested();
                assert(*I!=c); // no self loops...
                vertex_descriptor w=*I;
                _degs.unlink(w);
            }

            make_clique_and_detach(c, _g, *bags_i);
#ifndef NDEBUG // safety net.
            check(_g);
#endif

            redegree(NULL, _g, *bags_i, _degs);
            if(!_t){ untested();
                bags_i->clear();
            }

            _degs.unlink(c, min_ntd);

            assert(boost::degree(c, _g)==0);

            _degs.flush();
        }
        assert(boost::num_edges(_g)==0);
        assert(boost::num_vertices(_g));

#if 0 // elimination_ordering
        BOOST_AUTO(it, cdegs[0].begin());
        for(; it!=cdegs[0].end(); ++it){ untested();
            assert(i<elim_vertices.size());
            elim_vertices[i++] = get_vd(_g, *it);
        }
#endif

        //Build a treedecomposition.
        if(_t){
            if(!_iiv){ untested();
          //      i = num_vert;
            }
        }
        _i = i;

        trace1("MD done", upper_bound);
        _ub = upper_bound;

    } // do_it
public:
    void reset()
    { untested();
        assert(!_iiv); // incomplete, unneeded...?
        *_t = T_t();
        _degs = degs_type(_g);
        _i=0;
    }
    void tree_decomposition()
    {
        assert(_o); // ?!
        const degs_type& cdegs(_degs);
        auto const& B=cdegs[0];
        auto it=B.begin();
        // collect isolated vertices created during
        // do_it. (initially isolated vs have been taken care of
        // conditionally).
        for(; it!=B.end(); ++it){ untested();
            assert(_i<_o->size());
            (*_o)[_i++] = get_vd(_g, *it);
        }

        assert(_t);
        treedec::detail::skeleton_to_treedec(_g, *_t, _bags, *_o, _i);
    }
    vertices_size_type get_bagsize()
    {
        return _ub;
    }
    O_t& elimination_ordering()
    {
        trace3("MD elo", _visited.size(), _o->size(), _i);
        while(_i<_o->size()){
            trace1("pm", _i);
            auto v = _degs.pick_min(0, 0, true).first;
            unsigned int pos = get_pos(v, _g);
            if(_visited.size() && _visited[pos]){ untested();
                // ignore this vertex...
            }else{
                (*_o)[_i] = v;
                ++_i;
            }
        }

        assert(_o);
        return *_o;
    }
private:
    G_t& _g;
    T_t* _t;
    O_t* _o;
    bool _own_o;

    vertices_size_type _ub_in;
    vertices_size_type _ub;
    bool _iiv;
    size_t _i;
    std::vector<bag_type> _bags; // BUG. use _t;
    std::vector<bool> _visited;
    degs_type _degs;
}; // minDegree

// the fillIn heuristic.
template <typename G_t, typename T_t, typename O_t>
class fillIn{
public: // types
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename fill_chooser<G_t>::type fill_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;

public:
//    template<typename G_t>
    struct fill_update_cb : public graph_callback<G_t>{ //
        typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
        typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
        typedef typename fill_chooser<G_t>::type fill_type;

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

public: // implementation
    fillIn(G_t &G, O_t& O, bool ignore_isolated_vertices, vertices_size_type ub=-1):
        _g(G), _t(NULL), _o(&O), _ub_in(ub), _fill(_g)
    {
        _i = 0;

        unsigned int n_ = 0;

        //Mark isolated vertices as already visited.
        if(ignore_isolated_vertices){
            _visited.resize(boost::num_vertices(G));
            typename boost::graph_traits<G_t>::vertex_iterator vit, vend;
            for(boost::tie(vit, vend) = boost::vertices(_g); vit != vend; vit++){
                if(boost::degree(*vit, _g) == 0){
                    unsigned int pos = get_pos(*vit, _g);
                    _visited[pos] = true;
                }else{
                    n_++;
                }
            }
            _elim_vertices.resize(n_);
            trace3("fillIn with ignore", n_, boost::num_vertices(G), ub);
        }else{
            _elim_vertices.resize(boost::num_vertices(G));
            auto r=_elim_vertices.begin();
            typename boost::graph_traits<G_t>::vertex_iterator vit, vend;
            for(boost::tie(vit, vend) = boost::vertices(_g); vit!=vend; ++vit){
                if(boost::degree(*vit, _g) == 0){
                    *r = *vit;
                    ++r;
                    ++_i;
                }else{
                }
            }
        }

    }
    fillIn(G_t &G, T_t *T=NULL, O_t *O=NULL, vertices_size_type ub=-1):
        _g(G), _t(T), _o(O), _ub_in(ub), _fill(_g)
    {
        _i = 0;
        if (_o){ untested();
        }else{
        }

        // BUG: use _o, don't copy
        _elim_vertices.resize(boost::num_vertices(G));
    }

// private: //not yet
    void do_it()
    {
        bag_type bag_i;
        bag_type* bags_i = &bag_i;

        std::set<vertex_descriptor> refill_q;

        typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(_g);
        if(_t){
            _bags.resize(num_vert);
           // _elim_vertices.resize(num_vert);
        }else if(_o){
           // _elim_vertices.resize(num_vert);
        }


        fill_update_cb cb(&_fill, _g);

        unsigned int i = _i;
        unsigned int min_fill = -1;
        unsigned int upper_bound = 0; // computed, if T

        vertex_descriptor v;

        while(boost::num_edges(_g)){
            INTERRUPTION_POINT // BUG? use callbacks...
                //Find a vertex v such that least edges are missing for making the
                //neighbourhood of v a clique.
                //
            _fill.check();
            boost::tie(v, min_fill) = _fill.pick_min(0, -1, true);
            trace4("picked min", v, min_fill, _elim_vertices.size(), i);
            _fill.check();
            assert(is_valid(v,_g));

            BOOST_AUTO(deg, boost::degree(v, _g));

            // can happen... if there is a clique and an isolated node.
            // assert(deg);

            //Abort if the width of this decomposition would be larger than 'ub'.
            if(deg > _ub_in){
                assert(_t);
                _t->clear();
                throw exception_unsuccessful();
            }

            if(_t){
                assert(i<_bags.size());
                bags_i = &_bags[i];
                assert(i<_elim_vertices.size());
                _elim_vertices[i] = get_vd(_g, v);
            }else if(_o){
                assert(i<_elim_vertices.size());
                _elim_vertices[i] = get_vd(_g, v);
            }else{
            }

            _fill.mark_neighbors(v, min_fill);

            assert(!bags_i->size());


#ifndef NDEBUG
            size_t newedges = make_clique_and_detach(v, _g, *bags_i, &cb);
            if(newedges == min_fill){
            }else{ untested();
                assert(false); // for now.
                // something is terribly wrong.
                // or some extra-heuristics is active
            }
#else
            make_clique_and_detach(v, _g, *bags_i, &cb);
#endif
            _fill.unmark_neighbours(*bags_i);


            if(!_t){
                bags_i->clear();
                if(deg > upper_bound){
                    upper_bound = deg;
                }
            }else{
            }

            assert(boost::degree(v, _g)==0);
            ++i; // number of nodes in tree decomposition tree
        } // while(edges)
        std::cout << "c FI_ub " << upper_bound+1 << "\n";
        _i = i;

        // move to tree_decomposition...
        if(_t){
        }else{
            _upper_bound = upper_bound;
        }
        if(_o){
            // hack...
            *_o = _elim_vertices;
        }
    }
    vertices_size_type get_bagsize() const
    { untested();
        // if(!_done)do_it();?
        return _upper_bound;
    }
    T_t& tree_decomposition()
    {
        // _elim_vertices is not an elimination ordering. but sufficient here!
        trace3("fi::td", _visited.size(), _bags.size(), _i);
        assert(_t);
        treedec::detail::skeleton_to_treedec(_g, *_t, _bags, _elim_vertices, _i);
        return *_t;
    }
    O_t& elimination_ordering()
    {
        trace3("fill elo", _visited.size(), _o->size(), _i);
        while(_i<_o->size()){
            trace1("pm", _i);
            auto v = _fill.pick_min(0, 0, true).first;
            unsigned int pos = get_pos(v, _g);
            if(_visited.size() && _visited[pos]){ untested();
                // ignore this vertex...
            }else{
                (*_o)[_i] = v;
                ++_i;
            }
        }

        assert(_o);
        return *_o;
    }

private: // data
    G_t &_g;
    T_t *_t;
    O_t *_o;
    vertices_size_type _upper_bound;
    vertices_size_type _ub_in; // redundant?
    vertices_size_type _i; // redundant?
    fill_type _fill;
    std::vector<bag_type> _bags; // BUG: what's this? (use tree?!)
    O_t _elim_vertices; // TODO: external?!
    std::vector<bool> _visited;
}; // FillIn

} // namespace impl

} // namespace treedec

#endif // guard

// vim:ts=8:sw=4:et
