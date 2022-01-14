// Lukas Larisch, 2016
// Felix Salfelder, 2016-2017, 2021
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
//
// greedy base
//
#ifndef TREEDEC_GREEDY_BASE_HPP
#define TREEDEC_GREEDY_BASE_HPP

// FIXME, rearrange, maybe move all to bits?
#include "../graph_traits.hpp"
#include "../skeleton.hpp"
#include "../generic_elimination_search_overlay.hpp"
#include "../induced_subgraph.hpp"

#ifdef DEBUG_FILL
#include "../graph_util.hpp"
#endif

namespace treedec{

namespace impl{

template <typename G_t, typename O_t,
          template<class G, class...> class CFGT_t=algo::default_config>
class greedy_base : public ::treedec::algo::draft::algo1{
private: // forbidden
    greedy_base() = delete;
    greedy_base(greedy_base&&) = delete;
    greedy_base(greedy_base const&) = delete;
    // greedy_base(){unreachable();}
public:
    typedef typename directed_view_select<G_t>::type graph_type;
    using D_t=graph_type;
    typedef typename boost::graph_traits<D_t>::edges_size_type edges_size_type;
    typedef typename boost::graph_traits<D_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<D_t>::vertices_size_type vertices_size_type;
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
    typedef typename boost::property_map<D_t, boost::vertex_index_t>::type idmap_type;
    typedef std::vector<vertices_size_type> degree_type;
    typedef boost::iterator_property_map<vertices_size_type*,
        idmap_type, vertices_size_type, vertices_size_type&> degreemap_type;
    typedef treedec::draft::NUMBERING_1<D_t> numbering_type;
    struct sgm{
        sgm(sgm const& n)
          : _n(n._n) {
        }
        sgm(numbering_type const& n)
          : _n(n) {
        }
        // which one?
        bool operator()(vertex_descriptor v) const{
            return _n.is_not_numbered(v);
        }
        bool operator[](vertex_descriptor v) const{ untested();
            return _n.is_not_numbered(v);
        }
        sgm& operator=(const sgm& o){
            assert(&_n==&o._n);
            return *this;
        }
        numbering_type const& _n;
    };
    typedef sgm member_pred_type;
    // TODO: alternative subgraph with edge deletion
    typedef INDUCED_SUBGRAPH_1<D_t, member_pred_type, degreemap_type> subgraph_type;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename std::vector<vertex_descriptor> bag_t;

protected: // construct/destruct
    //greedy_base(G_t const &g, unsigned ub, bool ignore_isolated_vertices=false)
    //  : greedy_base(g, ub, ignore_isolated_vertices) { untested();
    //      trace1("got own g", &_g);
    //      _own_g = true;
    //  }
    template<class G_maybe_const>
    greedy_base(G_maybe_const &g, unsigned ub, bool ignore_isolated_vertices=false)
      : algo1("."),
        _g(g),
        _o(NULL), _own_o(true), _ub_in(ub),
        _iiv(ignore_isolated_vertices), _i(0),
        _min(0), _ub_tw(0),
        _current_N(NULL),
        _num_vert(boost::num_vertices(_g)),
        _num_edges(treedec::num_edges(g)),
        _idmap(boost::get(boost::vertex_index, _g)),
        _numbering(g, _idmap),
        _degree(boost::num_vertices(_g)),
        _degreemap(boost::make_iterator_property_map(_degree.data(),
                                                     _idmap,
                                                     vertices_size_type())),
        _subgraph(_g, member_pred_type(_numbering), _degreemap),
        _marker(boost::num_vertices(_g))
    {
        trace2("greedy_base", _num_vert, _num_edges);
#if 0
        auto vv = boost::vertices(g);
        for(; vv.first != vv.second; ++vv.first){ untested();
        }
        auto ee = boost::edges(g);
        for(; ee.first != ee.second; ++ee.first){ untested();
            auto s = boost::source(*ee.first, g);
            auto t = boost::target(*ee.second, g);
            std::cout << "(" << s << "," << t << "),";
        }
        std::cout << "\n";
#endif

        if(_own_o){
            _o = new O_t;
        }else{ untested();
        }

        // BUG. part of subgraph constructor?!
        // yes: not necessary if we use subgraph with edge deletion.
        auto p=boost::vertices(g);
        for(; p.first!=p.second; ++p.first){
            auto d=boost::out_degree(*p.first, _g);
            _degreemap[*p.first] = d;
        }

#if 0
        //the following seems to be unnecessary
        if(_t){ untested();
            incomplete();
// _bags.resize(_num_vert);
        }else{ untested();
        }
#endif

        _o->resize(_num_vert);
    }

    virtual ~greedy_base(){
        if(_own_g){ untested();
            trace1("delete own g", &_g);
//            delete &_g;
        }else{
        }
        if(_own_o){
            delete _o;
        }else{ untested();
        }
    }

protected: // implementation
    void initialize(){
        // yuck.
        auto p=boost::vertices(_g);
        assert(_numbering.total()==0);
        unsigned checksum=0;
        for(; p.first!=p.second; ++p.first){
            auto d=boost::out_degree(*p.first, _g); // BUG
            checksum += d;
            assert(_idmap[*p.first] < _degree.size());
            _degreemap[*p.first] = d;
            if(d == 0){
                if(!_iiv){
                    // eliminate isolated vertices first
                    (*_o)[_i++] = *p.first;
                    _numbering.put(*p.first);
                    _numbering.increment();
                }else{
                    // just ignore them, so they don't show up in the order
                    assert(_num_vert);
                    --_num_vert;
                    // _zeroes.push_back(*p.first); // HACK
                }
            }else{
            }
        }
        trace1("initialised", _numbering.total());
        assert(checksum==_num_edges*2);
    }

public:
    // dump a tree decomposition into t
    template<class T>
    void get_tree_decomposition(T& t){
        size_t numbags = _numbering.total();
        trace2("gt", numbags, _i);

        t = T(_i);
        treedec::set_bagsize(t, bagsize());
        assert(_o);
        _o->resize(_i);
        assert(_i == boost::num_vertices(_g));
        get_elimination_ordering(*_o);
#ifndef NDEBUG
        for(auto x: *_o){
            trace1("order", x);
        }
        auto p = boost::vertices(_g);
        for(; p.first!=p.second; ++p.first){
            if(_numbering.is_numbered(*p.first)){
                trace2("numbered", *p.first, _numbering.get_position(*p.first));
            }
        }
#endif

        if(!_iiv){
            assert(_o->size() == boost::num_vertices(_g));
            // assert(_numbering.total() == boost::num_vertices(_g)); no. numbering is bags only
        }else{ untested();
        }

        typedef treedec::draft::SKELETON<D_t, numbering_type, O_t> skeleton_type;
        skeleton_type skel(_g, _numbering, *_o);

        // assert(_numbering.total()==skel.size()); // not in develop.
        trace1("skel?", _numbering.total());
        treedec::detail::skeleton_helper<D_t, T, skeleton_type, numbering_type>
            S(_g, t, skel, _numbering);

        trace1("pre td?", boost::num_vertices(t));
        S.do_it();
        trace1("td?", boost::num_vertices(t));
    }

public:
    vertices_size_type bagsize(){
        return _ub_tw+1;
    }
    vertices_size_type get_bagsize(){
        return _ub_tw+1;
    }

    // later
    O_t& get_elimination_ordering() {
        for (auto x: _zeroes){ untested();
            // HACK
            _o->push_back(x); // HACK
        } // HACK
        _zeroes.resize(0); // HACK
        return *_o;
    }

    template<class O>
    void get_elimination_ordering(O& o) const{
        if(!_iiv){
            assert(_i==boost::num_vertices(_g));
            assert(_o->size()==boost::num_vertices(_g));
        //     assert(_i==_numbering.total()); no. _numbering numbers bags only.
            // o.resize(boost::num_vertices(_g));
            assign(o, *_o);
        }else{untested();
            o.resize(_i);
            incomplete(); //??
        }
        // incomplete(); inefficient perhaps

        // ???
//        auto p = boost::vertices(_g);
//        for(; p.first!=p.second; ++p.first){
//            if(_numbering.is_numbered(*p.first)){
//                auto pos = _numbering.get_position(*p.first);
//                assert(pos<o.size());
//                o[pos] = _idmap[*p.first];
//                assert(o[pos] == (*_o)[pos]);
//            }else{ untested();
//            }
//        }
        if(!_iiv){
            assert(_i==boost::num_vertices(_g));
        }else{untested();
        }
    }

//    virtual void initialize() = 0; ??
    virtual bool next(vertex_descriptor &){ unreachable(); return false; }
    virtual void eliminate(vertex_descriptor){ unreachable(); }
    virtual void postprocessing(){ untested(); }

#ifdef DEBUG_FILL
private: // debugging
    virtual size_t fill_cached_(vertex_descriptor) const{
        unreachable();
        return -1ul;
    }
    virtual bool fill_cached_is_lb_(vertex_descriptor) const{
        unreachable();
        return false;
    }
#endif

public:

    // greedy_base::
    void do_it(){
        check(_g);

        timer_on();

        if(!_num_vert){ untested();
            timer_off();
            return;
        }else{
        }

        assert(_o);
        O_t& o = *_o;

#ifndef NDEBUG
        check(_g);
#endif

        initialize();
        _o->resize(_num_vert);
//        assert(elim_vertices.size() == _num_vert);
        vertex_descriptor c = 0;

        // assert(_num_vert == vertices_left()); // no. have already removed isolated nodes.

        // auto cnt = _num_vert;

        trace3("greedy_base::do_it", _i, _num_vert, _numbering.total());
        trace2("main loop", _num_vert, c);
        while(next(c)){
            trace2("greedy. next is", _i, c);

#if 0 // maybe later.
            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){ untested();
                throw exception_unsuccessful();
            }else{ untested();
            }
#endif

            o[_i] = c;

            size_t degc = boost::out_degree(c, _subgraph);
            if(degc > _ub_tw){
                _ub_tw = boost::out_degree(c, _subgraph);
            }else{
            }

            // assert(bags_i);?!
#ifdef DEBUG_FILL
            typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
            marker_type debug_marker(_num_vert);
            std::vector<vertex_descriptor> c_neigh;
            {
                auto cr = boost::adjacent_vertices(c, _subgraph);
                for(; cr.first!=cr.second; ++cr.first){
                    auto n = *cr.first;
                    auto me = treedec::count_missing_edges(n, debug_marker, _subgraph);
                    trace4("DEBUG_FILL pre-elim", n, fill_cached_(n), fill_cached_is_lb_(n), me);
                    if(fill_cached_is_lb_(n)){
                        assert(me >= fill_cached_(n));
                    }else{
                        assert(me == fill_cached_(n));
                    }
                    c_neigh.push_back(*cr.first);
                }
                assert(c_neigh.size()==degc);
            }
#endif
            eliminate(c);

#ifdef DEBUG_FILL
            {
                for(auto n : c_neigh){
                    auto me = treedec::count_missing_edges(n, debug_marker, _subgraph);
                    trace4("DEBUG_FILL post-elim", n, fill_cached_(n), fill_cached_is_lb_(n), me);
                    if(fill_cached_is_lb_(n)){
                        assert(me >= fill_cached_(n));
                    }else{
                        assert(me == fill_cached_(n));
                    }
                }
                assert(c_neigh.size()==degc);
            }
#endif
            ++_i;
            assert(_numbering.total()==_i);
        }

        // BUG
        trace2("done loop", _i, _num_vert);
        assert(_numbering.total()==_i);

        // add one more node maybe.
        // should list remaining clique, not implemented
        postprocessing();

        timer_off();

    } // do_it

    void set_ignore_isolated(bool v=true){
        _iiv = v;
    }

protected:
    size_t vertices_left() const{
        return _num_vert - _numbering.total();
    }

protected:
    D_t _g;
    //T_t* _t;
    O_t* _o;
    bool _own_o;
    bool _own_g{false};

    vertices_size_type _ub_in;
    bool _iiv{false};
    size_t _i;
    unsigned _min;

    vertices_size_type _ub_tw;

    // bag_t bag_i;
    bag_t* _current_N; // BUG

    unsigned _num_vert;
    edges_size_type _num_edges; // hmm, somehow use _subgraph
    idmap_type _idmap;
    numbering_type _numbering;
    degree_type _degree;
    degreemap_type _degreemap;
    subgraph_type _subgraph;
    marker_type _marker; // here? recheck: need another marker in fill?
    std::vector<vertex_descriptor> _zeroes;
}; // greedy_base

} // impl

} // treedec

#endif
// vim:ts=8:sw=4:et
