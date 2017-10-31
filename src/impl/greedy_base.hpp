// Lukas Larisch, 2016
// Felix Salfelder, 2016-2017
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
// greedy base
//
#ifndef TD_GREEDY_BASE_HPP
#define TD_GREEDY_BASE_HPP

// FIXME, rearrange, maybe move all to bits?
#include "../graph_traits.hpp"
#include "../skeleton.hpp"
#include "../induced_subgraph.hpp"

namespace treedec{

namespace impl{

#if 0
// check: why is this in experimental?
template<class G>
struct directed_view_select{
    typedef treedec::draft::directed_view<G> type;
};
#endif

template <typename G_t, typename O_t,
          template<class G, class...> class CFGT_t=algo::default_config>
class greedy_base : public ::treedec::algo::draft::algo1{
private: // forbidden
    greedy_base(){unreachable();}
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
        sgm(numbering_type const& n)
          : _n(n)
        {
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
    greedy_base(G_t &g, unsigned ub, bool ignore_isolated_vertices=false)
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
        trace2("ghb", _num_vert, _num_edges);
        if(_own_o){
            _o = new O_t;
        }else{
        }

        // BUG. part of subgraph constructor?!
        // yes: not necessary if we use subgraph with edge deletion.
        auto p=boost::vertices(g);
        for(; p.first!=p.second; ++p.first){
            auto d=boost::out_degree(*p.first, g);
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
        if(_own_o){
            delete _o;
        }else{
        }
    }

protected: // implementation
    void initialize(){
        // yuck.
        auto p=boost::vertices(_g);
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
                    _numbering.increment(); // TODO. possibly unnecessary.
                }else{
                    // just ignore them, so they don't show up in the order
                    assert(_num_vert);
                    --_num_vert;
                }
            }else{
            }
        }
        assert(checksum==_num_edges*2);
    }

public:
	 // dump a tree decomposition into t
	 template<class T>
    void tree_decomposition(T& t){
        assert(_o);
        _o->resize(_i);
        elimination_ordering(*_o);
#ifndef NDEBUG
        for(auto x: *_o){
            trace1("order", x);
        }
        auto p=boost::vertices(_g);
        for(; p.first!=p.second; ++p.first){
            if(_numbering.is_numbered(*p.first)){
                trace2("number", *p.first, _numbering.get_position(*p.first));
            }
        }
#endif

        typedef treedec::draft::SKELETON<D_t, numbering_type, O_t> skeleton_type;
        skeleton_type skel(_g, _numbering, *_o);

        assert(_i==skel.size());
        treedec::detail::skeleton_helper<D_t, T, skeleton_type, numbering_type>
            S(_g, t, skel, _numbering);
        S.do_it();
    }

public:
    vertices_size_type get_bagsize(){
        return _ub_tw+1;
    }

    // later
    O_t& get_elimination_ordering() const { untested();
        return *_o;
    }

    template<class O>
    void get_elimination_ordering(O& o) const{
		 if(_i){
		 }else{untested();
		 }
		 // incomplete(); inefficient perhaps

        o.resize(_i);
        auto p=boost::vertices(_g);
        for(; p.first!=p.second; ++p.first){
            if(_numbering.is_numbered(*p.first)){
                auto pos=_numbering.get_position(*p.first);
                assert(pos<o.size());
                o[pos] = _idmap[*p.first];
                assert(o[pos] == (*_o)[pos]);
            }else{
            }
        }
    }

//    virtual void initialize() = 0; ??
    virtual bool next(vertex_descriptor &){ unreachable(); return false; }
    virtual void eliminate(vertex_descriptor){ unreachable(); }
    virtual void postprocessing(){ untested(); }

    // greedy_base::
    void do_it(){
        trace2("do_it", _i, _num_vert);
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

        vertex_descriptor c;

        while(next(c)){

#if 0 // maybe later.
            //Abort if the width of this decomposition would be larger than 'ub'.
            if(_min >= _ub_in){ untested();
                throw exception_unsuccessful();
            }else{
            }
#endif

            elim_vertices[_i] = c;

            if(boost::out_degree(c, _subgraph)>_ub_tw){
                _ub_tw = boost::out_degree(c, _subgraph);
            }else{
            }

            // assert(bags_i);?!

            eliminate(c);

#if 0
            if(!_t){ untested();
                _current_N->clear();
            }
#endif

            ++_i;
        }

        // BUG
        trace2("done loop", _i, _num_vert);

        // add one more node maybe.
        // should list remaining clique, not implemented
        postprocessing();

        timer_off();

    } // do_it

protected:
    D_t _g;
    //T_t* _t;
    O_t* _o;
    bool _own_o;

    vertices_size_type _ub_in;
    bool _iiv;
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
}; // greedy_base

} // impl

} // treedec

#endif
// vim:ts=8:sw=4:et
