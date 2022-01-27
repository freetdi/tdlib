// Felix Salfelder 2016-2017, 2021
// Lukas Larisch 2016
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

#ifndef TREEDEC_FILL_HPP
#define TREEDEC_FILL_HPP

#include <boost/graph/graph_traits.hpp>
#include <assert.h>

#include "trace.hpp"
#include "graph_util.hpp"

#if __cplusplus >= 201103L
# include <unordered_set>
#endif

static bool const obs_lazy_init=true;
static bool const lazy_init=false;
static bool const obs_catch_zeroes_in_decrement=true;
static bool const catch_zeroes_in_decrement=false;

namespace treedec{

namespace detail {
template<class G>
struct fill_config{
    typedef typename boost::graph_traits<G>::vertex_descriptor vd_type;
#if __cplusplus < 201103L
    typedef std::set<vd_type> bag_type;
#else
    typedef std::unordered_set<vd_type> bag_type;
#endif
    static void alloc_init(size_t)
    { untested();
    }
    static unsigned num_threads(){return 1;}
};

} // detail

namespace pending {

template<class G_t, class CFG=detail::fill_config<G_t> >
class FILL{
public: // types
    // typedef CFGT<G_t> CFG;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::edges_size_type edges_size_type;
    typedef edges_size_type size_type;
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
public: // types
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename CFG::bag_type bag_type;
    typedef typename bag_type::iterator bag_iterator;
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::const_type idmap_type;
private: // types
    class status_t{
    public: // types
        typedef size_type value_type;
    public:
        status_t() : _value(0), _lb(false) {}
    public:
        void operator=(size_type x) { untested();
            _value = x;
        }
        operator size_type&() {
            return _value;
        }
        operator size_type const&() const { untested();
            return _value;
        }
    public:
        size_type get_value() const { untested();
            assert(!_lb);
            return _value;
        }
        // hmm called from shift
        size_type& value() {
            return _value;
        }
        size_type const& value() const {
            return _value;
        }
        void set_lb(bool x=true) {
            // TODO: negative value indicates LB?
            _lb = x;
        }
        void set_lb(size_t v) { untested();
            _value = v;
            _lb = true;
        }
        void update_value(size_t v) { untested();
            // same ass assign?
            _value = v;
        }
        void set_value(size_t v) {
            _value = v;
            _lb = false;
        }
        bool is_lb() const{
//            return _value < 0;
            return _lb;
        }
        bool is_known() const{ untested();
            // return _value >= 0;
            return !_lb;
        }
    private:
        size_type _value; // fill value, -1==unknown. INTMAX==unknown?
        bool _lb; // the number is just a lower bound.
    };
    // FIXME: use vertex_size_type?
private:
    typedef boost::bucket_sorter<size_type, edges_size_type,
           boost::iterator_property_map<status_t*, idmap_type, size_type, size_type&>,
           idmap_type >
        container_type;
    // typedef typename container_type::iterator iterator;
    // typedef typename container_type::const_iterator const_iterator;
    typedef typename std::make_unsigned<
                          typename boost::graph_traits<G_t>::edges_size_type
                          >::type fill_t;
    typedef typename std::make_signed<fill_t>::type offset_t;
private:
    fill_t max_missing_edges() const {
        // BUG: ask CFG
        auto mm = boost::num_edges(*_g); // treedec::num_edges?
        return mm;
    }
public: // construct
    FILL(const G_t& g, unsigned nv /*bug*/, bool init_=true)
       : _g(g),
         //_degree(boost::num_vertices(g), 0),
         _vi(boost::get(boost::vertex_index, g)),
         _vals(nv),
         _max_fill(max_missing_edges()),
         _fill(nv, // length
               _max_fill+1, // number of buckets.
               boost::make_iterator_property_map(&_vals[0], _vi, size_type()),
               _vi),
         _neigh_marker(nv)
    {
        assert(_fill.size() == _max_fill+1);
        if(init_){
            // bug
            init(g);
        }else{
        }
    } // FILL
public:
    template<class GG>
    void init(GG const& g){
        size_t nv = boost::num_vertices(_g);
        trace2("FILL", nv, _max_fill);
        auto idmap=boost::get(boost::vertex_index, g);

        _init = true;
        bool foundzero=false;
        unsigned checksum=0;
        trace1("FILL", nv);
        for( auto p = boost::vertices(g); p.first!=p.second; ++p.first){
            ++checksum;
            auto pos=boost::get(idmap, *p.first);
            trace2("init", *p.first, pos);
            assert(_vals.size()>pos);
            (void) pos;
            auto v = *p.first;
            auto deg = boost::out_degree(v, g);

            if(deg==0){
//            }else if(deg==1){ untested();
            }else{
                size_t missing_edges=-1;

                if(foundzero){ untested();
                    incomplete();
                    // bypass init. already found a nt. clique
                    q_eval(v); //later.
                }else{
                    // BUG: does not work without neigh_marker.
                    missing_edges = treedec::count_missing_edges(v, _neigh_marker, g);
                    trace3("init reg", v, deg, missing_edges);
                    reg(v, missing_edges);
                    assert(_vals[pos]==missing_edges);
                }

                if (!missing_edges){
                    trace2("found isol", v, lazy_init);
                    // faster by a few percent. sometimes?
                    foundzero = lazy_init;
                }else{
                }
            }
        }
        assert(nv==checksum);
        _init = false;
    }
#if 0
    void update(const vertex_descriptor& v, size_t missing_edges)
    { untested();
        if(missing_edges<_max_fill){ untested();
        }else{ untested();
            missing_edges=_max_fill;
        }
        assert(treedec::is_valid(v,_g));

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos] = missing_edges;
        _fill.update(v);
        _vals[pos].queued = false; // (not yet) obsolete
        assert(!_vals[pos].is_unknown() || _init);
    }
#endif
    void reg(const vertex_descriptor& v, size_t missing_edges)
    {
        auto pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos].set_value(missing_edges);
        if(missing_edges<=_max_fill){
            _fill.push(v); // update reverse map
        }else{ untested();
            //std::cerr<<"found " << missing_edges << " for " << v << "\n";
            //incomplete();
            _vals[pos].set_value(_max_fill);
            _fill.push(v); // update reverse map
            _vals[pos].set_value(missing_edges);
        }
        assert(treedec::is_valid(v,_g));

        trace2("push", pos, missing_edges);
    } // reg
public:
    // called on 2 neighbours that are not 1 neighbour
    void decrement_fill(const vertex_descriptor v) {
        auto idmap=boost::get(boost::vertex_index, _g);
        auto pos=boost::get(idmap, v);
        if(_neigh_marker.is_marked(pos)){
            // it's a neighbour of c. don't touch.
            return;
//        else if(_vals[pos]==_max_fill) untested();
            // BUG.
        }else if(_vals[pos] > _max_fill){ untested();
            --_vals[pos]; // stay in same bucket.
        }else if(_vals[pos]){
            assert(_vals[pos]<=_max_fill);
            --_vals[pos];
            _fill.update(v);
        }else{ untested();
            // unreachable(); is reachable!
        }

        if(catch_zeroes_in_decrement && _vals[pos]==0){ untested();
            incomplete();
            // fill.remove(v);
            // }else{ untested();
            //     q_eval(v, _vals[pos].get_value()-1);
            //
        }else{
        }
    }

    void q_eval(const vertex_descriptor v){
        auto idmap = boost::get(boost::vertex_index, _g);
        auto pos = boost::get(idmap, v);
        _vals[pos].set_lb();
        // _fill.update_back(v);
    }
    void update(const vertex_descriptor v){ untested();
        auto idmap = boost::get(boost::vertex_index, _g);
        auto pos = boost::get(idmap, v);
        auto const& value = _vals[pos].value();
        if(value < _min_bucket){ untested();
            _min_bucket = value;
        }else{ untested();
        }
        if(v.is_lb()){ untested();
            _fill.update_back(v); // front?
        }else{ untested();
            _fill.update(v);
        }
    }
    void prefer(const vertex_descriptor& v) { untested();
        _fill.update_front(v);
    }


    void shift(vertex_descriptor v, offset_t offset) {
        auto idmap = boost::get(boost::vertex_index, _g);
        auto pos = boost::get(idmap, v);
        auto& value = _vals[pos].value();
//        trace3("shift", pos, offset, _vals[pos].is_lb());
        if(offset >= 0){
//            trace2("positive ", v, offset);
            value += offset;
        }else if(long(value) < -offset){ untested();
            trace2("zero pad ", v, offset);
            value = 0;
            _fill.update(v); // why?
        }else{
            trace2("other ", v, offset);
            assert(offset<0);
            value += offset;
            if(value>=_max_fill){ itested();
               // still gt ..
            }else{ itested();
                _fill.update(v);
            }
        }
//        trace3("shifted ", pos, value, _vals[pos].is_lb());

        if(is_lb(v)){
        }else{
            assert(_min_bucket <= value);
        }
    }
    bool is_lb(vertex_descriptor v) const{
        auto pos = boost::get(boost::get(boost::vertex_index, _g), v);
        return _vals[pos].is_lb();
    }
    size_type const& get_value(vertex_descriptor v) const{
        auto pos = boost::get(boost::get(boost::vertex_index, _g), v);
        assert(pos < _vals.size());
        return _vals[pos].value();
    }

    // hack.
//    size_type& value(vertex_descriptor v){ untested();
//        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
//        return _vals[pos].value();
//    }
public: // O(1) neighbor stuff.
    void mark(vertex_descriptor v){
        auto idmap = boost::get(boost::vertex_index, _g);
        auto pos = boost::get(idmap, v);
        _neigh_marker.mark(pos);
    }
    void clear_marker(){
        _neigh_marker.clear();
    }
    // for n \in neigbors(c):
    //   |X| = deg(n)-deg(c)
    //   queue new_fill(n) = old_fill(n) - old_fill(c) - |X|
    // (override in case n is incdent to a newly inserted edge)
    // pending::fill::
    void mark_neighbours(vertex_descriptor c, size_t /*cfill*/) {
        trace1("newFILL::mark_neighbours", c);
        assert(c < boost::num_vertices(_g));
        _neigh_marker.clear();
        auto idmap=boost::get(boost::vertex_index, _g);
        auto posc=boost::get(idmap, c);
        (void) posc;

        // auto degc=boost::out_degree(c, _g);
        typename boost::graph_traits<G_t>::adjacency_iterator n, nEnd;

//        assert(treedec::is_valid(c,_g)); // incomplete
        assert(c < boost::num_vertices(_g));
        auto p = boost::adjacent_vertices(c, _g);

        // mark and propagate to neighbours
        for(; p.first!=p.second; ++p.first){
            auto idmap = boost::get(boost::vertex_index, _g);
            auto const& n=*p.first;
            trace1("newFILL::mark_neighbour", n);
            auto pos = boost::get(idmap, n);
            _neigh_marker.mark(pos);

        }
    } // pending::fill::mark_neighbours
private:
//    typename container_type::const_stack eval_queue() const{ untested();
//        return _fill[EVALQ_BUCKET];
////        return _eval_queue;
//    }

    template<class B>
    vertex_descriptor find_min_in_bucket(B const& b) {
        fill_t min = -1;
        auto idmap = boost::get(boost::vertex_index, _g);
        assert(!b.empty());
        vertex_descriptor ret = *b.begin();

        for(vertex_descriptor v : b){ untested();
            auto pos = boost::get(idmap, v);
            auto value = _vals[pos].value();
            if(value < min) { untested();
                min = value;
                ret = v;
            }
        }

        auto pos = boost::get(idmap, ret);
        if(_vals[pos].is_lb()) { untested();
            incomplete();
        }else{ untested();
        }

        return ret;

    }

    template<class B>
    typename B::const_iterator find_in_bucket(B const& b, size_t req_fill){
        // look out for a node with fill in bucket b.
        // nodes in bucket b are
        //  - of fill req_fill => use that
        //  - of fill > req_fill => update (move to bucket)
        //  - only lb => recompute.

        auto idmap = boost::get(boost::vertex_index, _g);

        while(!b.empty()){
            auto f = b.front();
            auto pos = boost::get(idmap, f);

            if(_vals[pos].is_lb()) { untested();
                // BUG: does not work without neigh_marker.
                auto me = treedec::count_missing_edges(f, _neigh_marker, _g);
                // auto me = treedec::count_missing_edges(f, _g);
                trace2("fill done count", f, me);
                if(me<=req_fill){ untested();
                    assert(me==req_fill);
#ifdef DEBUG
         //           _vals[pos].set_lb(false);
#endif
                    trace2("lb in Q confirmed", pos, me);
                    return b.begin();
                }else if(me>_max_fill){ untested();
                    trace4("lb in Q max", pos, me, _vals[pos].value(), _vals[pos].is_lb() );
                    _vals[pos].set_value(_max_fill);
                    _fill.update(f);
                    _vals[pos].set_value(me);
                    trace5("lb in Q max", pos, me, _vals[pos].value(), _vals[pos].is_lb(), _max_fill );

#ifndef NDEBUG
                    assert(_max_fill < _fill.size());
                    for(auto x : _fill[_max_fill]){ untested();
                        trace1("found", x);
                    }
#endif
                }else{ untested();
                    // req_fill < me <= _max_fill
                    trace2("lb in Q new", pos, me);
                    _vals[pos].set_value(me);
                    _fill.update(f);
                }
            }else if(_vals[pos].value()<=req_fill){
                trace3("found ex in Q", pos, _vals[pos], req_fill);
                assert(_vals[pos].value() == req_fill);
                return b.begin();
            //}else if(_vals[pos].value()>_max_fill){ untested();
            //    return b.begin();
            }else{
                trace3("rebucket", f, _vals[pos], req_fill);

                auto value = _vals[pos].value();

                if(value > _max_fill){ untested();
                    _vals[pos].set_value(_max_fill);
                    _fill.update(f);
                    _vals[pos].set_value(value);
                }else{
                    _fill.update(f);
                }
            }
        }
        return b.end();
    }
public: // picking
    // vertex_descriptor pick(unsigned fill)
    // { untested();
    //     return *_fill[fill].begin();
    // }
    // pick a minimum fill vertex within fill range [lower, upper]
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower=0,
            unsigned upper=-1u, bool erase=false)
    {
        auto idmap=boost::get(boost::vertex_index, _g);
        if(upper!=-1u){
            // incomplete();
        }else{
        }
        trace4("pickmin", erase, _max_fill, lower, upper);

        if(upper>_max_fill){
            // incomplete();
            upper = _max_fill;
        }else{
        }


        assert(lower==0); // for now.

#ifdef DEBUG_FILL
        for(unsigned b=0; b<=_max_fill; ++b){ itested();
            assert(b < _fill.size());
            auto bucket_min = _fill[b];
            for(auto x : bucket_min){
                auto me = treedec::count_missing_edges(x, _g); // _neigh_marker, _g);
                trace4("check", b, x, _vals[x], me);
                assert(b<=me);
                if(b<=_vals[x]){
                }else if(b<=me){ untested();
                }else{ untested();
                    assert(false);
                }

            }
        }

        unsigned min_fill=0;
        for(; min_fill<_min_bucket; ++min_fill){ untested();
            assert(_fill[min_fill].empty());
        }
#else
        unsigned min_fill = _min_bucket;
#endif

        bool found = false;
        vertex_descriptor c;
        assert(min_fill < _fill.size());
        assert(max_fill() < _fill.size());

//        trace2("================ next elim?", min_fill, _max_fill);
        while (!found) { itested();
            assert(min_fill<_max_fill);
            assert(min_fill < _fill.size());

            auto bucket_min = _fill[min_fill];
            auto it = find_in_bucket(bucket_min, min_fill);
            if(it == bucket_min.end()){
                // next bucket.
                ++min_fill;
                if(min_fill==max_fill()){ untested();
                    auto bucket_min = _fill[min_fill];
                    assert(!bucket_min.empty());

                    c = find_min_in_bucket(bucket_min);
                    found = true;
                }else{
                }
            }else{
                c = *it;
                found = true;
            }
        }
        auto pos = boost::get(idmap, c);
        auto value = _vals[pos].value();
        if(value != min_fill){
            std::cerr << "picking out of range " << c
                      << " " << value << " " << min_fill << "\n";
        }else{
        }
        trace2("================ next elim", c, _vals[pos].value());

        assert(treedec::is_valid(c, _g));

        // assert(pos<boost::num_vertices(_g)); no, if _g is a supergraph...
        (void)pos;
        if(!erase){ untested();
        }else if(value <= max_fill()){
            _fill.remove(c);
        }else if(_vals[pos].is_lb()){ untested();
            // is this needed?
            _vals[pos].set_value(max_fill());
            _fill.remove(c);
            _vals[pos].set_lb(value);
        }else{ untested();
            // is this needed?
            _vals[pos].set_value(max_fill());
            _fill.remove(c);
            _vals[pos].set_value(value);
        }

        // not really. why?!
        // assert(treedec::count_missing_edges(p.first,_g) == p.second);
        return std::make_pair(c, value);
    } // pick_min

#if 0
    size_t num_nodes() const
    { untested();
        unsigned N=0;
        for(const_iterator i=_fill.begin(); i!=_fill.end(); ++i) { itested();
            N+=i->size();
        }
        return N;
    }
#endif
    void remove(vertex_descriptor v) {
       return _fill.remove(v);
    }
    marker_type const& marked() const{
       return _neigh_marker;
    }
    size_t max_fill() const{
        return _max_fill;
    }

#if 0
    bag_type const& operator[](size_t x) const
    { untested();
        return _fill[x];
    }
    size_t size() const
    { untested();
        return _fill.size();
    }
#endif

#if 0
public:
    bag_type& operator[](size_t x)
    { untested();
        return _fill[x];
    }
#endif

private:
    bool _init; // initializing.
    const G_t& _g;
    idmap_type _vi;
//private: // later.
    std::vector<status_t> _vals;
    fill_t _min_bucket{0};
    fill_t _max_fill;
    container_type _fill;

//    mutable std::set<vertex_descriptor> _eval_queue;
    typedef std::vector<vertex_descriptor> eq_t;
//    mutable eq_t _eval_queue;

    marker_type _neigh_marker; // used in pick_min and in eliminate
}; // FILL

} // pending


template<class G_t>
struct fill_chooser{
    typedef typename treedec::pending::FILL<G_t> type;
};

} // treedec

#include "fill_obsolete.hpp"

#endif //guard
// vim:ts=8:sw=4:et
