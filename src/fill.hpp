// Felix Salfelder 2016
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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//

#ifndef TD_FILL_HPP
#define TD_FILL_HPP

#include <boost/graph/graph_traits.hpp>
#include <assert.h>

#include "trace.hpp"
#include "graph_util.hpp"

#if __cplusplus >= 201103L
# include <unordered_set>
#endif

static bool const lazy_init=true;
static bool const catch_zeroes_in_decrement=true;

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
    {
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
            _lb = x;
        }
        void set_lb(size_t v) { untested();
            _value = v;
            _lb = true;
        }
        void set_value(size_t v) {
            _value = v;
            _lb = false;
        }
        bool is_lb() const{
            return _lb;
        }
        bool is_known() const{ return !_lb; }
    private:
        size_type _value; // fill value, -1==unknown.
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
    typedef typename boost::graph_traits<G_t>::vertices_size_type fill_t;
private:
    size_t max_missing_edges() const {
        size_t nv=_vals.size();

        long mm=nv*sqrt(nv);
        //CFG::message(bLOG, "FI max missing: %l\n", mm)
        std::cerr << "FI max missing " << mm << "\n";
        return mm;
    }
public: // construct
    FILL(const G_t& g, unsigned nv /*bug*/)
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
        trace2("FILL", nv, _max_fill);
        auto idmap=boost::get(boost::vertex_index, _g);

        _init = true;
        bool foundzero=false;
        auto p=boost::vertices(g);
        unsigned checksum=0;
        trace2("FILL", nv, p.first!=p.second);
        for(; p.first!=p.second; ++p.first){
            ++checksum;
            auto pos=boost::get(idmap, *p.first);
            trace2("init", *p.first, pos);
            assert(_vals.size()>pos);
            (void) pos;
            auto n=*p.first;
            auto deg=boost::out_degree(n, g);

            if(deg==0){ untested();
//            }else if(deg==1){ untested();
            }else{
                size_t missing_edges=-1;

                if(foundzero){ untested();
                    incomplete();
                    // bypass init. already found a nt. clique
                    q_eval(*p.first); //later.
                }else{
                    missing_edges = treedec::count_missing_edges(*p.first, _neigh_marker, _g);
//                    trace2("init reg", *p.first, missing_edges);
                    reg(*p.first, missing_edges);
                    assert(_vals[pos]==missing_edges ||
                           _vals[pos]==_max_fill);
                }

                if (!missing_edges){ itested();
                    trace2("found isol", *p.first, lazy_init);
                    // faster by a few percent. sometimes?
                    foundzero = lazy_init;
                }else{
                }
            }
        }
        assert(nv==checksum);
        _init = false;
    } // FILL
public:
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
        if(missing_edges<_max_fill){
        }else{ untested();
            //std::cerr<<"found " << missing_edges << " for " << v << "\n";
            //incomplete();
            missing_edges=_max_fill;
        }
        assert(treedec::is_valid(v,_g));

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        trace2("push", pos, missing_edges);
        _vals[pos].set_value(missing_edges);
        _fill.push(v); // updates reverse map.
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
        }else if(_vals[pos]){
            assert(_vals[pos]<=_max_fill);
            --_vals[pos];
            _fill.update(v);
        }else{
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
        auto idmap=boost::get(boost::vertex_index, _g);
        auto pos=boost::get(idmap, v);
        _vals[pos].set_lb();
    }
    void update(const vertex_descriptor v){
        _fill.update(v);
    }

    void shift(vertex_descriptor v, long /*?*/ offset) {
        auto idmap=boost::get(boost::vertex_index, _g);
        auto pos=boost::get(idmap, v);
        auto& value=_vals[pos].value();
        if(offset >= 0){
            trace2("positive ", v, offset);
            value += offset;
        }else if(long(value) < -offset){
            trace2("zero pad ", v, offset);
            value = 0;
            _fill.update(v);
        }else{
            trace2("other ", v, offset);
            assert(offset<0);
            value += offset;
            if(value>_max_fill){
                value = _max_fill;
                q_eval(v); // set lb
            }else{
            }
            _fill.update(v);
        }

        assert(value<(1<<24));
    }
    size_type const& get_value(vertex_descriptor v){
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        return _vals[pos].value();
    }

    // hack.
//    size_type& value(vertex_descriptor v){
//        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
//        return _vals[pos].value();
//    }
public: // O(1) neighbor stuff.
    // for n \in neigbors(c):
    //   |X| = deg(n)-deg(c)
    //   queue new_fill(n) = old_fill(n) - old_fill(c) - |X|
    // (override in case n is incdent to a newly inserted edge)
    void mark_neighbours(vertex_descriptor c, size_t /*cfill*/)
    {
        assert(c < boost::num_vertices(_g));
        _neigh_marker.clear();
        auto idmap=boost::get(boost::vertex_index, _g);
        auto posc=boost::get(idmap, c);
        (void) posc;

        // auto degc=boost::out_degree(c, _g);
        typename boost::graph_traits<G_t>::adjacency_iterator n, nEnd;

//        assert(treedec::is_valid(c,_g)); // incomplete
        assert(c < boost::num_vertices(_g));
        auto p=boost::adjacent_vertices(c, _g);

        // mark and propagate to neighbours
        for(; p.first!=p.second; ++p.first){
            auto idmap=boost::get(boost::vertex_index, _g);
            auto const& n=*p.first;
            auto pos=boost::get(idmap, n);
            _neigh_marker.mark(pos);

//            auto old_fill=_vals[pos].get_value();
//            assert(old_fill>=0);

#if 0 // does not work
            auto degn=boost::out_degree(n, _g);
            long X=long(degn)-long(degc);
            long offset=-X-long(cfill);

//            _vals[pos].shift(offset);
            shift(n, offset);

            if(offset<0){
                // _fill.update_back(n); not yet
                _fill.update(n);
            }else{
                //?  _fill.update(n);
            }
#endif
        }
    } // mark_neighbours
private:
//    typename container_type::const_stack eval_queue() const{ untested();
//        return _fill[EVALQ_BUCKET];
////        return _eval_queue;
//    }

    template<class B>
    typename B::const_iterator find_in_bucket(B const& b, size_t req_fill){
        // look out for a node with fill in bucket b.
        // nodes in bucket b are
        //  - of fill req_fill => use that
        //  - of fill > req_fill => update (move to bucket)
        //  - only lb => recompute.

        while(!b.empty()){
            auto f=b.front();
            auto idmap=boost::get(boost::vertex_index, _g);
            auto pos=boost::get(idmap, f);

            if(_vals[pos].value()==_max_fill){
                incomplete();
            }else if(_vals[pos].is_lb()){
                auto me = treedec::count_missing_edges(f, _neigh_marker, _g);
                if(me==req_fill){
#ifdef DEBUG
                    _vals[pos].set_lb(false);
#endif
                    return b.begin();
                }else{
                    // rebucket
                    if(me<=req_fill){
#ifdef DEBUG
                        std::cerr << pos << ": missing:" << me << " req " << req_fill 
                           << " queued " << _vals[pos]  << "\n";
#endif
                        assert("wrong fill" && 0);
                    }
                    if(me>_max_fill){
                        me = _max_fill;
                        _vals[pos].set_lb();
                    }
                    _vals[pos].set_value(me); // unsets lb
                    _fill.update(f);
                }
            }else{
                return b.begin();
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
            //incomplete();
            upper = _max_fill;
        }else{
        }


        assert(lower==0); // for now.

#ifndef NDEBUG
        for(unsigned b=0; b<=_max_fill; ++b){
            auto bucket_min = _fill[b];
            for(auto x : bucket_min){
                auto me = treedec::count_missing_edges(x, _neigh_marker, _g);
                if(b<=_vals[x]){
                }else if(b<=me){
                }else{
                    trace4("check", b, x, _vals[x], me);
                    assert(false);
                }

            }
        }
#endif

        unsigned min_fill=0; // lower? maybe later?
        bool found=false;
        vertex_descriptor c;

        while (!found) {
            assert(min_fill<_max_fill);

            auto bucket_min = _fill[min_fill];
            auto it=find_in_bucket(bucket_min, min_fill);
            if(it==bucket_min.end()){
                // next bucket.
                ++min_fill;
                if(min_fill==_max_fill){ untested();
                    c = bucket_min.front();
                    found = true;
                }else{
                }
            }else{
                c = *it;
                found = true;
            }
        }
        trace2("================ next elim", c, min_fill);
//        std::cerr << "ELIM " << c << " " << min_fill << "\n";

        assert(treedec::is_valid(c, _g));

        auto pos = boost::get(idmap, c);
        (void)pos;
        if(erase){
            _fill.remove(c);
        }else{ untested();
        }

        // not really. why?!
        // assert(treedec::count_missing_edges(p.first,_g) == p.second);
        return std::make_pair(c, min_fill);
    }

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
    marker_type const& marked() const{
       return _neigh_marker;
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
    size_t _max_fill;
    container_type _fill;

//    mutable std::set<vertex_descriptor> _eval_queue;
    typedef std::vector<vertex_descriptor> eq_t;
//    mutable eq_t _eval_queue;

    marker_type _neigh_marker; // used in pick_min and in eliminate
}; // FILL

} // pending

namespace obsolete {

// old hash based fill.
template<class G_t, class CFG=detail::fill_config<G_t> >
class FILL{
    FILL(const FILL&)
    { untested();
    }
public: // types
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename CFG::bag_type bag_type;
    typedef typename bag_type::iterator bag_iterator;
    //typedef std::vector<bag_type> container_type;
    typedef std::set<std::pair<size_t, vertex_descriptor> > container_type;
    // typedef typename container_type::iterator iterator;
    // typedef typename container_type::const_iterator const_iterator;
    typedef typename boost::graph_traits<G_t>::vertices_size_type fill_t;

private:
    size_t max_missing_edges() const
    {
        // isolated nodes?
        size_t nv=boost::num_vertices(_g);
        return (nv-1)*(nv-2);
    }

public: // construct
    FILL(const G_t& g): _g(g)
    {
        _init = true;
        _vals.resize(boost::num_vertices(g));
     //   CFG::alloc_init(boost::num_vertices(g));
        vertex_iterator vIt, vEnd;
        bool foundzero=false;
        for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *vIt);
            (void) pos;
            if(boost::out_degree(*vIt, g)){
                size_t missing_edges=-1;

                if(foundzero){
                    q_eval(*vIt); //later.
                    assert(_vals[pos]==-1);
                    assert(contains(_eval_queue, *vIt));
                }else{
                    missing_edges = treedec::count_missing_edges(*vIt, _g);
                    reg(*vIt, missing_edges);
                    assert(_vals[pos]==missing_edges);
                }

                if (!missing_edges){
                    // faster by a few percent. sometimes?
                    foundzero = lazy_init;
                }
            }else{
                //skip isolated vertices
            }
        }
        _init = false;
    }
public: // check
    void check(){
#ifdef EXCESSIVE_CHECK
        for(auto p : _fill) {
            assert(treedec::count_missing_edges(p.second, _g) == p.first);
        }
#endif
    }

public: // queueing
    void unlink(const vertex_descriptor& v, size_t f)
    {
        assert(f>=0);
        int n=_fill.erase(std::make_pair(f,v));
        (void)n;
        assert(n==1 || _init);

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos].value = -1; // dequeued...
        _vals[pos].queued = false;
    }
    void unlink(const vertex_descriptor& v)
    {
        assert(treedec::is_valid(v,_g));
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        unlink(v, _vals[pos]);
    }

public:
    void reg(const vertex_descriptor& v, size_t missing_edges)
    {
        assert(treedec::is_valid(v,_g));
        bool n=_fill.insert(std::make_pair(missing_edges,v)).second;
        assert(n);
        (void)n;

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos].value = missing_edges;
        _vals[pos].queued = false;
        assert(!_vals[pos].is_unknown() || _init);
    }
public:
    void reg(const vertex_descriptor v)
    {
        size_t missing_edges=treedec::count_missing_edges(v,_g);
        reg(v, missing_edges);
    }
    void q_decrement(const vertex_descriptor v)
    {
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        if(_vals[pos].is_neighbour()){
            // it's a neighbour. don't touch.
            return;
        }else if(_vals[pos].is_unknown()){ untested();
            // unknown. don't touch.
            assert(_vals[pos].queued);
            assert(contains(_eval_queue, v));
            return;
        }else{
            // getting here, because v is a common neighbor of a new edge.
            // and not a neighbor...
            assert(_vals[pos].value>0);
            q_eval(v, _vals[pos].value-1);

            if(catch_zeroes_in_decrement &&_vals[pos].value==0){
                // stash for fill, in bucket #0
                reg(v, 0);
                // extra hack: leave in queue, but mark unqueued.
                // (it's not possible to remove from queue)
                _vals[pos].queued = false;
            }

            assert(contains(_eval_queue, v));
            return;
        }
    }

    // queue for later evaluation.
    // default value may be specified in case evalution turns out to be
    // not necessary
    void q_eval(const vertex_descriptor v, int def=-1)
    {
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        assert(def>=-1);

        if(def==-1 && _vals[pos].is_unknown()){
            assert(_init || contains(_eval_queue, v));
            return;
        }else if(!_vals[pos].queued){
            unlink(v, _vals[pos].value);
            _eval_queue.push_back(v);
            _vals[pos].queued = true;
            assert(contains(_eval_queue, v));
        }else{
            // hmm queued w/default already, update default...
            assert(contains(_eval_queue, v));
        }
        _vals[pos].value = def;
    }

public: // O(1) neighbor stuff.
    // for n \in neigbors(c):
    //   |X| = deg(n)-deg(c)
    //   queue new_fill(n) = old_fill(n) - old_fill(c) - |X|
    // (override in case n is incdent to a newly inserted edge)
    void mark_neighbors(vertex_descriptor c, size_t cfill)
    {
        typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
        unsigned int posc = boost::get(boost::get(boost::vertex_index, _g), c);
        (void) posc;
        vertices_size_type degc=boost::out_degree(c, _g);
        typename boost::graph_traits<G_t>::adjacency_iterator n, nEnd;

        boost::tie(n, nEnd) = boost::adjacent_vertices(c, _g);
        for(; n!=nEnd; ++n){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *n);
            _vals[pos].set_neighbour();

            if(_vals[pos].is_unknown()){
                // neighbor fill unknown. leave it like that.
                assert(_vals[pos].queued);
                assert(contains(_eval_queue, *n));
                continue;
            }
            long old_fill=_vals[pos].value;
            assert(old_fill>=0);

            vertices_size_type degn=boost::out_degree(*n, _g);
            if(degn>=degc){
                long X = degn - degc;
                long new_fill = old_fill - cfill - X;
                if(new_fill < 0){
                    // new fill is wrong.
                    // there must be edges to fill adjacent to n...
                    q_eval(*n);
                }else{
                    q_eval(*n, new_fill);
                }
            }else{
                q_eval(*n);
            }
        }
    }
    // faster?!
    template<class N>
    void unmark_neighbours(N const& neighs)
    {
        typename N::const_iterator i=neighs.begin();
        typename N::const_iterator e=neighs.end();
        for(; i!=e; ++i){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *i);
            _vals[pos].unset_neighbour();
        }
    }

public: // picking
    // vertex_descriptor pick(unsigned fill)
    // {
    //     return *_fill[fill].begin();
    // }
    // pick a minimum fill vertex within fill range [lower, upper]
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower=0,
            unsigned upper=-1u, bool erase=false)
    {
        if(upper!=-1u){
//            incomplete();
        }
        BOOST_AUTO(fp, _fill.begin());
        if(_fill.empty() || fp->first){

#if 1 // slightly slower...
        typename eq_t::const_iterator qi = _eval_queue.begin();
        typename eq_t::const_iterator qe = _eval_queue.end();
        assert(qe!=qi || !_fill.empty());
        for(; qi!=qe; ++qi){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *qi);
            size_t missing_edges = _vals[pos].value;

            if(!_vals[pos].queued){
                // taken out of queue, because fill==0.
                assert(missing_edges==0);
                // ignore
                continue;
            }

            if(_vals[pos].is_unknown()){
                // unknown...
                missing_edges = treedec::count_missing_edges(*qi, _g);
            }
            assert(missing_edges>=0);
            reg(*qi, missing_edges);
            assert(_vals[pos] == missing_edges);
            // assert(!contains(_eval_queue, *qi));
        }
        _eval_queue.clear();
#else // faster...? (try fifo?!)
        while(!_eval_queue.empty()){
            vertex_descriptor v=_eval_queue.back();
            _eval_queue.pop_back();
            unsigned pos = boost::get(boost::get(boost::vertex_index, _g), v);
            unsigned missing_edges = _vals[pos].value;

            if(!_vals[pos].queued){
                // taken out of queue, because fill==0.
                assert(missing_edges==0);
                // ignore
                continue;
            }

            if(missing_edges == -1u){
                // unknown...
                missing_edges = treedec::count_missing_edges(v, _g);
            }
            if(missing_edges){
            }
            else if(erase){
                 // shortcut...
                _vals[pos].queued = false;
                return std::make_pair(v, 0);
            }
            reg(v, missing_edges);
        }
#endif
        }
        else{
            // no need to process q. it's already there
        }

        assert(!_fill.empty());

        assert(lower==0); // for now.

        BOOST_AUTO(b, _fill.begin());
        assert(treedec::is_valid(b->second, _g));

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), b->second);
        (void)pos;
        assert(!_vals[pos].is_unknown());
        assert(_vals[pos]==b->first);

        BOOST_AUTO(p, std::make_pair(b->second, b->first));
        assert(treedec::is_valid(p.first, _g));
        if(erase){
            unlink(p.first, p.second);
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), p.first);
            _vals[pos].set_value(0);
            assert(!_vals[pos].queued);
        }else{ untested();
        }

        assert(treedec::count_missing_edges(p.first,_g) == p.second);
        return p;
    }

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

#if 0
    bag_type const& operator[](size_t x) const
    {
        return _fill[x];
    }
    size_t size() const
    {
        return _fill.size();
    }
#endif

#if 0
public:
    bag_type& operator[](size_t x)
    {
        return _fill[x];
    }
#endif

private:
    class status_t{
    public:
        status_t() : value(0), queued(false), _n(false) {}
    public:
        size_t value; // fill value, -1==unknown.
        bool queued; // it is not in a bucket right now.
                     // !queued means it's in bucket _fill[value].
        void set_neighbour()
        {
            assert(!_n);
            _n=true;
        }
        void unset_neighbour()
        {
            _n=false;
        }
        bool is_neighbour() const
        {
            return(_n);
        }
        bool is_unknown() const
        {
            return(value==(size_t)-1);
        }
        size_t get_value() const
        {
            assert(!is_unknown());
            return value;
        }
        void set_value(size_t v)
        {
            value = v;
        }
        bool operator==(size_t v) const{return value==v;}
    private:
        bool _n; // it's a neighbour of the vertex currently processed
    };
private:

    bool _init; // initializing.
    const G_t& _g;
//private: // later.
    container_type _fill;
    std::vector<status_t> _vals;

//    mutable std::set<vertex_descriptor> _eval_queue;
    typedef std::vector<vertex_descriptor> eq_t;
    mutable eq_t _eval_queue;
}; // FILL

} // obsolete

template<class G_t>
struct fill_chooser{
    typedef typename treedec::pending::FILL<G_t> type;
};

}

#endif //guard
// vim:ts=8:sw=4:et
