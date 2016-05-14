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

#include "graph.hpp"

#if __cplusplus >= 201103L
# include <unordered_set>
#endif

static bool lazy_init=true;

// temporary.
#ifndef untested
#define untested()
#endif
#ifndef itested
#define itested()
#endif
#ifndef incomplete
#define incomplete()
#endif


namespace misc {

namespace detail {
// FIXME: not here
template<class G>
struct fill_config{
    typedef typename boost::graph_traits<G>::vertex_descriptor vd_type;
#if __cplusplus < 201103L
    typedef std::set<vd_type> bag_type;
#else
    typedef std::unordered_set<vd_type> bag_type;
#endif
    // typedef stx::btree_set<vd_type> bag_type;
    static void alloc_init(size_t){
    }
    static unsigned num_threads(){return 1;}
};
} // detail
} // misc


namespace misc{ // treedec::misc? hmmm.

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
        for(auto& v:_vals){ itested();
            v.value = 0;
            v.queued = false;
            v.neighbour = false;
        }
     //   CFG::alloc_init(boost::num_vertices(g));
        vertex_iterator vIt, vEnd;
        bool foundzero=false;
        for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *vIt);
            if(boost::degree(*vIt, g)){
                size_t missing_edges=-1;

                if(foundzero){ untested();
                    q_eval(*vIt); //later.
                    assert(_vals[pos]==-1);
                    assert(contains(_eval_queue, *vIt));
                }else{
                    missing_edges = treedec::count_missing_edges(*vIt,_g);
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
        assert(noboost::is_valid(v,_g));
        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        unlink(v, _vals[pos]);
    }

public:
    void reg(const vertex_descriptor& v, size_t missing_edges)
    {
        assert(noboost::is_valid(v,_g));
        bool n=_fill.insert(std::make_pair(missing_edges,v)).second;
        assert(n);
        (void)n;

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
        _vals[pos].value = missing_edges;
        _vals[pos].queued = false;
        assert(_vals[pos].value!=-1 || _init);
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
        if(_vals[pos].neighbour){ untested();
            // it's a neighbour. don't touch.
            return;
        }else if(_vals[pos]==-1){ itested();
            // unknown. don't touch.
            assert(_vals[pos].queued);
            assert(contains(_eval_queue, v));
            return;
        }else{untested();
            // getting here, because v is a common neighbor of a new edge.
            // and not a neighbor...
            assert(_vals[pos].value>0);
            q_eval(v, _vals[pos].value-1);

            if(_vals[pos].value==0){ untested();
                reg(v, 0);
                _vals[pos].queued = false;
            }
//TODO: what if missing_edges==0?

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
        unsigned deg = boost::degree(v, _g);
        assert(def>=-1);

        //std::cerr << _init << "q_eval " << pos << " val: " << _vals[pos].value
                  //<< " def " << def << " deg " << deg << "\n";

        if(_vals[pos].queued){
            assert(_init || contains(_eval_queue, v));
        }

        if(def==-1 && _vals[pos].value==-1){ untested();
            assert(_init || contains(_eval_queue, v));
            return;
        }else if(!_vals[pos].queued){ untested();
            unlink(v, _vals[pos].value);
            _eval_queue.push_back(v);
            _vals[pos].queued = true;
            assert(contains(_eval_queue, v));
        }else{ untested();
            // hmm queued w/default already, update default...
            assert(contains(_eval_queue, v));
        }
        _vals[pos].value = def;
    }

public: // optimize.

    // for n \in neigbors(c):
    //   |X| = deg(n)-deg(c)
    //   queue new_fill(n) = old_fill(n) - old_fill(c) - |X|
    // (override in case n is incdent to a newly inserted edge)
    void q_neighbors(vertex_descriptor c, size_t cfill)
    {
        unsigned int posc = boost::get(boost::get(boost::vertex_index, _g), c);
        BOOST_AUTO(degc, boost::degree(c, _g));
        typename boost::graph_traits<G_t>::adjacency_iterator n, nEnd;

        boost::tie(n, nEnd) = boost::adjacent_vertices(c, _g);
        for(; n!=nEnd; ++n){
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *n);
            long old_fill=_vals[pos].value;
            assert(!_vals[pos].neighbour);
            _vals[pos].neighbour = true;

            if(_vals[pos].value==-1){ untested();
                // neighbor fill unknown. leave it like that.
                assert(contains(_eval_queue, *n));
                continue;
            }

            assert(old_fill>=0);
            // unlink(*n);
            BOOST_AUTO(degn, boost::degree(*n, _g));
            if(degn>=degc){
                long X = degn - degc;
                long new_fill = old_fill - cfill - X;
                if(new_fill < 0){ untested();
                    new_fill = -1;
                }
                q_eval(*n, new_fill);
            }else{ untested();
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
            _vals[pos].neighbour = false;
        }
    }

public: // picking
    // vertex_descriptor pick(unsigned fill)
    // {
    //     return *_fill[fill].begin();
    // }
    // pick a minimum fill vertex within fill range [lower, upper]
    std::pair<vertex_descriptor, fill_t> pick_min(unsigned lower=0, unsigned upper=-1, bool erase=false)
    {
        BOOST_AUTO(fp, _fill.begin());
        if(_fill.empty() || fp->first){

#if 0 // slightly slower...
        typename eq_t::const_iterator qi = _eval_queue.begin();
        typename eq_t::const_iterator qe = _eval_queue.end();
        for(; qi!=qe; ++qi){ untested();
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), *qi);
            size_t missing_edges = _vals[pos].value;

            if(!_vals[pos].queued){ untested();
                // taken out of queue, because fill==0.
                assert(missing_edges==0);
                // ignore
                continue;
            }else if(_vals[pos]==0) { untested();

            }

            if(missing_edges == -1){
                // unknown...
                missing_edges = treedec::count_missing_edges(*qi, _g);
            }else{
                assert(missing_edges == treedec::count_missing_edges(*qi, _g));
            }
            if(!missing_edges && erase){
                 // shortcut...?
                 // see below...
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
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), v);
            size_t missing_edges = _vals[pos].value;

            if(!_vals[pos].queued){ untested();
                // taken out of queue, because fill==0.
                assert(missing_edges==0);
                // ignore
                continue;
            }

            if(missing_edges == -1){
                // unknown...
                missing_edges = treedec::count_missing_edges(v, _g);
            }else{
                assert(missing_edges == treedec::count_missing_edges(v, _g));
            }
            if(missing_edges){
            }else if(erase){
                 // shortcut...
                _vals[pos].queued = false;
                return std::make_pair(v, 0);
            }else{ untested();
            }
            reg(v, missing_edges);
        }
#endif
        }else{ untested();
            // no need to process q. it's already there
        }

        assert(!_fill.empty());

        assert(lower==0); // for now.

        BOOST_AUTO(b, _fill.begin());
        assert(noboost::is_valid(b->second, _g));

        unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), b->second);
        assert(_vals[pos].value!=-1); (void)pos;
        assert(_vals[pos]==b->first);

        BOOST_AUTO(p, std::make_pair(b->second, b->first));
        assert(noboost::is_valid(p.first, _g));
        if(erase){
            unlink(p.first, p.second);
            unsigned int pos = boost::get(boost::get(boost::vertex_index, _g), p.first);
            _vals[pos].value = 0;
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
    typedef struct{
        size_t value; // fill value, -1==unknown.
        bool queued; // it is not in a bucket right now.
                     // !queued means it's in bucket _fill[value].
        bool neighbour; // it's a neighbour of the vertex currently processed
        bool operator==(size_t v) const{return value==v;}
    } status_t;
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

} //namespace misc

namespace noboost{

template<class G_t>
struct fill_chooser{
    typedef typename misc::FILL<G_t> type;
    typedef type fill_type; // transition? don't use.
};

}

#endif //guard

// vim:ts=8:sw=4:et
