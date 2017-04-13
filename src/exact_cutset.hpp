// # define COUNTERS
// Lukas Larisch, 2014 - 2015
// Felix Salfelder 2016
//
// (c) 2014-2015 Goethe-Universität Frankfurt
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

/*
 * Offers functionality to compute a tree decomposition of exact width.
 *
 * Provides following functions:
 *
 * - void exact_cutset(G_t &G, T_t &T, int lb)
 * - void exact_cutset(G_t &G, T_t &T)
 *
 */

// TODO:
// - statically allocated treedec
// - use treedec::bag properly, not bagdraft::bag
// - move everthing to treedec
// - leaf upper bounds
// - parallel ...

#ifndef TD_EXACT_CUTSET
#define TD_EXACT_CUTSET

// #define EXCUT_USE_DELTAC

#include "algo.hpp"
#include "iter.hpp"
#include "lower_bounds.hpp"
#include "treedec.hpp"
#include "overlay.hpp"

#include <stx/btree_set>
#include <boost/container/flat_set.hpp>
#include <boost/graph/adjacency_matrix.hpp>
/// HACK
// #define EC_COMP_SET stx::btree_set
// #define EC_COMP_SET std::set
// #define EC_COMP_SET boost::container::flat_set
#define EC_COMP_SET std::vector


typedef BOOL EXCUT_BOOL;

namespace treedec{

template<class S, class T>
void check_same(S const& s, T const& t)
{
    assert(s.size()==t.size());
    BOOST_AUTO(si, s.begin());
    BOOST_AUTO(ti, t.begin());
    for(; si!=s.end(); ++si){
        assert(*si==*ti);
        ++ti;
    }
}


namespace draft{

template<class G>
bool operator!=(typename treedec::VECTOR_TD<G>::value_type& a,
                typename treedec::VECTOR_TD<G>::value_type& b)
{ untested();
    return &a != &b;
}

}//draft
} // treedec

namespace treedec{
namespace bagdraft{
    // FIXME: should be just treedec::bag. does not work yet
template<class TV>
inline typename TV::bag_type const& bag(
       TV const&,
       typename TV::const_vertex_descriptor x)
{
    return x->second;
}
template<class TV>
inline typename TV::bag_type& bag(
        TV &,
        typename TV::const_vertex_descriptor x)
{
//    return x->second;
    return const_cast<typename TV::T_vertex_descriptor>(x)->second;
}
template<class TV>
inline typename TV::bag_type& bag(
        TV &,
        typename TV::vertex_descriptor x)
{
    return x->second;
    return const_cast<typename TV::vertex_descriptor>(x)->second;
}
template<class TV>
inline typename TV::bag_type const& bag(
        TV const&,
        typename TV::vertex_descriptor x)
{
    return x->second;
    return const_cast<typename TV::vertex_descriptor>(x)->second;
}
} // bagdraft

template<class G>
unsigned get_pos(typename treedec::VECTOR_TD<G>::const_vertex_descriptor v,
                 treedec::VECTOR_TD<G> const& _g)
{
    size_t s = sizeof(typename treedec::VECTOR_TD<G>::value_type);
    return (intptr_t(v) - intptr_t(*_g.begin()))/s;
}

namespace detail{ //

template<class G>
class excut_control;
template<class G>
class excut_worker : public VECTOR_TD<G>
    { //
public: // common types
    typedef typename boost::graph_traits<G>::vertex_descriptor vd;
    typedef typename std::vector<vd>::iterator A;// bag_iter
    typedef std::vector<EXCUT_BOOL> mask_t;

public:
    typedef typename VECTOR_TD<G>::vertex_iterator T_vertex_iterator;
    typedef typename VECTOR_TD<G>::bag_type bag_type;
    typedef typename VECTOR_TD<G>::value_type value_type;
    typedef typename VECTOR_TD<G>::vertex_descriptor vertex_descriptor;
    typedef typename VECTOR_TD<G>::vertex_descriptor T_vertex_descriptor;
//    typedef value_type *vertex_descriptor;
    typedef std::pair<vertex_iterator_G, vertex_iterator_G> VRP;
    typedef typename bag_type::iterator bag_iterator;
    typedef typename bag_type::const_iterator bag_const_iterator;
// TODO don't use detail
    typedef std::pair<bag_const_iterator, bag_const_iterator> BRP; // ?!
    typedef typename ::detail::bfs_iter<G, VRP, EXCUT_BOOL>::scratch_type nrs;
    typedef typename ::detail::bfs_iter<G, VRP, EXCUT_BOOL>::scratch_type ors;
    typedef typename ::detail::bfs_iter<G, BRP, EXCUT_BOOL>::bfs_range bfs_range;
    // not yet.

    typedef typename ::detail::components_iter<G, VRP, EXCUT_BOOL>::scratch_type crs;
    typedef typename boost::graph_traits<G>::vertex_descriptor G_vertex_descriptor;
    typedef EC_COMP_SET<G_vertex_descriptor> co_t;
    typedef value_type const* const_vertex_descriptor;
    typedef const_vertex_descriptor td_cvd;
    typedef T_vertex_descriptor td_vd;

    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> excg_graph;
    // todo. table with integers...
    class job_scratch{
    public:
        ors nr; // neighbourhood range scratch
        crs cr; // candidate range scratch.
        std::vector<EXCUT_BOOL> cc; // another visited array...
    };

    class cjob_t{ //
    public:
        cjob_t(G const& g, bag_type const& b)
            : _candidates_b(g,0), _candidates_e(g, 0), _cutbag(&b)
        {
        }

    public:
        typename excut_worker<G>::bfs_range::first_type _candidates_b;
        typename excut_worker<G>::bfs_range::second_type _candidates_e;

        td_cvd parent;
        td_vd cut_ext;

        job_scratch scratch;
        std::vector<EXCUT_BOOL> visited; // some working array
        std::vector<EXCUT_BOOL> cc_mask; // full conn component mask
#ifdef COUNTERS
        unsigned cc_size;                // size of cc
#endif
        std::vector<G_vertex_descriptor> cand_cache; // candidate subset generation scratch
        typename bag_type::const_iterator cut_prefix_end;
        T_vertex_iterator cut_extip1; // hack

        bag_type const& cutbag() const{
            assert(_cutbag);
            return *_cutbag;
        }
        void attach_cutbag(bag_type const& b)
        {
            _cutbag = &b;
        }
    private:
        bag_type const* _cutbag;
    };
public: // explore
    template<class CRB, class CRI>
    bool explore_cutsets(
            CRB const& cut_red_bag,
             mask_t& ccm, CRI celt, unsigned cmps, nrs* nrs,
             td_vd cut_ext);

    template<class CRB, class CRI, class CC, class NeSc>
    bool q_explore_cutsets(CRB const & cut_red_bag,
            td_vd cut_ext,
            CC&, CRI celt, unsigned cmps, NeSc* );
private: //  common implementation
    bool try_candidate_set(cjob_t&, bool neighbours) throw();
    bool work_candidates(cjob_t& comp_job);
    template<class B>
    bfs_range compile_candidates_range(unsigned, B const& cut, const G&,
        ors* orsp, mask_t& mask);
    template<class cr, class cmp>
    bool viceatovin(td_vd cut_ext, cr& cut_red_bag, cmp& cc_scratch,
        unsigned& n);
//    T_vertex_descriptor
    template<class cvd, class cr, class tt, class NeSc>
    unsigned
    do_component(cvd& cut_ext,
                 cr& cmps_range,
                 job_scratch& scratch,
                 NeSc* nrs,
                 tt const&
#ifdef COUNTERS
                 , unsigned number_unvisited
//                 , cjob_t& comp_job
                 , bool neighs
#endif
                 ) throw();
#ifdef EXCUT_USE_DELTAC
    void add_induced_edges(unsigned cliquesize, excg_graph&);
    void add_induced_edges_square(unsigned cliquesize, excg_graph&);
#endif
private: // common stuff
     cjob_t* new_cj(unsigned bs, bag_type const& bag)
     {
         cjob_t*J;
         if( _cjs_trash.empty()) {
             J = new cjob_t(_g, bag);
             J->visited.resize(boost::num_vertices(_g));
//             J->visited_cc.resize(boost::num_vertices(_g));
             J->cc_mask.resize(boost::num_vertices(_g));
             J->scratch.cc.resize(boost::num_vertices(_g)); // replace ordered?
             J->scratch.nr.resize(bs);
             J->cand_cache.resize(bs-1);
         } else{
             J = _cjs_trash.top();
             J->attach_cutbag(bag);
             _cjs_trash.pop();
         }
         return J;
     }
     void recycle(cjob_t* trash)
     {
         _cjs_trash.push(trash);
     }

    // serial excut worker
public: // construct
    excut_worker(const G& g, size_t b)
        : VECTOR_TD<G>(g,b), _g(g), _bagsize(b)
#ifdef EXCUT_USE_DELTAC
            ,_rorigin(boost::num_vertices(_g))
#endif
    {
    }
    VECTOR_TD<G>& v()
    { untested();
        return *this;
    }
    ~excut_worker()
    {
         while(!_cjs_trash.empty()) {
             delete _cjs_trash.top();
             _cjs_trash.pop();
         }
    }
private: // common data
    const G& _g;
    unsigned _bagsize;
    /*static*/ std::stack<cjob_t*> _cjs_trash;
#ifdef EXCUT_USE_DELTAC
    excg_graph _cg;
    std::vector<G_vertex_descriptor> _origin;
    std::vector<uint16_t> _rorigin;
#endif
};

//serial
template<class G>
class excut_control{ //
public:
    typedef excut_worker<G> W;
    typedef typename W::value_type value_type;
    typedef typename W::vertex_descriptor vertex_descriptor;
    typedef typename W::const_vertex_descriptor const_vertex_descriptor;
    typedef typename W::vertex_descriptor T_vertex_descriptor;
    typedef typename W::bag_type bag_type;
    typedef T_vertex_descriptor td_vd;
public:
    excut_control(const G& g, size_t bagsize)
        : _g(g), _bagsize(bagsize), _results(g, bagsize)
    {
        _results.reserve(2*boost::num_vertices(g));
    }
    template<class C_t, class CC, class NeSc>
    void q_root_cutset(
            td_vd root, CC& ccm,
            C_t c, unsigned cmps, NeSc* nrs)
    {
        auto& cut_red_bag=bagdraft::bag(_results, root);

        // incomplete repeat "catch small" code..
        _success = _results.q_explore_cutsets(cut_red_bag, root, ccm, c, cmps, nrs);
    }
    void run() const{}
    bool join() const
    {
        return _success;
    }
private:
    const G& _g;
    size_t _bagsize;
public:  // HACK

    excut_worker<G> _results;
    bool _success;

}; // excut_control


// cut_red = vertices in cut_ext adjacent to a vertex in comp_ordered
// return false if too many
template<class G>
template<class cr, class cmp>
bool excut_worker<G>::viceatovin(td_vd cut_ext, cr& cut_red_bag, cmp& cc_mask,
        unsigned& n)
{
    BOOST_AUTO(const &ce_bag, bagdraft::bag(*this, cut_ext));
    cut_red_bag.resize(ce_bag.size());
    for(BOOST_AUTO(sIt, ce_bag.begin()); sIt!=ce_bag.end(); ++sIt)
    {
        typename boost::graph_traits<G>::adjacency_iterator nIt, nEnd;
        typename boost::graph_traits<G>::vertex_descriptor II = *sIt;

        for(boost::tie(nIt, nEnd)=boost::adjacent_vertices(II, _g); nIt!=nEnd; ++nIt){
            if(!cc_mask[get_pos(*nIt, _g)]){
                if(n+1 == _bagsize){
                    return false;
                }
                cut_red_bag[n++] = *sIt;
                break;
            }
        }
    }
    return true;
}

#ifdef EXCUT_USE_DELTAC
template<class G>
void excut_worker<G>::add_induced_edges(unsigned cliquesize, excg_graph& h)
{
    unsigned nv=_origin.size();
    unsigned compsize=nv-cliquesize;

    for(unsigned i=0; i<nv; ++i){ itested();
        boost::add_vertex(h);
    }

    for(unsigned i=0; i<compsize; ++i){ itested();
        auto A=boost::adjacent_vertices(_origin[i], _g);

        for(; A.first!=A.second; ++A.first){
            assert(_rorigin[*A.first]<nv);
            boost::add_edge(i, _rorigin[*A.first], h);
        }

    }
    for(unsigned i=nv-cliquesize; i<nv; ++i){
        unsigned j=i+1;
        for(; j<nv; ++j){
            boost::add_edge(i, j, h);
        }
    }
}
// add an edge to h if origin[edge] is an edge in g
// also insert clique
template<class G>
void excut_worker<G>::add_induced_edges_square(unsigned cliquesize, excg_graph& h)
{
    unsigned nv=_origin.size();
    while(boost::num_vertices(h)>nv){
        boost::remove_vertex(boost::num_vertices(h)-1, h);
    }
    while(boost::num_vertices(h)<nv){
        boost::add_vertex(h);
    }
    for(unsigned i=0; i<nv-1; ++i){ itested();
        unsigned j=i+1;
        for(; j<nv-cliquesize; ++j){
            if(boost::edge(_origin[i], _origin[j], _g).second){
                boost::add_edge(i, j, h);
            }
        }
    }
    for(unsigned i=nv-cliquesize; i<nv; ++i){
        unsigned j=i+1;
        for(; j<nv; ++j){
            boost::add_edge(i, j, h);
        }
    }
}
#endif // DELTA_C

template <typename G_t, class BV_t, class CRI, class NeSc>
bool explore_cutsets(G_t const &G,
         typename BV_t::const_vertex_descriptor cut_base,
         CRI cbegin, CRI cend, unsigned cmps, NeSc* nrs,
         BV_t &results,
         unsigned bagsize);

// see what cut_base cuts off, connect a TD to parent.
template <typename G_t, class BV_t, class CRI, class NeSc>
bool q_explore_cutsets(G_t const &G,
         typename BV_t::const_vertex_descriptor cut_base,
         typename BV_t::const_vertex_descriptor parent,
         CRI cbegin, CRI cend, unsigned cmps, NeSc* nrs,
         BV_t &results,
         unsigned bagsize);


// try to seperate component
// q if possible
template<class G>
template<class cvd, class cr, class tt, class NeSc>
//typename excut_worker<G>::T_vertex_descriptor
unsigned
excut_worker<G>::do_component(cvd& cut_ext,
        cr& comp_range,
        job_scratch& scratch,
        NeSc* nrs,
        tt const&
#ifdef COUNTERS
        , unsigned number_unvisited
//        , cjob_t& comp_job
        , bool neighs
#endif
        ) throw()
{

#ifdef EXCUT_USE_DELTAC
    _origin.resize(0);
    // excg_graph _cg;
#endif
    bool red_successful=true;

    // necessary?
    scratch.cc.assign(boost::num_vertices(_g), true);

    auto some_element=*comp_range.first;
    unsigned cos=0;

    // recompute component mask... really needed?
    for(; comp_range.first!=comp_range.second; ++comp_range.first){
        auto pos=get_pos(*comp_range.first, _g);
#ifdef EXCUT_USE_DELTAC
        _origin.push_back(pos);
        _rorigin[pos]=cos;
#endif
        scratch.cc[pos]=false;
        ++cos;
    }
#ifdef COUNTERS
    assert(cos<=number_unvisited);
    if(neighs){
    }else if(cos==number_unvisited){
        // there's only one component...
        // hmmm
        return false;
    }
#endif
    td_vd cut_ext2=boost::add_vertex(*this);
    auto& cut_red_bag=bagdraft::bag(*this, cut_ext2);

    unsigned redsize=0;
    // cut_red <- vertices in cut_ext adjacent to a vertex in cut-off component
    // (unsuccessful, if more than k.)
    // inefficient. merge into add_induced_edges?

    // this only works if scratch.cc exactly masks the component.
    red_successful = viceatovin(cut_ext, cut_red_bag, scratch.cc, redsize);

    if(!red_successful){
        return 0;
    }else{
        cut_red_bag.resize(redsize);
    }

#ifdef EXCUT_USE_DELTAC
//    if(3+redsize<_bagsize)
    assert(redsize<_bagsize);
    if(cos<_bagsize){
        // compon too small
    }else if(cos>5*_bagsize){
        // copy too expensive
    }else if(redsize+1==_bagsize){
        for(auto i : cut_red_bag){
            _rorigin[i] = _origin.size();
            _origin.push_back(i);
        }
        assert(_origin.size()==redsize+cos);

//        trace2("clear?", boost::num_vertices(_cg), boost::num_edges(_cg));
        assert(boost::num_edges(_cg)==0);
//        _cg.clear(); // ouch. expensive.
        add_induced_edges(redsize, _cg);

        unsigned lb_deltaC = treedec::lb::deltaC_least_c(_cg) + 1;
        if(lb_deltaC>_bagsize){ itested();
            return 0;
        }
    }else{
        // deltaC does not help a lot
    }
#endif

    bool visited_all=false; // for now.
    bool small_enough = cut_red_bag.size() + cos <= _bagsize;

//    compute upper bound
//    if small enough
//        use it.
    bool success = false;

    if(small_enough){
        auto leaf=cut_ext2;
        BOOST_AUTO(& target, bagdraft::bag(*this, leaf));
        target.clear();
        target.push_back(some_element);
        trace2("small done", cut_red_bag.size(), cos);
        success = true;;
    }else if (visited_all){
        incomplete();
        // reuse stack
#if IN_PLACE
        savestuff = comp_job.things;
        success = explore_cutsets_in_place(comp_job);
        comp_job.things = savestuff;
#endif
    }else{
        // extend stack, then call
        // work_candidates(comp_job);
        assert( &cut_red_bag == &bagdraft::bag(*this, cut_ext2)); //?
        success = explore_cutsets(cut_red_bag, scratch.cc, some_element, cos, nrs, cut_ext2);
    }

    if(success){
        boost::add_edge(cut_ext2, cut_ext, *this);
        assert(cos);
        return cos;
    }else{
        // free the vertex.
        boost::remove_vertex(cut_ext2, *this);
        // undo visited mask?
        return 0;
    }
} // components.


// only used from q_root. merge!
template <class G2>
template <class CRB, class CRI, class CC, class NeSc>
bool excut_worker<G2>::q_explore_cutsets(
         CRB const & cutred_bag,
         td_vd parent,
         CC& cc_mask,
         CRI celt, unsigned cmps, NeSc* nrs)
{
    typedef T_vertex_descriptor td_vd;

    // assert(cut.size()); no!
    assert(cutred_bag.size() < _bagsize);
    if(cutred_bag.size() + cmps <= _bagsize){ untested();
        incomplete(); // can still happen in q_root_cutset. cleanup...
        // use leaf trick...
        trace2("small", cutred_bag.size(), cmps);
        td_vd leaf=boost::add_vertex(*this);
        BOOST_AUTO(& target, bagdraft::bag(*this, leaf));
        target.push_back(celt);
        boost::add_edge(leaf, parent, *this); // HERE??!
        trace2("small done", cutred_bag.size(), cmps);
        return true;
    }else{
        trace2("ec from q", cutred_bag.size(), celt );
#ifdef EXTRABAG
        /// extra bag at root. prbably a bug.
        td_vd cut_ext=boost::add_vertex(*this);
        auto& cutred_bag2 = bagdraft::bag(*this, cut_ext);
        cutred_bag2=cutred_bag; // one element (root case only here)

//         assert( &cut_ext_bag == &bagdraft::bag(cut_ext, *this)); //?
        bool success=explore_cutsets(cutred_bag2, cc_mask, celt, cmps, nrs, cut_ext);
        if(success){
            boost::add_edge(cut_ext, parent, *this);
            return true;
        }else{
            boost::remove_vertex(cut_ext, *this);
            return false;
        }
#else
        return explore_cutsets(cutred_bag, cc_mask, celt, cmps, nrs, parent);
#endif
    }
}

// reorder component = [cbegin, cend)
// neighbors of cut, then rest.
template<class G_t>
template<class B>
typename excut_worker<G_t>::bfs_range
excut_worker<G_t>::compile_candidates_range(unsigned, B const& cut, const G_t& G,
        typename excut_worker<G_t>::ors* orsp,
        typename excut_worker<G_t>::mask_t& mask)
{
    assert(cut.size());
//    trace2("compile cand", cut.size(), _bagsize);
    return make_bfs_range(cut.begin(), cut.end(), G, &mask, orsp);
//    return OR; // does not work (why?)
}

#ifndef NDEBUG
template<class CB, class C>
void bagfillassert(CB& ext_bag, CB const& cut, C cand)
{
    assert(contains(ext_bag, cand));
    for(auto i: cut){
        assert(contains(ext_bag, i));
    }
}
template<class CB, class C, class E>
void bagfillassert(CB& ext_bag, CB const& cut, C candi, E cande)
{
    for(;candi!=cande;++candi){
        assert(contains(ext_bag, *candi));
    }
    for(auto i: cut){
        assert(contains(ext_bag, i));
    }
}
#endif

template<class G>
bool excut_worker<G>::try_candidate_set(cjob_t& comp_job,
                                        bool neighbours) throw()
{
    auto const& cand_start = comp_job.cut_prefix_end;
//    trace2("try candi", cutbag.size(), *comp_job.cut_prefix_end);

    nrs* nrsp=&comp_job.scratch.nr;
    BOOST_AUTO(&visited, comp_job.visited);
    BOOST_AUTO(const& mask_cc, comp_job.cc_mask); // not initialized yet?! OBSOLETE
    BOOST_AUTO(topp1, comp_job.cut_extip1);
    BOOST_AUTO(&cr_scratch, comp_job.scratch.cr);
    BOOST_AUTO(cut_ext, comp_job.cut_ext);
    auto& cutextbag=bagdraft::bag(*this, cut_ext);
    assert(cutextbag.size()<=_bagsize);

    unsigned all_successful = 1;
    // bagdraft::bag(cut_ext, results).resize(0);
    //      cut_ext= root + candidate, retain ordering!
    assert(cutextbag.size());
    auto const& cand_end=cutextbag.cend();
//    trace2("try_candidate_set set", cand_end-ic, *cand_start);

#ifndef NDEBUG
    auto cutbag = comp_job.cutbag();
    bagfillassert(cutextbag, cutbag, cand_start, cand_end);
#endif

    visited = mask_cc;
    for(auto ic=cand_start; ic!= cand_end; ++ic){
        BOOST_AUTO(pos, get_pos(*ic, _g));
        visited[pos] = true;
    }

//    assert_connected(*comp_job.cut_prefix_end, mask_cc, _g);

    auto const& nvis(visited);

    // compute neighbours of candidates...
    auto N=make_neighbourhood_range(comp_job.cut_prefix_end, cand_end, _g, nvis);
    assert(N.first!=N.second);

    // build components neighboring the candidate set.
    BOOST_AUTO(cmps_range,
            make_components_range(
                N.first, N.second,
                _g, &visited, &cr_scratch, BOOL()));

    // iterate through components
    // break if one fails.
    assert( topp1 <= this->end() );
    BOOST_AUTO(&scratch, comp_job.scratch);
    unsigned cmps_sum=0;
    unsigned cmps_num=0;
    unsigned compno=0;
#ifdef COUNTERS
    auto number_of_candidates=cand_end-cand_start;
#ifndef NDEBUG
    assert(count_unmasked(mask_cc) == count_unmasked(visited) + number_of_candidates);
#endif
    unsigned number_unvisited = comp_job.cc_size - number_of_candidates;
#endif

    // big components first? (how?!)
    for(; cmps_range.first != cmps_range.second; ++cmps_range.first){
        ++compno;

        assert(cmps_range.first != cmps_range.second);

        BOOST_AUTO(comp_range, *(cmps_range.first));
        assert(comp_range.first != comp_range.second);
        all_successful = do_component(cut_ext, comp_range,
                scratch, nrsp, &cut_ext
#ifdef COUNTERS
                , number_unvisited
//                , comp_job // FIXME: *use* this...
                , neighbours
#endif
                );
#ifndef COUNTERS
        (void) neighbours;
#endif
        if(!all_successful) break;

        cmps_sum += all_successful;
        ++cmps_num;
    }

    if(all_successful){
        // boost::add_edge(cut_ext, cut_red, results);
        return true;
    }else{
#if 0 // later?
        for(auto ic=cand_start; ic!= cand_end; ++ic){
            BOOST_AUTO(pos, get_pos(*ic, _g));
        }
#endif
        // delete anything the above do_component calls have created.
        // ought to be all children of cut_ext
        // how to do this with boost::??
        this->erase(topp1, this->end()); // keep cut_ext
        /// --> other threads must clean up. how?!
        // try next candidate
        return false;
    }
}

template<class G>
bool excut_worker<G>::work_candidates(cjob_t& comp_job)
{
    bool success=false;

    assert(comp_job._candidates_b!=comp_job._candidates_e);
    BOOST_AUTO(&candi, comp_job._candidates_b);
    BOOST_AUTO(&cande, comp_job._candidates_e);
    assert(candi!=cande);

    auto const& cutbag(comp_job.cutbag());
    auto& candidate_cache=comp_job.cand_cache;
    // BOOST_AUTO(cut_exti, comp_job.cut_exti);
//    BOOST_AUTO(cut_ext, *cut_exti);
    BOOST_AUTO(cut_ext, comp_job.cut_ext);
    unsigned number_of_candidates = _bagsize - cutbag.size();
    // unsigned number_of_candidates = 1;
    assert(number_of_candidates);
    candidate_cache.resize(0);

    auto& cutextbag=bagdraft::bag(*this, cut_ext);

    cutextbag = cutbag; // uuh!!
    auto prefix_size=cutextbag.size();
    auto firstcandpos=cutextbag.end();
    auto seek=cutextbag.end();
    comp_job.cut_prefix_end=cutextbag.end();

    unsigned done=0;
    // prepare a prefix
    while(done<number_of_candidates-1 && candi != cande){
        candidate_cache.push_back(*candi);
        cutextbag.push_back(*candi);
        ++done;
        ++candi;
    }
    assert(cutextbag.size()<=_bagsize);

#ifndef NDEBUG
        bagfillassert(cutextbag, cutbag, candidate_cache.begin(), candidate_cache.end());
#endif
    bool neighbours=true;

    // more candidates?
    cutextbag.resize(prefix_size+number_of_candidates);
    assert(cutextbag.size()<=_bagsize);
    for(;candi != cande; ++candi){
        *firstcandpos=*candi;
        // subsets without newly added candidate.
        //
        auto subsets=make_subsets_range(candidate_cache.begin(),
                                        candidate_cache.end(),
                                        number_of_candidates-1,
                                        number_of_candidates-1);
        auto Si = subsets.first;
        auto Se = subsets.second;
        for(; Si!=Se; ++Si){
            // iterate subsets...
            auto subset=*Si;

            seek=firstcandpos; ++seek;
            for(; subset.first!=subset.second; ++(subset.first) ){
                assert(number_of_candidates>1);
                *seek=*(subset.first);
                ++seek;
            }
            assert(seek==cutextbag.end());

            if(number_of_candidates==1){
                assert(cutextbag.back() == *candi);
            }
            success = try_candidate_set(comp_job, neighbours);
            if(success){
                return true;
            }

            // keep prefix. next subset...
            // cutextbag.resize(prefix_size+1); //??
        }
        // use in next loop.
        candidate_cache.push_back(*candi);
#if COUNTERS
        if(neighbours){
            neighbours = candi.is_neighbour();
        }else{
        }
#endif
    }

    return success;
}

// recursively explore cutsets
// try candidates until all components are successful
template<class G>
template<class CRB, class CRI>
bool excut_worker<G>::explore_cutsets(CRB const& cut_ext_bag,
         mask_t& cc_mask_in, CRI /*celt*/, unsigned cmps, nrs* nrsa, td_vd cut_ext)
{
    assert( &cut_ext_bag == &bagdraft::bag(*this, cut_ext)); // no. root?
//    assert(count_unmasked(cc_mask_in) == cmps);
//    assert(count_range(cbegin, cend) == cmps);
    (void) cmps;
    assert(cut_ext_bag.size() + cmps > _bagsize);

    BOOST_AUTO(topp1, boost::vertices(*this).second);
    assert(topp1 == this->end());

    // rest = component \without N
    //      = component \without ( neighbors(cut_red_bag) \intersect component )
    //      = component \without neighbours(cut_red_bag)
    //
    cjob_t &comp_job = *new_cj(_bagsize, cut_ext_bag);

    std::vector<EXCUT_BOOL>& visited(comp_job.visited);
    std::vector<EXCUT_BOOL>& mask_cc(comp_job.cc_mask);
    assert(visited.size() == boost::num_vertices(_g));

    visited.assign(boost::num_vertices(_g), true);
    // TODO: dummy compile_candidates...
    mask_cc = cc_mask_in;
#ifdef COUNTERS
    comp_job.cc_size = cmps;
#endif
    // visited_cc=cc_mask_in;
    auto P = compile_candidates_range(cmps, cut_ext_bag, _g, nrsa, cc_mask_in);
    assert(P.first!=P.second);
    comp_job._candidates_b = MOVE(P.first);
    comp_job._candidates_e = MOVE(P.second);

    // find a candidate...
    bool success = false;

//    comp_job.cutbag = &cut_ext_bag;

    comp_job.cut_ext = cut_ext;
    comp_job.cut_extip1 = topp1;
//    assert(*top == cut_ext);

//    assert_connected(*comp_job._candidates_b, mask_cc, _g);

    success = work_candidates(comp_job);
    recycle(&comp_job);

    return success;
}

} //detail

} // treedec
namespace treedec{

// FIXME: proper types. use graph.hpp etc.
template<class G>
std::pair<unsigned,unsigned> find_max_degree_vertex(G const& g)
{
    assert(boost::num_vertices(g));
    auto V=boost::vertices(g);
    auto cand=*V.first;
    unsigned degree=boost::degree(*V.first, g);
    ++V.first;

    for(; V.first!=V.second; ++V.first){
        auto degree_thisone=boost::degree(*V.first, g);
        if(degree_thisone > degree){
            cand=*V.first;
            degree = degree_thisone;
        }
    }
    return std::make_pair(cand, degree);
}

// compute a tree decomposition of bagsize bs, if it exists
namespace draft{

template <typename G_t,
        template<class G_> class config=algo::default_config>
class exact_cutset { // baseclass?
public:
    exact_cutset(G_t const& g)
        : _g(g) {}
    ~exact_cutset() {
    }
public:
    template<class T_t>
    bool try_it(T_t &T, unsigned bs);
    template<class T_t>
    void do_it(T_t &T, unsigned& bs){
        while(!try_it(T, bs)){
            bs++;
        }
    }
private:
    G_t const& _g;
};

template <typename G_t, template<class G_> class config>
template<class T_t>
bool exact_cutset<G_t, config>::try_it(T_t &T, unsigned bagsize)
{
    trace1("exact_cutset", bagsize);
    assert_connected(_g);

    typedef typename boost::graph_traits<T_t>::vertex_descriptor vertex_descriptor_T;
//    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    if(boost::num_vertices(_g) == 0){
        boost::add_vertex(T);
        return true;
    }else{
        incomplete();
    }

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(_g);

    if(boost::num_vertices(_g) == 1){
        vertex_descriptor_T t=boost::add_vertex(T);
        insert(bag(t, T), *vIt);
        if(bagsize <= 1){
            return true;
        }
        return false;
    }

    treedec::detail::excut_control<G_t> c(_g, bagsize);

    VECTOR_TD<G_t>& results=c._results;

    trace0("init newone");
    if(bagsize<=1){
        return false;
    }

    typedef typename treedec::detail::excut_worker<G_t>::T_vertex_descriptor T_vertex_descriptor;

    assert(boost::num_vertices(results)==0);
    T_vertex_descriptor root=boost::add_vertex(results);
    bagdraft::bag(results, root).reserve(bagsize);

    auto M = find_max_degree_vertex(_g);
    bagdraft::bag(results, root).push_back(M.first);

    typedef typename treedec::detail::excut_control<G_t>::bag_type B;
    typedef typename B::const_iterator bag_const_iterator;
    auto nrsp=new_bfs_range_scratch(_g,
                std::pair<bag_const_iterator, bag_const_iterator>(), EXCUT_BOOL());

    unsigned nv=boost::num_vertices(_g);

    std::vector<EXCUT_BOOL> ccm(nv, false);
    ccm[M.first] = true;

    assert(ccm.size() == boost::num_vertices(_g));
    c.q_root_cutset(root, ccm, M.first, nv-1, nrsp);
    c.run();
    bool success=c.join();
    delete nrsp;

    if(success){
    }else{
        results.erase(results.begin(), results.end());
    }

    if(success){
    }else{
        return false;
    }

    trace2("done excut. write back ", boost::num_vertices(T), &_g);
    // TODO: just use results as tree.
    assert(!boost::num_vertices(T));

    for(BOOST_AUTO(ii, results.begin()); ii!=results.end(); ++ii){
        boost::add_vertex(T);
    }
    unsigned bag_index=0;
    std::vector<BOOL> visited(boost::num_vertices(_g));
    typedef typename boost::graph_traits<G_t>::vertex_iterator vit_G;
    typedef std::pair<vit_G, vit_G> VRP;
    typename ::detail::components_iter<G_t, VRP, EXCUT_BOOL>::scratch_type crscr;
    visited.assign(boost::num_vertices(_g), false);

    auto rv=boost::vertices(results);

    for(BOOST_AUTO(ii, rv.first); ii!=rv.second; ++ii){
        auto& source_bag=ii->second;

        auto parent_it=boost::adjacent_vertices(*ii, results).first;
        unsigned parent_pos = get_pos(*parent_it, results);
        if(boost::degree(*ii,results)){
            assert(boost::degree(*ii,results) == 1);
            assert(bag_index);
            assert(bag_index<boost::num_vertices(T));
            boost::add_edge(bag_index, parent_pos, T);
        }else{
            assert(!bag_index);
            // root...
        }
        auto& target_bag=bag(get_pos(*ii, results), T);
///// =============  LEAFTRICK =========
        if(bag_index && source_bag.size()==1){
            auto const& parentB=bag(get_pos(*parent_it, results), T);
           // trace2("leaftrick", bag_index, *target_bag.begin());

            // mark parent bag visited.
            for( auto jj : parentB){
//                trace1("", jj);
        //        target_bag.insert(jj);
                auto pos=get_pos(jj, _g);
                visited[pos]=true;
            }
            // put connected component of leaf vertex into final bag
            auto N=make_components_range(source_bag.begin(), source_bag.end(),
                    _g, &visited, &crscr, BOOL());

            auto C=*(N.first);
            for(;C.first!=C.second; ++(C.first)){
                auto lbv=*C.first;
                auto pos=get_pos(lbv, _g);
                (void)pos;
                // assert(visited[pos]);
                insert(target_bag, lbv); //pos?
            }
            auto S=target_bag;
            for(auto ti=S.begin(); ti!=S.end(); ++ti){
                auto lN = boost::adjacent_vertices(*ti, _g);
                for(; lN.first!=lN.second; ++(lN.first)){
                    if(contains(parentB, *(lN.first))){
                        insert(target_bag, *(lN.first));
                    }
                }
            }
        }else{
            // no leaftrick. move??
            for(BOOST_AUTO(jj, source_bag.begin()); jj!=source_bag.end(); ++jj){
                insert(target_bag, *jj); //pos?
            }
        }

        ++bag_index;
    }

    // bagsize = treedec::get_bagsize(T);
    // boost::put(T, treedec::bagsize, bagsize);
    return true;
} // do_it

} // draft

template <typename G_t, typename T_t>
bool exact_cutset(G_t const &G, T_t &T, int tw)
{
    draft::exact_cutset<G_t> a(G);
    unsigned bs=tw+1;
    return a.try_it(T, bs);
}

// compute a tree decomposition
// return the bagsize.
template <typename G_t, typename T_t>
unsigned exact_cutset(G_t &G, T_t &T){ untested();
    int lb_bs = 0;

    draft::exact_cutset<G_t> a;
    while(!a.do_it(G, T, lb_bs)){
        ++lb_bs;
    }
    return lb_bs;
}

} //namespace treedec

#endif //TD_EXACT_CUTSET
// vim:ts=8:sw=4:et
