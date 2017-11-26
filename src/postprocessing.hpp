// Lukas Larisch, 2014 - 2017
// Felix Salfelder 2016
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
 * Offers functionality to possibly reduce the width of a tree decomposition of a given graph.
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void MSVS(G_t &G, T_t &T)
 * - void minimalChordal(G_t G, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &old_elimination_ordering,
 *                              std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &new_elimination_ordering)
 *
 */

#ifndef TREEDEC_POSTPROCESSING_HPP
#define TREEDEC_POSTPROCESSING_HPP

#include <boost/graph/adjacency_list.hpp>
#include "elimination_orderings.hpp"
#include "network_flow.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "graph.hpp"
#include "overlay.hpp"
#include "treedec.hpp"

#define get_pos(a,b) ( boost::get(boost::vertex_index, b, a) )

namespace treedec{

//Check a modified induced subgraph of the bag 'bag(t_desc, T)' for
//possible improvement.
template <typename B_t, typename S_t, typename vd_t>
bool is_improvement_bag(B_t const &H,
                        std::vector<BOOL> &disabled,
                        S_t &X,
                        S_t &Y,
                        vd_t a,
                        vd_t b)
{

    typedef typename boost::graph_traits<B_t> graph_traits;
    typedef typename graph_traits::adjacency_iterator adjacency_iterator_H;

    disabled.assign(boost::num_vertices(H), false);

    //Test for completeness.
    typename graph_traits::vertices_size_type nv=boost::num_vertices(H);
    typename graph_traits::edges_size_type ne=boost::num_edges(H);
    if(nv*(nv-1u) == 2*ne){
        //There is no remaining non-edge.
        X.clear();
        Y.clear();
        return false;
    }

    assert(a!=b);
    assert(!boost::edge(a, b, H).second);

    //Collect the neighbours of x and y, resulting in the sets X and Y.
    adjacency_iterator_H nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(a, H); nIt!=nEnd; ++nIt){
        X.push_back(*nIt);
    }

    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(b, H); nIt!=nEnd; ++nIt){
        Y.push_back(*nIt);
    }

    unsigned int pos1 = get_pos(a, H);
    unsigned int pos2 = get_pos(b, H);
    disabled[pos1] = true;
    disabled[pos2] = true;

    return true;
}

/* MinimalSeparatingVertexSet(MSVS)-algorithm
 *
 * Try to find a minimal separator S in a graph H, that
 *    (1) contains the induced subgraph of a maximum-sized bag B(t) of T,
 *    (2) has an edge {x,y}, if {x,y} is a subset of B(t') for some neighbour t' of t in T.
 * If no separator can be found for none of the maximum-sized bags, the algorithm stops. Otherwise,
 * the tree decomposition T is refined according to S.
 *
 */
namespace impl{

template <typename G_t, typename T_t>
class MSVS{
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<T_t>::vertex_descriptor bag_descriptor;
    typedef typename boost::graph_traits<T_t>::vertex_iterator bag_iterator;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor G_vertex_descriptor;
    typedef typename graph_traits<G_t>::immutable_type immutable_type;
    typedef typename boost::graph_traits<immutable_type>::vertex_descriptor imm_vertex_descriptor;

public:
    MSVS(G_t const &G, T_t &T)
        : _g(G), _t(T)
    {
        assert(is_valid_treedecomposition(_g, _t));
    }

    void do_it(){
        std::vector<BOOL> disabled;
        std::vector<BOOL> disabled_;
        vertices_size_type bagsize = treedec::get_bagsize(_t);
        std::set<vertex_descriptor> S;
        std::vector<imm_vertex_descriptor> X, Y;

        immutable_type H; // malloc/free, where?

        std::set<typename immutable_type::vertex_descriptor> S_;
        std::vector<G_vertex_descriptor> vdMap_, vdMap;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component;

        //TODO: (terribly) inefficient:
        std::vector<std::set<vertex_descriptor> > union_S_W_i;

        std::vector<bag_descriptor> newN;
        std::vector<bag_descriptor> oldN;
        typename treedec_traits<T_t>::bag_type intersection;

        while(true){
            bagsize = treedec::get_bagsize(_t);

            callback(bagsize);

            //Check all maximum sized bags, whether they can be improved or not. Take the first improvable.
            H.clear();
            X.clear();
            Y.clear();
            disabled.resize(0);
            disabled_.resize(0);
            vdMap.resize(0);

            bag_iterator tIt, tEnd;
            bag_descriptor refinement_vertex;
            immutable_type const* HI=NULL;
            bool status=false;

            for(boost::tie(tIt, tEnd) = boost::vertices(_t); tIt!=tEnd; ++tIt){
                callback(0);
                if(bag(*tIt, _t).size() == bagsize){
                    disabled_.resize(0);
                    vdMap_.resize(0);

                    /* draft:: */ is_in_neighbour_bd<vertex_descriptor, T_t> cb(_t, *tIt);
                    BOOST_AUTO(mybag, bag(*tIt, _t));
                    HI = &treedec::draft::immutable_clone(_g, H, mybag.begin(), mybag.end(), mybag.size(), &vdMap_, &cb);
                    status = is_improvement_bag<immutable_type, std::vector<imm_vertex_descriptor>, long unsigned>
                                                                                  (*HI, disabled_, X, Y, cb.a, cb.b);

                    if(status){
                        refinement_vertex = *tIt;
                        disabled = MOVE(disabled_);
                        vdMap = MOVE(vdMap_);
                        assert(vdMap.size() == boost::num_vertices(*HI));
                        break;
                    }
                }
            }

            if(!status){
                return;
            }
            assert(HI);

#ifndef NDEBUG
            std::vector<BOOL>::const_iterator x=disabled.begin();
            unsigned num_dis=0;
            for(; x!=disabled.end(); ++x){
                if(*x) ++num_dis;
            }
            assert(num_dis==2);
#endif

            //Compute a seperating set S.
            S_.clear();
            seperate_vertices(H, disabled, 2, X, Y, S_);

            //S consists of vertex descriptors of H. Use vd_map to map these to descriptors of _g.
            S.clear();
            map_descriptors(S_, S, *HI, vdMap);

            //Mark the vertices of the seperator as visited (used for computing connected components).
            std::vector<BOOL> visited(boost::num_vertices(H), false);
            BOOST_AUTO(sIt, S_.begin());
            for(; sIt!=S_.end(); ++sIt){
                auto pos=boost::get(boost::vertex_index, *HI, *sIt);
                visited[pos] = true;
            }

            //Convert the descriptors stored in S to a bag and replace
            //the bag of 'refinement_vertex' with this bag.
            typename treedec_traits<T_t>::bag_type B;
            treedec::map_descriptors_to_bags<G_t>(S, B);
            typename treedec_traits<T_t>::bag_type old_bag = bag(refinement_vertex, _t);
            bag(refinement_vertex, _t) = MOVE(B);

            //Store the connected components of H[V(H)\S] in 'components'.
            typedef typename boost::graph_traits<immutable_type>::vertex_descriptor HI_vertex_descriptor;
            std::vector<std::set<HI_vertex_descriptor> > components;
            treedec::get_components_provided_map(*HI, components, visited);

            //Store the (neighbours of 'refinement_vertex' in _t) in 'oldN'.
            typename boost::graph_traits<T_t>::adjacency_iterator t_nIt, t_nEnd;
            oldN.resize(boost::out_degree(refinement_vertex, _t));
            unsigned int c = 0;
            for(boost::tie(t_nIt, t_nEnd) = boost::adjacent_vertices(refinement_vertex, _t); t_nIt != t_nEnd; t_nIt++){
                oldN[c++] = *t_nIt;
            }

            boost::clear_vertex(refinement_vertex, _t);

            //'refinement_vertex' gets |connected_components|-many neighbours, that are new vertices in _t.
            //The bag of neighbours i will be (connected_components[i] v S).
            union_S_W_i.resize(components.size());

            //TODO: proper container!
            BOOST_AUTO(ii, union_S_W_i.begin());
            for(; ii != union_S_W_i.end(); ++ii ){
                ii->clear();
            }

            newN.resize(components.size());

            for(unsigned int i = 0; i < components.size(); i++){
                component.clear();
                treedec::map_descriptors(components[i], component, *HI, vdMap);

                std::set_union(S.begin(), S.end(),
                               component.begin(), component.end(),
                               std::inserter(union_S_W_i[i], union_S_W_i[i].begin()));

                newN[i] = boost::add_vertex(_t);
                typename treedec_traits<T_t>::bag_type uB;
                treedec::map_descriptors_to_bags<G_t>(union_S_W_i[i], uB);
                bag(newN[i], _t) = MOVE(uB);

                assert(!boost::edge(refinement_vertex, newN[i], _t).second);
                assert(!boost::edge(newN[i], refinement_vertex, _t).second);

                treedec::add_edge(refinement_vertex, newN[i], _t);

                 assert(boost::edge(refinement_vertex, newN[i], _t).second);
                assert(boost::edge(newN[i], refinement_vertex, _t).second);
            }

            //Let intersection_i be the intersection of the old bag of 'refinement_vertex' with
            //the bag of the old neighbour i. Add an edge between i and a new neighbour j of
            //'refinement_vertex', if the bag of j includes intersection_i. This can only be done
            //with exactly one new neighbour of 'refinement_vertex'.
            for(unsigned int i = 0; i <  oldN.size(); i++){
                intersection.clear();
                std::set_intersection(old_bag.begin(), old_bag.end(),
                                      bag(oldN[i], _t).begin(),
                                      bag(oldN[i], _t).end(),
                                      std::inserter(intersection, intersection.begin()));

                for(unsigned int j = 0; j < newN.size(); j++){
                    if(std::includes(bag(newN[j], _t).begin(), bag(newN[j], _t).end(),
                                     intersection.begin(), intersection.end()))
                    {
                        boost::add_edge(newN[j], oldN[i], _t);
                        break;
                    }
                }
            }
        }
    }
    protected:
        virtual void callback(vertices_size_type){}
    private:
        G_t const& _g;
        T_t& _t;
}; //MSVS

} // impl

template <typename G_t, typename T_t>
void MSVS(G_t const &G, T_t &T)
{
    impl::MSVS<G_t, T_t> A(G, T);
    A.do_it();
}

template <typename G_t, class O_t, class E_t>
bool is_candidate_edge(E_t &edge, unsigned int i,
                       O_t &elimination_ordering, G_t &M)
{
    //Position i in 'elimination_ordering_' will store the 'elimination date' of vertex i
    std::vector<unsigned int> elimination_ordering_(boost::num_vertices(M));
    for(unsigned int t = 0; t < elimination_ordering.size(); t++){
        auto pos=boost::get(boost::vertex_index, M, elimination_ordering[t]);
        elimination_ordering_[pos] = t;
    }

    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(edge.first, M); nIt != nEnd; nIt++){
        unsigned int pos = get_pos(*nIt, M);
        if(elimination_ordering_[pos] > i && boost::edge(edge.second, *nIt, M).second
       && !boost::edge(*nIt, elimination_ordering[i], M).second)
        {
            return false;
        }
    }

    return true;
}

template <typename G_t, typename E_t>
inline void delete_edges(G_t &G, E_t &edges){
    for(unsigned int i = 0; i < edges.size(); i++){
        boost::remove_edge(edges[i].first, edges[i].second, G);
    }
}

namespace impl {

template <typename G_t, typename O_t, template<class GG, class ...> class CFGT>
class minimalChordal{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef std::pair<vertex_descriptor, vertex_descriptor> edge_type;
public:
    minimalChordal(G_t& g, O_t& o)
        : _g(g), _o(o)
    {}

public:
    void do_it();
    O_t const& ordering() const {return _no;}

private:
    G_t& _g;
    O_t const& _o;
    O_t _no;
};

template <typename G_t, typename O_t, template<class GG, class ...> class CFGT>
inline void impl::minimalChordal<G_t, O_t, CFGT>::do_it()
{
    _no.resize(_o.size());
    //Make 'G' a filled-in graph according to '_o'. This operation stores
    //all new edges in F.
    std::vector<std::set<vertex_descriptor> > C;
    std::vector<std::vector<std::pair<vertex_descriptor, vertex_descriptor> > > F;
    treedec::make_filled_graph(_g, _o, C, F);

    for(int i = _o.size()-1; i >= 0; i--){
        //Checks if F[i][j] is an candidate edge. If this is the case, F[i][j] will be stored in
        //'candidate'. The endpoints of F[i][j] will be stored in 'incident'.
        std::vector<edge_type> candidate;
        std::set<vertex_descriptor> incident;
        for(unsigned int j = 0; j < F[i].size(); j++){
            if(treedec::is_candidate_edge(F[i][j], i, _o, _g)){
                candidate.push_back(F[i][j]);
                incident.insert(F[i][j].first);
                incident.insert(F[i][j].second);
            }
        }
        if(candidate.size() != 0){
            //If there is some candidate edge, create W_i := (incident, G[incident] \ candidate)
            //and run the LEX_M algorithm. The algorithm will possibly return some edges not in W_I that
            //have to be added to W_i to make the graph chordal. 
            G_t W_i;
            typename std::vector<vertex_descriptor> vdMap;
            treedec::induced_subgraph_omit_edges(W_i, _g, incident, candidate, vdMap);

            std::vector<edge_type> keep_fill_;
            treedec::LEX_M_fill_in(W_i, keep_fill_);

            //Translate descriptors of W_i to descriptors of G.
            std::vector<edge_type> keep_fill(keep_fill_.size());
            for(unsigned int j = 0; j < keep_fill_.size(); j++){
                auto pos1=boost::get(boost::vertex_index, W_i, keep_fill_[j].first);
                auto pos2=boost::get(boost::vertex_index, W_i, keep_fill_[j].second);
                keep_fill[j].first = vdMap[pos1];
                keep_fill[j].second = vdMap[pos2];
            }

            //Delete all candidate edges that can be deleted in G according to LEX_M_fill_in.
            for(unsigned int j = 0; j < candidate.size(); j++){
                for(unsigned int k = 0; k < keep_fill.size(); k++){
                    if(  (candidate[j].first == keep_fill[k].first
                       && candidate[j].second == keep_fill[k].second)
                       ||(candidate[j].first == keep_fill[k].second
                       && candidate[j].second == keep_fill[k].first))
                    {
                        candidate.erase(candidate.begin()+j);
                        break;
                    }
                }
            }

            treedec::delete_edges(_g, candidate);
        }
    }
    treedec::LEX_M_minimal_ordering(_g, _no);
}

} // impl

/* minimalChordal-algorithm
 *
 * Compute possibly redundant fill-in-edges and runs LEX-M to check,
 * if the graph after removal of a fill-in-edge is chordal.
 * Finally, compute a perfect elimination ordering,
 * possibly of lower width than '_o'.
 */
template <typename G_t, class O_t>
inline void minimalChordal(G_t &G,
     O_t& old_elimination_ordering,
     O_t& new_elimination_ordering)
{
    ::treedec::impl::minimalChordal<G_t, O_t, algo::default_config> A(G, old_elimination_ordering);
    A.do_it();
    new_elimination_ordering = A.ordering();
}

} // treedec

#undef get_pos

#endif //ifdef TREEDEC_POSTPROCESSING_HPP

// vim:ts=8:sw=4:et:
