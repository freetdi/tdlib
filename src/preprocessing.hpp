// Lukas Larisch, 2014 - 2016
// Felix Salfelder 2016
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


/*
 * Offers functionality to preprocess a graph, such that after
 * the repeated application of reduction rules, which in case the
 * input graph has tree-width at most 3 allow us to determine it's tree-width exactly
 * and in addition compute the corresponding tree decomposition. If the tree-width
 * is larger, the reduction rules return a possibly smaller instance of the same
 * tree-width as the original graph, a 'partial' tree decomposition and a lower bound
 * with respect to tree-width, such that
 * further algorithms can be applied to the resulting graph.
 *
 *
 * These functions are most likely to be interesting for outside use:
 *
 * -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags)
 * -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &lb)
 *
 */

#ifndef TD_PREPROCESSING
#define TD_PREPROCESSING

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>

#include "graph.hpp"
#include "misc.hpp"

namespace treedec{

namespace impl{

// Check if there exists a degree-0-vertex.
template <typename G_t, typename B_t>
void Islet(G_t &G, B_t &bags, int &low)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            typename treedec_traits<T_t>::vd_type vd=get_vd(G, *vIt);
            typename treedec_traits<T_t>::bag_type emptybag;

            bags.push_back(boost::make_tuple(vd, MOVE(emptybag)));

            low = (low > 0)? low : 0;
        }
    }
}

template <typename G_t, typename T_t>
void Islet(G_t &G, T_t &bags)
{
    int low = -1;
    Islet(G, bags, low);
}

/* (Islet,) Twig and Series rules. */
template <typename G_t, typename B_t, typename DEGS>
void eliminate_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G,
         B_t &bags, int &low, DEGS &degs)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typename treedec_traits<T_t>::bag_type bag;

    unlink_1_neighbourhood(v, G, degs);
    degs.unlink(v);

    unsigned deg = boost::degree(v, G);
    typename outedge_set<G_t>::type xbag;
    treedec::make_clique_and_detach(v, G, xbag);
    assert(!boost::degree(v, G));

    redegree(NULL, G, xbag, degs);

    bags.push_back(boost::make_tuple(v, MOVE(xbag)));

    low = (low > (int)deg)? low : deg;
}

//Apply the Triangle rule if applicable (checks if there exists a
//degree-3-vertex, such that at least one edge exists in its neighbourhood).
//return true, if degs has been modified.
template <typename G_t, typename T_t, typename DEGS>
bool Triangle(G_t &G,
              typename boost::graph_traits<G_t>::vertex_descriptor v,
              T_t &bags, int &low, DEGS &degs)
{
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;

    std::vector<vertex_descriptor> N(3);
    adjacency_iterator f=boost::adjacent_vertices(v, G).first;
    N[0] = *f;
    N[1] = *(++f);
    N[2] = *(++f);

    if(boost::edge(N[0], N[1], G).second
     || boost::edge(N[0], N[2], G).second
     || boost::edge(N[1], N[2], G).second)
    {
        unlink_1_neighbourhood(v, G, degs);
        degs.unlink(v, 3);
        typename outedge_set<G_t>::type xbag;
        treedec::make_clique_and_detach(v, G, xbag);
        redegree(NULL, G, xbag, degs);

        bags.push_back( boost::make_tuple(v, MOVE(xbag)));

        low = (low > 3)? low : 3;
        return true;
    }

    return false;
}

//Applies the Buddy rule if applicable (checks if there exists two degree-3-vertices,
//such that they share their neighbours)
template <typename G_t, typename DEGS>
bool Buddy(G_t &G,
           typename boost::graph_traits<G_t>::vertex_descriptor v,
           typename boost::graph_traits<G_t>::vertex_descriptor w,
           std::vector<boost::tuple<
             typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
             typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
           > > &bags, int &low, DEGS &degs)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typedef typename treedec_traits<T_t>::vd_type vd_type;

    if(check_twins(v, w, G)){
        typename treedec_traits<T_t>::bag_type N1, N2;
        assign_neighbours(N1, v, G);
        unlink_1_neighbourhood(v, G, degs);
        degs.unlink(v, 3);
        degs.unlink(w, 3);
        typename outedge_set<G_t>::type xbag;
        treedec::make_clique_and_detach(v, G, xbag);
        boost::clear_vertex(w, G);
        redegree(NULL, G, xbag, degs);

        vd_type vd1 = get_vd(G, v);
        vd_type vd2 = get_vd(G, w);

        bags.push_back(boost::make_tuple(vd1, xbag));
        bags.push_back(boost::make_tuple(vd2, MOVE(xbag)));

        low = (low > 3)? low : 3;
        return true;
    }
    return false;
}

template<class T>
inline void rearrange_cube(T* N, T x)
{
    if(N[0] == x){
        N[0] = N[2];
    }
    else if(N[1] == x){
        N[1] = N[2];
    }
    else{}
}

//Applies the Cube rule if applicable.
template <typename G_t, typename DEGS>
bool Cube(G_t &G,
          typename boost::graph_traits<G_t>::vertex_descriptor vertex,
          std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
            typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
           > > &bags, int &low, DEGS &degs)
{
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename treedec_chooser<G_t>::type T_t;
    typedef typename treedec_traits<T_t>::vd_type vd_type;

    vertex_descriptor x, a, b, c;
    x = vertex;

    adjacency_iterator nIt, nEnd;
    adjacency_iterator f=boost::adjacent_vertices(x, G).first;
    a = *f;
    b = *(++f);
    c = *(++f);

    if(boost::degree(a, G)!=3){
        return false;
    }else if(boost::degree(b, G)!=3){
        return false;
    }else if(boost::degree(c, G)!=3){
        return false;
    }

    vertex_descriptor N[9];
    vertex_descriptor* Na=&N[0];
    vertex_descriptor* Nb=&N[3];
    vertex_descriptor* Nc=&N[6];

    f = boost::adjacent_vertices(a, G).first;
    Na[0] = *f; Na[1] = *(++f); Na[2] = *(++f);
    f = boost::adjacent_vertices(b, G).first;
    Nb[0] = *f; Nb[1] = *(++f); Nb[2] = *(++f);
    f = boost::adjacent_vertices(c, G).first;
    Nc[0] = *f; Nc[1] = *(++f); Nc[2] = *(++f);

    rearrange_cube(Na, x);
    rearrange_cube(Nb, x);
    rearrange_cube(Nc, x);

    typename treedec_traits<T_t>::bag_type bag;
    assign_neighbours(bag, a, b, c, G);

    if(bag.size() != 4){
        return false;
    }

    typename boost::graph_traits<G_t>::vertex_descriptor u, v, w;

    if(Na[0] == Nb[0]){      u = Na[0]; v = Na[1]; w = Nb[1]; }
    else if(Na[0] == Nb[1]){ u = Na[0]; v = Na[1]; w = Nb[0]; }
    else if(Na[1] == Nb[0]){ u = Na[1]; v = Na[0]; w = Nb[1]; }
    else if(Na[1] == Nb[1]){ u = Na[1]; v = Na[0]; w = Nb[0]; }
    else{ return false; }

    if((Nc[0] == v && Nc[1] == w) || (Nc[1] == v && Nc[0] == w)){
        bag.clear();
        vd_type vdx = get_vd(G, x);
        vd_type vda = get_vd(G, a);
        vd_type vdb = get_vd(G, b);
        vd_type vdc = get_vd(G, c);
        vd_type vdu = get_vd(G, u);
        vd_type vdv = get_vd(G, v);
        vd_type vdw = get_vd(G, w);

        bag.insert(vdu); bag.insert(vdv); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vda, MOVE(bag)));
        bag.clear();

        bag.insert(vdw); bag.insert(vdv); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vdc, MOVE(bag)));
        bag.clear();

        bag.insert(vdw); bag.insert(vdu); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vdb, MOVE(bag)));

        degs.unlink(a, 3);
        degs.unlink(b, 3);
        degs.unlink(c, 3);
        degs.unlink(u);
        degs.unlink(v);
        degs.unlink(w);

        boost::clear_vertex(a, G);
        boost::clear_vertex(b, G);
        boost::clear_vertex(c, G);

        boost::add_edge(u, v, G);
        boost::add_edge(u, w, G);
        boost::add_edge(u, x, G);
        boost::add_edge(v, w, G);
        boost::add_edge(v, x, G);
        boost::add_edge(w, x, G);

        degs.reg(u);
        degs.reg(v);
        degs.reg(w);

        low = (low > 3)? low : 3;
        return true;
    }
    return false;
}

//Applies the Simplicial rule, if possible (checks if there exists a vertex,
//such that its neighbours induce a clique).
template <typename G_t, typename DEGS>
bool Simplicial(G_t &G,
                typename boost::graph_traits<G_t>::vertex_descriptor v,
                std::vector<boost::tuple<
                  typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
                  typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
                > > &bags, int &low, DEGS &degs)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typedef typename treedec_traits<T_t>::vd_type vd_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    //The neighbourhood of v is a clique, if no "edge miss" occures.
    bool isClique = true;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(v, G); nIt1 != nEnd; nIt1++){
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            if(!boost::edge(*nIt1, *nIt2, G).second){
                isClique = false;
                goto DOUBLE_BREAK;
            }
        }
    }

    DOUBLE_BREAK:

    if(isClique){
        vd_type vd = get_vd(G, v);

        unlink_1_neighbourhood(v, G, degs);
        degs.unlink(v);
        bag_type xbag;
        treedec::make_clique_and_detach(v, G, xbag);
        redegree(NULL, G, xbag, degs);

        if (unsigned(low) < xbag.size()){ untested();
            low = xbag.size();
        }else{ untested();
        }

        bags.push_back(boost::make_tuple(vd, MOVE(xbag)));
        return true;
    }
    return false;
}


//Applies the Almost Simplicial rule if possible.
template <typename G_t, typename DEGS>
bool AlmostSimplicial(G_t &G,
                      typename boost::graph_traits<G_t>::vertex_descriptor v,
                      std::vector<boost::tuple<
                        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
                        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
                      > > &bags, int &low, DEGS &degs)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename treedec_traits<T_t>::vd_type vd_type;

    bool isAlmostSimplicial = true;
    bool specialNeighbourFound = false;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    vertex_descriptor cand1, cand2, specialNeighbour;
    unsigned int missingEdgesCount;

    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(v, G); nIt1 != nEnd; ++nIt1){
        nIt2 = nIt1;
        nIt2++;
        missingEdgesCount = 0;
        for(; nIt2 != nEnd; nIt2++){
            if(specialNeighbourFound && (*nIt1 == specialNeighbour || *nIt2 == specialNeighbour)){
                continue;
            }

            if(!boost::edge(*nIt1, *nIt2, G).second){
                if(specialNeighbourFound){
                    //#special neighbours > 1.
                    isAlmostSimplicial = false;
                    goto DOUBLE_BREAK;
                }
                //*nIt1 or *nIt2 is a special neighbour.
                cand1 = *nIt1;
                cand2 = *nIt2;
                missingEdgesCount++;
            }
        }

        if(missingEdgesCount > 0){
            if(missingEdgesCount == 1){
                //cand2 has to be the special neighbour.
                specialNeighbour = cand2;
            }
            else{
                //cand1 has to be the special neighbour.
                specialNeighbour = cand1;
            }
            specialNeighbourFound = true;
        }
    }

    DOUBLE_BREAK:

    if(isAlmostSimplicial){
        vertices_size_type deg_v = boost::degree(v, G);
        assert(deg_v);
        assert(low>=0);

        if(deg_v <= (vertices_size_type)low){
            vd_type vd = get_vd(G, v);

            unlink_1_neighbourhood(v, G, degs);
            degs.unlink(v);
            typename outedge_set<G_t>::type xbag;
            treedec::make_clique_and_detach(v, G, xbag);
            redegree(NULL, G, xbag, degs);
            bags.push_back(boost::make_tuple(vd, MOVE(xbag)));

            return true;
        }
        else if(vertices_size_type(low) < deg_v-1u){
            low = deg_v-1;
            return true;
        }
        else{
            return false;
        }
    }
    return false;
}

//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version stores the resulting bags in a vector and does not call further algorithms.
template <typename G_t>
void preprocessing(G_t &G, std::vector< boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename treedec_chooser<G_t>::type T_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename deg_chooser<G_t>::type degs_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(G);

    if(num_vert == 0){
        return;
    }

    degs_type degs(G);
    const degs_type& cdegs(degs);

    //Islet rule
    assert(cdegs.size());
    if(!cdegs[0].empty()){
        auto const& B=cdegs[0];
        auto I=B.begin();
        auto E=B.end();
        for(;I!=E; ++I){
            vertex_descriptor v=*I;
            auto t=boost::make_tuple(v, bag_type());
            bags.push_back(t);
        }
        low = (low > 0)? low : 0;
    }

    unsigned min_ntd = 1;
    while(boost::num_edges(G) > 0){
        if(min_ntd>1){
            --min_ntd;
        }

        vertex_descriptor v;
        boost::tie(v, min_ntd) = degs.pick_min(min_ntd, num_vert);

        if(min_ntd <= 2){
            //degree {1,2}-rules
            eliminate_vertex(v, G, bags, low, degs);
            continue;
        }else if(min_ntd==3){
            //degree 3-rules
            auto const& B=cdegs[3];
            auto it1=B.begin();
            unsigned cnt=0;
            for(; cnt<4; ++cnt){
                if(it1==B.end()){
                    break;
                }else{
                    ++it1;
                }
            }
            it1=B.begin();
            for(; it1!=B.end(); ++it1){
                //Triangle
                if(Triangle(G, *it1, bags, low, degs)){
                    goto NEXT_ITER;
                }else{
                    // graph is unchanged.
                }
                //Buddy
                auto it2=it1;
                ++it2;
                for(; it2!=B.end(); ++it2){
                    if(Buddy(G, *it1, *it2, bags, low, degs)){
                        goto NEXT_ITER;
                    }else{
                        // graph is unchanged.
                    }
                }
                if(cnt!=4){
                    // less than 4.
                    // not enough for cube rule
                }else if(Cube(G, *it1, bags, low, degs)){
                    goto NEXT_ITER;
                }else{
                    // graph is unchanged.
                }
            }
            goto ARBITRARY_DEGREE;
        }
        else{
            ARBITRARY_DEGREE:

            low = (low >= 4)? low : 4;

            for(unsigned int i = min_ntd; i < num_vert; ++i){
                auto const& B=cdegs[i];
                auto it=B.begin();
                for(; it != B.end(); ++it){
                    if(Simplicial(G, *it, bags, low, degs)){
                        goto NEXT_ITER;
                    }
                    if(AlmostSimplicial(G, *it, bags, low, degs)){
                        goto NEXT_ITER;
                    }
                }
            }
        }
        return;
NEXT_ITER:
        ;
    }
}

} //namespace impl


template <typename G_t, typename BV_t>
void preprocessing(G_t &G, BV_t &bags, int &low)
{
    impl::preprocessing(G, bags, low);
}

template <typename G_t, typename BV_t>
void preprocessing(G_t &G, BV_t &bags)
{
    int low = -1;
    preprocessing(G, bags, low);
}



} //namespace treedec

#endif //TD_PREPROCESSING

// vim:ts=8:sw=4:et
