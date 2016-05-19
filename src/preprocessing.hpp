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
// Offers functionality to preprocess a graph, such that after
// the repeated application of reduction rules, which in case the
// input graph has tree-width at most 3 allow us to determine it's tree-width exactly
// and in addition compute the corresponding tree decomposition. If the tree-width
// is larger, the reduction rules return a possibly smaller instance of the same
// tree-width as the original graph, a 'partial' tree decomposition and a lower bound
// with respect to tree-width, such that
// further algorithms can be applied to the resulting graph.
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, tree_dec_node> tree_dec_t;
//
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> graph_t;
//
//

/*
These functions are most likely to be interesting for outside use:

   -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags)
   -void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &lb)
*/

#ifndef TD_PREPROCESSING_EXP
#define TD_PREPROCESSING_EXP

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include "degree.hpp"
#include "graph.hpp"
#include "misc.hpp"

namespace treedec{

namespace impl{

// Check if there exists a degree-0-vertex.
template <typename G_t>
void Islet(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            typename noboost::treedec_traits<T_t>::vd_type vd=noboost::get_vd(G, *vIt);
            typename noboost::treedec_traits<T_t>::bag_type emptybag;

            bags.push_back(
                    boost::tuple<vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd, emptybag));

            low = (low > 0)? low : 0;
        }
    }
}

template <typename G_t>
void Islet(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags)
{
    int low = -1;
    Islet(G, bags, low);
}

/* (Islet,) Twig and Series rules. */
template <typename G_t, typename DEGS>
void eliminate_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
//    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;
    typename noboost::treedec_traits<T_t>::bag_type bag;

    assign_neighbours(bag, v, G);
    unlink_1_neighbourhood(v, G, degs);
    degs.unlink(v);

    unsigned deg = boost::degree(v, G);
    typename outedge_set<G_t>::type xbag;
    treedec::make_clique_and_detach(v, G, xbag);
    xbag.clear(); // provide interface with clear included? (not urgent)
    assert(!boost::degree(v, G));

    redegree(NULL, G, bag, degs);

    bags.push_back(boost::make_tuple(v, bag));

    low = (low > (int)deg)? low : deg;
}

//Applies the Triangle rule if applicable (checks if there exists a degree-3-vertex,
//such that at least one edge exists in its neighbourhood).
template <typename G_t, typename DEGS>
bool Triangle(G_t &G,
              typename boost::graph_traits<G_t>::vertex_descriptor v,
              std::vector<boost::tuple<
                typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
                typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
              > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typename noboost::treedec_traits<T_t>::bag_type bag;

    assign_neighbours(bag, v, G);

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(3);
    BOOST_AUTO(f, boost::adjacent_vertices(v, G).first);
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
        xbag.clear(); // provide interface with clear included? (not urgent)
        redegree(NULL, G, bag, degs);

        bags.push_back( boost::make_tuple(v, bag));

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
             typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
             typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
           > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

    typename noboost::treedec_traits<T_t>::bag_type N1, N2;
    assign_neighbours(N1, v, G);
    assign_neighbours(N2, w, G);

    if(N1 == N2){
        unlink_1_neighbourhood(v, G, degs);
        degs.unlink(v, 3);
        degs.unlink(w, 3);
        noboost::make_clique(boost::adjacent_vertices(v, G), G);
        boost::clear_vertex(v, G);
        boost::clear_vertex(w, G);
        redegree(NULL, G, N1, degs);

        vd_type vd1 = noboost::get_vd(G, v);
        vd_type vd2 = noboost::get_vd(G, w);

        bags.push_back(boost::tuple<vd_type, bag_type>(vd1, N1));
        bags.push_back(boost::tuple<vd_type, bag_type>(vd2, N2));

        low = (low > 3)? low : 3;
        return true;
    }
    return false;
}

//Applies the Cube rule if applicable.
template <typename G_t, typename DEGS>
bool Cube(G_t &G,
          typename boost::graph_traits<G_t>::vertex_descriptor vertex,
          std::vector<boost::tuple<
            typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
            typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
           > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

    typename boost::graph_traits<G_t>::vertex_descriptor x, a, b, c;
    x = vertex;

    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Na(3), Nb(3), Nc(3);

    a = *(boost::adjacent_vertices(x, G).first);
    b = *(++boost::adjacent_vertices(x, G).first);
    c = *(++(++boost::adjacent_vertices(x, G).first));

    if(boost::degree(a, G) != 3 || boost::degree(b, G) != 3 || boost::degree(c, G) != 3){
        return false;
    }

    Na[0] = *(boost::adjacent_vertices(a, G).first);
    Na[1] = *(++boost::adjacent_vertices(a, G).first);
    Na[2] = *(++(++boost::adjacent_vertices(a, G).first));
    Nb[0] = *(boost::adjacent_vertices(b, G).first);
    Nb[1] = *(++boost::adjacent_vertices(b, G).first);
    Nb[2] = *(++(++boost::adjacent_vertices(b, G).first));
    Nc[0] = *(boost::adjacent_vertices(c, G).first);
    Nc[1] = *(++boost::adjacent_vertices(c, G).first);
    Nc[2] = *(++(++boost::adjacent_vertices(c, G).first));

    if(Na[0] == x){
        Na[0] = Na[2];
    }
    else if(Na[1] == x){
        Na[1] = Na[2];
    }
    else{}
    if(Nb[0] == x){
        Nb[0] = Nb[2];
    }
    else if(Na[1] == x){
        Nb[1] = Nb[2];
    }
    else{}
    if(Nc[0] == x){
        Nc[0] = Nc[2];
    }
    else if(Nc[1] == x){
        Nc[1] = Nc[2];
    }
    else{}

    typename noboost::treedec_traits<T_t>::bag_type bag;
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
        vd_type vdx = noboost::get_vd(G, x);
        vd_type vda = noboost::get_vd(G, a);
        vd_type vdb = noboost::get_vd(G, b);
        vd_type vdc = noboost::get_vd(G, c);
        vd_type vdu = noboost::get_vd(G, u);
        vd_type vdv = noboost::get_vd(G, v);
        vd_type vdw = noboost::get_vd(G, w);

        bag.insert(vdu); bag.insert(vdv); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vda, bag));
        bag.clear();

        bag.insert(vdw); bag.insert(vdv); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vdc, bag));
        bag.clear();

        bag.insert(vdw); bag.insert(vdu); bag.insert(vdx);
        bags.push_back(boost::make_tuple(vdb, bag));

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
        boost::add_edge(v, w, G);
        boost::add_edge(u, x, G);
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
                  typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
                  typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
                > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

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
        bag_type bag;
        assign_neighbours(bag, v, G);

        vd_type vd = noboost::get_vd(G, v);
        bags.push_back(boost::make_tuple(vd, bag));

        unlink_1_neighbourhood(v, G, degs);
        degs.unlink(v);
        typename outedge_set<G_t>::type xbag;
        treedec::make_clique_and_detach(v, G, xbag);
        xbag.clear(); // provide interface with clear included? (not urgent)
        redegree(NULL, G, bag, degs);

        low = (low > (int)bag.size())? low : (int)bag.size();
        return true;
    }
    return false;
}


//Applies the Almost Simplicial rule if possible.
template <typename G_t, typename DEGS>
bool AlmostSimplicial(G_t &G,
                      typename boost::graph_traits<G_t>::vertex_descriptor v,
                      std::vector<boost::tuple<
                        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
                        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
                      > > &bags, int &low, DEGS &degs)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

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
        int deg_v = (int)boost::degree(v, G);

        if(deg_v <= low){
            bag_type bag;
            assign_neighbours(bag, v, G);
            vd_type vd = noboost::get_vd(G, v);

            bags.push_back(boost::tuple<vd_type, bag_type>(vd, bag));

            unlink_1_neighbourhood(v, G, degs);
            degs.unlink(v);
            typename outedge_set<G_t>::type xbag;
            treedec::make_clique_and_detach(v, G, xbag);
            xbag.clear(); // provide interface with clear included? (not urgent)
            redegree(NULL, G, bag, degs);

            return true;
        }
        else if(low < deg_v-1){
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
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename noboost::deg_chooser<G_t>::type degs_type;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

    // incomplete?
    //bag_type bag_i;
    //bag_type* bags_i=&bag_i;

    degs_type degs(G);
    const degs_type& cdegs(degs);
    typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(G);

    //Islet rule
    if(!cdegs[0].empty()){
        bag_type emptybag;
        BOOST_AUTO(I, cdegs[0].begin());
        BOOST_AUTO(E, cdegs[0].end());
        for(; I != E; I++){
            bags.push_back(boost::tuple<vd_type, bag_type>(*I, emptybag));
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

        bool reduction_complete = true;

        //degree {1,2}-rules
        if(min_ntd <= 2){
            eliminate_vertex(v, G, bags, low, degs);
            reduction_complete = false;
        }
        //degree 3-rules
        else if(min_ntd == 3){
            BOOST_AUTO(it1, cdegs[3].begin());
            for(; it1!=cdegs[3].end(); ++it1){
                //Triangle
                if(Triangle(G, *it1, bags, low, degs)){
                    reduction_complete = false;
                    goto NEXT_ITER;
                }
                //Buddy
                BOOST_AUTO(it2, it1);
                it2++;
                for(; it2 != cdegs[3].end(); ++it2){
                    if(Buddy(G, *it1, *it2, bags, low, degs)){
                        reduction_complete = false;
                        goto NEXT_ITER;
                    }
                }
                //Cube
                if(cdegs[3].size() >= 4 && Cube(G, *it1, bags, low, degs)){
                    reduction_complete = false;
                    goto NEXT_ITER;
                }
            }
            goto ARBITRARY_DEGREE;
        }
        else{
            ARBITRARY_DEGREE:

            low = (low >= 4)? low : 4;

            for(unsigned int i = min_ntd; i < num_vert; ++i){
                BOOST_AUTO(it, cdegs[i].begin());
                for(; it != cdegs[i].end(); ++it){
                    if(Simplicial(G, *it, bags, low, degs)){
                        reduction_complete = false;
                        goto NEXT_ITER;
                    }
                    if(AlmostSimplicial(G, *it, bags, low, degs)){
                        reduction_complete = false;
                        goto NEXT_ITER;
                    }
                }
            }
        }
        NEXT_ITER:
        if(reduction_complete){
            return;
        }
    }
}

} //namespace impl


template <typename G_t>
void preprocessing(G_t &G, std::vector< boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
             > > &bags, int &low)
{
    impl::preprocessing(G, bags, low);
}

template <typename G_t>
void preprocessing(G_t &G, std::vector< boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
             > > &bags)
{
    int low = -1;
    preprocessing(G, bags, low);
}



} //namespace treedec

#endif //TD_PREPROCESSING_EXP

// vim:ts=8:sw=4:et
