// Lukas Larisch, 2014 - 2016
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
   -void preprocessing_glue_bags(std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, T_t &T)
*/

#ifndef TD_PREPROCESSING
#define TD_PREPROCESSING

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include "TD_noboost.hpp"
#include "TD_misc.hpp"

namespace treedec{

//Checks if there exists a degree-0-vertex.
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

//Checks if there exists a degree-0-vertex.
template <typename G_t>
void Islet(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            typename noboost::treedec_traits<T_t>::vd_type vd=noboost::get_vd(G, *vIt);
            typename noboost::treedec_traits<T_t>::bag_type emptybag;

            bags.push_back(
                    boost::tuple<
                    typename noboost::treedec_traits<T_t>::vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd, emptybag));
        }
    }
}

//Checks if there exists a degree-1-vertex.
template <typename G_t>
bool Twig(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 1){
            typename noboost::treedec_traits<T_t>::vd_type vd=noboost::get_vd(G, *vIt);
            typename noboost::treedec_traits<T_t>::bag_type bag;

            noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(*vIt, G), G);

            bags.push_back(
                    boost::tuple<
                    typename noboost::treedec_traits<T_t>::vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > 1)? low : 1;
            return true;
        }
    }
    return false;
}

//Checks if there exists a degree-2-vertex.
template <typename G_t>
bool Series(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 2){
            typename noboost::treedec_traits<T_t>::vd_type vd=noboost::get_vd(G, *vIt);
            typename noboost::treedec_traits<T_t>::bag_type bag;

            noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(*vIt, G), G);
            noboost::make_clique(boost::adjacent_vertices(*vIt, G), G);

            bags.push_back(
                    boost::tuple<
                    typename noboost::treedec_traits<T_t>::vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > 2)? low : 2;
            return true;
        }
    }
    return false;
}

//Checks if there exists a degree-3-vertex, such that at least one edge exists in its neighbourhood (Triangle).
template <typename G_t>
bool Triangle(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 3){
            typename noboost::treedec_traits<T_t>::bag_type bag;

            noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(*vIt, G), G);

            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(3);
            N[0] = *(boost::adjacent_vertices(*vIt, G).first);
            N[1] = *(++boost::adjacent_vertices(*vIt, G).first);
            N[2] = *(++(++boost::adjacent_vertices(*vIt, G).first));

            if(boost::edge(N[0], N[1], G).second
            || boost::edge(N[0], N[2], G).second
            || boost::edge(N[1], N[2], G).second)
            {

                noboost::make_clique(boost::adjacent_vertices(*vIt, G), G);

                typename noboost::treedec_traits<T_t>::vd_type vd=noboost::get_vd(G, *vIt);

                bags.push_back(
                    boost::tuple<
                    typename noboost::treedec_traits<T_t>::vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd, bag));

                boost::clear_vertex(*vIt, G);

                low = (low > 3)? low : 3;
                return true;
            }
        }
    }
    return false;
}

//Checks if there exists two degree-3-vertices, such that they share their neighbours (Buddies).
template <typename G_t>
bool Buddy(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type vd_type;
    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vEnd;

    for(boost::tie(vIt1, vEnd) = boost::vertices(G); vIt1 != vEnd; vIt1++){
        if(boost::degree(*vIt1, G) != 3){
            continue;
        }

        typename noboost::treedec_traits<T_t>::bag_type N1;

        noboost::fetch_neighbourhood(N1, boost::adjacent_vertices(*vIt1, G), G);

        typename boost::graph_traits<G_t>::vertex_iterator vIt2;
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(boost::degree(*vIt2, G) != 3){
                continue;
            }

            typename noboost::treedec_traits<T_t>::bag_type N2;
            noboost::fetch_neighbourhood(N2, boost::adjacent_vertices(*vIt2, G), G);

            if(N1 == N2){
                std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(3);
                N[0] = *(boost::adjacent_vertices(*vIt1, G).first);
                N[1] = *(++boost::adjacent_vertices(*vIt1, G).first);
                N[2] = *(++(++boost::adjacent_vertices(*vIt1, G).first));

                noboost::make_clique(boost::adjacent_vertices(*vIt1, G), G);

                vd_type vd1=noboost::get_vd(G, *vIt1);
                vd_type vd2=noboost::get_vd(G, *vIt2);

                bags.push_back(
                    boost::tuple<vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd1, N1));

                bags.push_back(
                    boost::tuple<vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vd2, N2));

                boost::clear_vertex(*vIt1, G);
                boost::clear_vertex(*vIt2, G);

                low = (low > 3)? low : 3;
                return true;
            }
        }
    }
    return false;
}

//Checks if the Cube rule is applicable.
template <typename G_t>
bool Cube(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) != 3){
            continue;
        }

        typename boost::graph_traits<G_t>::vertex_descriptor x,a,b,c;
        x = *vIt;

        typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx(3);
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Na(3), Nb(3), Nc(3);

        Nx[0] = *(boost::adjacent_vertices(x, G).first);
        Nx[1] = *(++boost::adjacent_vertices(x, G).first);
        Nx[2] = *(++(++boost::adjacent_vertices(x, G).first));

        a = Nx[0];
        b = Nx[1];
        c = Nx[2];

        if(boost::degree(a, G) != 3 || boost::degree(b, G) != 3 || boost::degree(c, G) != 3){
            continue;
        }

        typename noboost::treedec_traits<T_t>::bag_type bag;

        Na[0] = *(boost::adjacent_vertices(a, G).first);
        Na[1] = *(++boost::adjacent_vertices(a, G).first);
        Na[2] = *(++(++boost::adjacent_vertices(a, G).first));
        Nb[0] = *(boost::adjacent_vertices(b, G).first);
        Nb[1] = *(++boost::adjacent_vertices(b, G).first);
        Nb[2] = *(++(++boost::adjacent_vertices(b, G).first));
        Nc[0] = *(boost::adjacent_vertices(c, G).first);
        Nc[1] = *(++boost::adjacent_vertices(c, G).first);
        Nc[2] = *(++(++boost::adjacent_vertices(c, G).first));

        noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(a, G), G);
        noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(b, G), G);
        noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(c, G), G);

        if(bag.size() != 4){
            bag.clear();
            continue;
        }

        typename boost::graph_traits<G_t>::vertex_descriptor u,v,w;

        if(Na[0] == Nb[0]){      u = Na[0]; v = Na[1]; w = Nb[1]; }
        else if(Na[0] == Nb[1]){ u = Na[0]; v = Na[1]; w = Nb[0]; }
        else if(Na[1] == Nb[0]){ u = Na[1]; v = Na[0]; w = Nb[1]; }
        else if(Na[1] == Nb[1]){ u = Na[1]; v = Na[0]; w = Nb[0]; }
        else{ continue; }

        if((Nc[0] == v && Nc[1] == w) || (Nc[1] == v && Nc[0] == w)){
            bag.clear();
            vd_type vdu=noboost::get_vd(G, u);
            vd_type vdv=noboost::get_vd(G, v);
            vd_type vdw=noboost::get_vd(G, w);
            vd_type vdx=noboost::get_vd(G, x);
            vd_type vda=noboost::get_vd(G, a);
            vd_type vdb=noboost::get_vd(G, b);
            vd_type vdc=noboost::get_vd(G, c);

            bag.insert(vdu); bag.insert(vdv); bag.insert(vdx);

            bags.push_back(
                boost::tuple<vd_type,
                typename noboost::treedec_traits<T_t>::bag_type
                >(vda, bag));

            bag.clear();

            bag.insert(vdw); bag.insert(vdv); bag.insert(vdx);

                bags.push_back(
                    boost::tuple<
                    typename noboost::treedec_traits<T_t>::vd_type,
                    typename noboost::treedec_traits<T_t>::bag_type
                    >(vdc, bag));

            bag.clear();

            bag.insert(vdw); bag.insert(vdu); bag.insert(vdx);

            bags.push_back(
                boost::tuple<
                typename noboost::treedec_traits<T_t>::vd_type,
                typename noboost::treedec_traits<T_t>::bag_type
                >(vdb, bag));

            boost::clear_vertex(a, G);
            boost::clear_vertex(b, G);
            boost::clear_vertex(c, G);

            boost::add_edge(u, v, G);
            boost::add_edge(u, w, G);
            boost::add_edge(v, w, G);
            boost::add_edge(u, x, G);
            boost::add_edge(v, x, G);
            boost::add_edge(w, x, G);

            low = (low > 3)? low : 3;
            return true;
        }
    }
    return false;
}

//Checks if there exists a vertex, such that its neighbours induce a clique (Simplicial).
template <typename G_t>
bool Simplicial(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            continue;
        }

        //The neighbourhood of *vIt is a clique, if no "edge miss" occures.
        bool isClique = true;

        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
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
            typename noboost::treedec_traits<T_t>::bag_type bag;
            noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(*vIt, G), G);

            vd_type vd=noboost::get_vd(G, *vIt);
            bags.push_back(
                boost::tuple<vd_type,
                typename noboost::treedec_traits<T_t>::bag_type
                >(vd, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > (int)bag.size())? low : (int)bag.size();
            return true;
        }
    }
    return false;
}


//Checks if there exists an almost simplicial vertex in G.
template <typename G_t>
bool AlmostSimplicial(G_t &G, std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::vd_type vd_type;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            continue;
        }

        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(boost::degree(*vIt, G)+1);
        N[0] = *vIt;

        unsigned int c = 0;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            N[++c] = *nIt;
        }

        //N except one vertex now potentially is a clique.
        bool isAlmostSimplicial = true;
        bool specialNeighbourFound = false;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator nIt1, nIt2;
        typename boost::graph_traits<G_t>::vertex_descriptor cand1, cand2, specialNeighbour;
        unsigned int missingEdgesCount;

        for(nIt1 = N.begin(); nIt1 != N.end(); nIt1++){
            nIt2 = nIt1;
            nIt2++;
            missingEdgesCount = 0;
            for(; nIt2 != N.end(); nIt2++){
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
            //Adding the edges, if specialNeighbourFound is true, N is a clique and *vIt is a simplicial vertex.
            if(specialNeighbourFound){
                for(unsigned int i = 0; i < N.size(); i++){
                    if(N[i] != specialNeighbour){
                        boost::add_edge(specialNeighbour, N[i], G);
                    }
                }
            }

            typename noboost::treedec_traits<T_t>::bag_type bag;
            noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(*vIt, G), G);

            vd_type vd=noboost::get_vd(G, *vIt);
            bags.push_back(
                boost::tuple<vd_type,
                typename noboost::treedec_traits<T_t>::bag_type
                >(vd, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > (int)bag.size())? low : (int)bag.size();
            return true;
        }
    }
    return false;
}

//Glues a single bag with the current tree decomposition T according to subset relation.
template<typename T_t>
void glue_bag_preprocessing(
        typename noboost::treedec_traits<T_t>::bag_type &bag,
        typename noboost::treedec_traits<T_t>::vd_type preprocessed_node,
        T_t &T)
{
    if(boost::num_vertices(T) == 0){
        bag.insert(preprocessed_node);
        typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
        noboost::bag(t_dec_node, T) = MOVE(bag);

        return;
    }

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        if(std::includes(noboost::bag(*vIt, T).begin(),
                         noboost::bag(*vIt, T).end(),
                         bag.begin(), bag.end()))
        {
            bag.insert(preprocessed_node);
            typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
            noboost::bag(t_dec_node, T) = MOVE(bag);

            boost::add_edge(*vIt, t_dec_node, T);
            return;
        }
    }

    //Case for a disconnected graph.
    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
    bag.insert(preprocessed_node);
    noboost::bag(t_dec_node, T) = MOVE(bag);
    boost::tie(vIt, vEnd) = boost::vertices(T);
    boost::add_edge(*vIt, t_dec_node, T);

}

//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version stores the resulting bags in a vector and does not call further algorithms.
template <typename G_t>
void _preprocessing(G_t &G, std::vector< boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    if(Twig(G, bags, low) || Series(G, bags, low) || Triangle(G, bags, low)
    || Buddy(G, bags, low) || Cube(G, bags, low) || Simplicial(G, bags, low)
    || AlmostSimplicial(G, bags, low)){

        _preprocessing(G, bags, low);
    }
    else if(boost::num_edges(G) != 0){
        low = (low > 4)? low : 4;
    }
}

template <typename G_t>
void preprocessing(G_t &G, std::vector< boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
              > > &bags, int &low)
{
    Islet(G, bags);
    _preprocessing(G, bags, low);
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

//Glues bags with the current tree decomposition.
template<typename T_t>
void preprocessing_glue_bags(std::vector< boost::tuple<
        typename noboost::treedec_traits<T_t>::vd_type,
        typename noboost::treedec_traits<T_t>::bag_type
             > > &bags, T_t &T)
{
    for(unsigned int i = bags.size(); i > 0; i--){
        typename noboost::treedec_traits<T_t>::vd_type first = boost::get<0>(bags[i-1]);
        typename noboost::treedec_traits<T_t>::bag_type& second = boost::get<1>(bags[i-1]);
        glue_bag_preprocessing(second, first, T);
    }
}

} //namespace treedec

#endif //TD_PREPROCESSING

// vim:ts=8:sw=4:et
