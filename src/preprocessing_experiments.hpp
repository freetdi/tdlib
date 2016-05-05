// Lukas Larisch, 2014 - 2015
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
// Offers functionality to preprocess a graph, such that after
// the repeated application of reduction rules, which in case the
// input graph has tree-width at most 3 allow us to determine it's tree-width exactly
// and in addition compute the corresponding tree decomposition. If the tree-width
// is larger, the reduction rules return a possibly smaller instance of the same
// tree-width as the original graph, a partial tree decomposition and a lower bound
// with respect to tree-width, such that
// further algorithms can be applied to the resulting graph.
// Currently, the minDegree-heuristic will be applied to the resulting graph,
// if the input graph can't be fully preprocessed.
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, tree_dec_node> tree_dec_t;
//
// These functions are most likely to be interesting for outside use:
//
// void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags)
// void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &lb)
// void preprocessing_glue_bags(std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, T_t &T)
//

#ifndef TD_PREPROCESSING_EXP
#define TD_PREPROCESSING_EXP

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include "degree.hpp"
#include "graph.hpp"
#include "misc.hpp"

namespace treedec{

namespace exp{

/* Islet, Twig and Series rules. */
template <typename G_t>
void eliminate_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

    typename noboost::treedec_traits<T_t>::bag_type bag;

    noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(v, G), G);
    unsigned int deg = noboost::eliminate_vertex(v, G);

    bags.push_back(
             boost::tuple<
              typename noboost::treedec_traits<T_t>::vd_type,
              typename noboost::treedec_traits<T_t>::bag_type
             >(v, bag));

    low = (low > deg)? low : deg;
}

//Checks if there exists a degree-3-vertex, such that at least one edge exists in its neighbourhood (Triangle).
template <typename G_t>
bool Triangle(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

    typename noboost::treedec_traits<T_t>::bag_type bag;

    noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(v, G), G);

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(3);
    N[0] = *(boost::adjacent_vertices(v, G).first);
    N[1] = *(++boost::adjacent_vertices(v, G).first);
    N[2] = *(++(++boost::adjacent_vertices(v, G).first));

    if(boost::edge(N[0], N[1], G).second
    || boost::edge(N[0], N[2], G).second
    || boost::edge(N[1], N[2], G).second)
    {
        noboost::make_clique(boost::adjacent_vertices(v, G), G);

        bags.push_back(
           boost::tuple<
        typename noboost::treedec_traits<T_t>::vd_type,
        typename noboost::treedec_traits<T_t>::bag_type
         >(v, bag));

        boost::clear_vertex(v, G);

        low = (low > 3)? low : 3;
        return true;
    }

    return false;
}

//Checks if there exists two degree-3-vertices, such that they share their neighbours (Buddies).
template <typename G_t>
bool Buddy(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v,
                   typename boost::graph_traits<G_t>::vertex_descriptor w,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

    typename noboost::treedec_traits<T_t>::bag_type N1, N2;
    noboost::fetch_neighbourhood(N1, boost::adjacent_vertices(v, G), G);
    noboost::fetch_neighbourhood(N2, boost::adjacent_vertices(w, G), G);

    if(N1 == N2){
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(3);
        N[0] = *(boost::adjacent_vertices(v, G).first);
        N[1] = *(++boost::adjacent_vertices(v, G).first);
        N[2] = *(++(++boost::adjacent_vertices(v, G).first));

        noboost::make_clique(boost::adjacent_vertices(v, G), G);

        bags.push_back(
           boost::tuple<
        typename noboost::treedec_traits<T_t>::vd_type,
        typename noboost::treedec_traits<T_t>::bag_type
         >(v, N1));

        bags.push_back(
            boost::tuple<
        typename noboost::treedec_traits<T_t>::vd_type,
        typename noboost::treedec_traits<T_t>::bag_type
         >(w, N2));

        boost::clear_vertex(v, G);
        boost::clear_vertex(w, G);

        low = (low > 3)? low : 3;
        return true;
    }
    return false;
}

//Checks if the Cube rule is applicable.
template <typename G_t>
bool Cube(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor vertex,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

    typename boost::graph_traits<G_t>::vertex_descriptor x, a, b, c;
    x = vertex;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx(3);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Na(3), Nb(3), Nc(3);

    Nx[0] = *(boost::adjacent_vertices(x, G).first);
    Nx[1] = *(++boost::adjacent_vertices(x, G).first);
    Nx[2] = *(++(++boost::adjacent_vertices(x, G).first));

    a = Nx[0];
    b = Nx[1];
    c = Nx[2];

    if(boost::degree(a, G) != 3 || boost::degree(b, G) != 3 || boost::degree(c, G) != 3){
        return false;
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
        return false;
    }

    typename boost::graph_traits<G_t>::vertex_descriptor u, v, w;

    if(Na[0] == Nb[0]){      u = Na[0]; v = Na[1]; w = Nb[1]; }
    else if(Na[0] == Nb[1]){ u = Na[0]; v = Na[1]; w = Nb[0]; }
    else if(Na[1] == Nb[0]){ u = Na[1]; v = Na[0]; w = Nb[1]; }
    else if(Na[1] == Nb[1]){ u = Na[1]; v = Na[0]; w = Nb[0]; }
    else{
        return false;
    }

    if((Nc[0] == v && Nc[1] == w) || (Nc[1] == v && Nc[0] == w)){
        bag.insert(u); bag.insert(v); bag.insert(x);

        bags.push_back(
            boost::tuple<
            typename noboost::treedec_traits<T_t>::vd_type,
            typename noboost::treedec_traits<T_t>::bag_type
            >(a, bag));

        bag.clear();

        bag.insert(w); bag.insert(v); bag.insert(x);

        bags.push_back(
          boost::tuple<
          typename noboost::treedec_traits<T_t>::vd_type,
          typename noboost::treedec_traits<T_t>::bag_type
         >(c, bag));

        bag.clear();

        bag.insert(w); bag.insert(u); bag.insert(x);

        bags.push_back(
          boost::tuple<
          typename noboost::treedec_traits<T_t>::vd_type,
          typename noboost::treedec_traits<T_t>::bag_type
         >(b, bag));

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
    return false;
}

//Checks if there exists a vertex, such that its neighbours induce a clique (Simplicial).
template <typename G_t>
bool Simplicial(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

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
        typename noboost::treedec_traits<T_t>::bag_type bag;
        noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(v, G), G);

        bags.push_back(
            boost::tuple<
            typename noboost::treedec_traits<T_t>::vd_type,
            typename noboost::treedec_traits<T_t>::bag_type
           >(v, bag));

        boost::clear_vertex(v, G);

        low = (low > (int)bag.size())? low : (int)bag.size();
        return true;
    }
    return false;
}


//Checks if there exists an almost simplicial vertex in G.
template <typename G_t>
bool AlmostSimplicial(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v,
          std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N(boost::degree(v, G)+1);
    N[0] = v;

    unsigned int c = 0;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
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
        noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(v, G), G);

        bags.push_back(
            boost::tuple<
            typename noboost::treedec_traits<T_t>::vd_type,
            typename noboost::treedec_traits<T_t>::bag_type
           >(v, bag));

        boost::clear_vertex(v, G);

        low = (low > (int)bag.size())? low : (int)bag.size();
        return true;
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

#ifndef REDEGREE
#define REDEGREE

//register a 1-neigborhood to DEGS
template<class U, class G_t, class B, class D>
void redegree(U, G_t &G, B& neighborhood, D& degree)
{
    BOOST_AUTO(I, neighborhood.begin());
    BOOST_AUTO(E, neighborhood.end());

    for(; I != E ; ++I){
        size_t deg = boost::degree(*I, G);
        degree.reg(*I, deg);
    }
}

#endif

//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version stores the resulting bags in a vector and does not call further algorithms.
template <typename G_t>
void _preprocessing(G_t &G, std::vector< boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > &bags, int &low)
{

    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename noboost::deg_chooser<G_t>::type degs_type;
    typedef typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type bag_type;
    bag_type bag_i;
    bag_type* bags_i=&bag_i;

    degs_type degs(G);
    typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(G);

    unsigned min_ntd = 1;
    while(boost::num_edges(G) > 0){
        if(min_ntd>1){
            --min_ntd;
        }

        vertex_descriptor v;
        
        boost::tie(v, min_ntd) = degs.pick_min(min_ntd, num_vert);
        noboost::make_clique_and_hijack(v, G, (void*)NULL, *bags_i);

        if(min_ntd <= 2){
            typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; ++nIt){
                degs.unlink(*nIt);
            }
            eliminate_vertex(v, G, bags, low);
            degs.unlink(v);
            redegree(NULL, G, *bags_i, degs);
        }
        else if(min_ntd == 3){
            for(typename std::unordered_set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator it
                    = degs[3].begin(); it != degs[3].end(); ++it)
            {
                if(Triangle(G, *it, bags, low)){
                    goto END;
                }
            }
            for(typename std::unordered_set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator it
                    = degs[3].begin(); it != degs[3].end(); ++it)
            {
                typename std::unordered_set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator it2 = it;
                ++it2;
                for(; it2 != degs[3].end(); ++it){
                    if(Buddy(G, *it, *it2, bags, low)){
                        goto END;
                    }
                }
            }
            for(typename std::unordered_set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator it
               = degs[3].begin(); it != degs[3].end(); ++it)
            {
                if(Cube(G, *it, bags, low)){
                    goto END;
                }
            }
            goto ARBITRARY_DEGREE;
        }
        else{
            ARBITRARY_DEGREE:

            size_t c = min_ntd;
            while(c < num_vert){
                for(typename std::unordered_set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator it
                    = degs[c].begin(); it != degs[c].end(); ++it)
                {
                    if(Simplicial(G, *it, bags, low) || AlmostSimplicial(G, *it, bags, low)){
                        goto END;
                    }
                }
            }
            return;
        }
        END:
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
        glue_bag_preprocessing(bags[i-1].get<1>(), bags[i-1].get<0>(), T);
    }
}

}

} //namespace treedec

#endif //TD_PREPROCESSING

// vim:ts=8:sw=4:et
