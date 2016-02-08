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
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, tree_dec_node> tree_dec_t;
//
// Vertices of the input graph have to provide the attribute 'id', e.g.:
//
// struct Vertex
// {
//  unsigned int id;
// };
// typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, Vertex> TD_graph_t;
//
//
//
// These functions are most likely to be interesting for outside use:
//
// void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags)
// void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &lb)
// void preprocessing_glue_bags(std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, T_t &T)
//

#ifndef TD_PREPROCESSING
#define TD_PREPROCESSING

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include "TD_misc.hpp"

namespace treedec{

//checks if there exists a degree-0-vertex
template <typename G_t>
void Islet(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(G); vertexIt != vertexEnd; vertexIt++){
        if(boost::out_degree(*vertexIt, G) == 0){
            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt].id, std::set<unsigned int>()));

            low = (low > 0)? low : 0;
        }
    }
}

//checks if there exists a degree-0-vertex
template <typename G_t>
void Islet(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags){
    typename boost::graph_traits<G_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(G); vertexIt != vertexEnd; vertexIt++){
        if(boost::out_degree(*vertexIt, G) == 0)
            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt].id, std::set<unsigned int>()));
    }
}

//checks if there exists a degree-1-vertex
template <typename G_t>
bool Twig(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(G); vertexIt != vertexEnd; vertexIt++){
        if(boost::out_degree(*vertexIt, G) == 1){
            std::set<unsigned int> bag;

            typename boost::graph_traits<G_t>::adjacency_iterator neighbourIt, neighbourEnd;
            for(boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(*vertexIt, G); neighbourIt != neighbourEnd; neighbourIt++)
                bag.insert(G[*neighbourIt].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt].id, bag));

            boost::clear_vertex(*vertexIt, G);

            low = (low > 1)? low : 1;
            return true;
        }
    }
    return false;
}

//checks if there exists a degree-2-vertex
template <typename G_t>
bool Series(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(G); vertexIt != vertexEnd; vertexIt++){
        if(boost::out_degree(*vertexIt, G) == 2){
            std::set<unsigned int> bag;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor > N;

            typename boost::graph_traits<G_t>::adjacency_iterator  neighbourIt, neighbourEnd;
            for(boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(*vertexIt, G); neighbourIt != neighbourEnd; neighbourIt++){
                bag.insert(G[*neighbourIt].id);
                N.push_back(*neighbourIt);
            }

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt].id, bag));

            boost::clear_vertex(*vertexIt, G);
            boost::add_edge(N.at(0), N.at(1), G);

            low = (low > 2)? low : 2;
            return true;
        }
    }
    return false;
}

//checks if there exists a degree-3-vertex, such that at least one edge exists in its neighbourhood (Triangle)
template <typename G_t>
bool Triangle(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(G); vertexIt != vertexEnd; vertexIt++){
        if(boost::out_degree(*vertexIt, G) == 3){
            std::set<unsigned int> bag;
            std::vector<typename boost::graph_traits<G_t>::adjacency_iterator> N;

            typename boost::graph_traits<G_t>::adjacency_iterator neighbourIt, neighbourEnd;
            for(boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(*vertexIt, G); neighbourIt != neighbourEnd; neighbourIt++){
                bag.insert(G[*neighbourIt].id);
                N.push_back(neighbourIt);
            }

            if(boost::edge(*N.at(0), *N.at(1), G).second || boost::edge(*N.at(0), *N.at(2), G).second || boost::edge(*N.at(1), *N.at(2), G).second){
                boost::add_edge(*N.at(0), *N.at(1), G);
                boost::add_edge(*N.at(0), *N.at(2), G);
                boost::add_edge(*N.at(1), *N.at(2), G);

                bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt].id, bag));

                boost::clear_vertex(*vertexIt, G);

                low = (low > 3)? low : 3;
                return true;
            }
        }
    }
    return false;
}

//checks if there exists two degree-3-vertices, such that they share their neighbours (Buddies)
template <typename G_t>
bool Buddy(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator  vertexIt1, vertexEnd;

    for(boost::tie(vertexIt1, vertexEnd) = boost::vertices(G); vertexIt1 != vertexEnd; vertexIt1++){
        if(boost::out_degree(*vertexIt1, G) != 3)
            continue;

        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> dN;
        std::set<unsigned int> N1;

        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd1;
        for(boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(*vertexIt1, G); nIt1 != nEnd1; nIt1++){
            dN.push_back(*nIt1);
            N1.insert(G[*nIt1].id);
        }

        typename boost::graph_traits<G_t>::vertex_iterator vertexIt2;
        vertexIt2 = vertexIt1;
        vertexIt2++;
        for(; vertexIt2 != vertexEnd; vertexIt2++){
            if(boost::out_degree(*vertexIt2, G) != 3)
                continue;

            std::set<unsigned int> N2;
            typename boost::graph_traits<G_t>::adjacency_iterator nIt2, nEnd2;
            for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(*vertexIt2, G); nIt2 != nEnd2; nIt2++){
                N2.insert(G[*nIt2].id);
            }
            if(N1 == N2){
                boost::add_edge(dN.at(0), dN.at(1), G);
                boost::add_edge(dN.at(0), dN.at(2), G);
                boost::add_edge(dN.at(1), dN.at(2), G);

                bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt1].id, N1));
                bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vertexIt2].id, N2));

                boost::clear_vertex(*vertexIt1, G);
                boost::clear_vertex(*vertexIt2, G);

                low = (low > 3)? low : 3;
                return true;
            }
        }
    }
    return false;
}

//checks if the Cube rule is applicable
template <typename G_t>
bool Cube(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) != 3)
            continue;

        typename boost::graph_traits<G_t>::vertex_descriptor x,a,b,c;
        x = *vIt;

        typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx, Na, Nb, Nc;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
            Nx.push_back(*nIt);

        a = Nx[0];
        b = Nx[1];
        c = Nx[2];

        if(boost::out_degree(a, G) != 3 || boost::out_degree(b, G) != 3 || boost::out_degree(c, G) != 3)
            continue;

        std::set<unsigned int> bag;

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(a, G); nIt != nEnd; nIt++){
            if(*nIt != x){
                Na.push_back(*nIt);
                bag.insert(G[*nIt].id);
            }
        }

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(b, G); nIt != nEnd; nIt++){
            if(*nIt != x){
                Nb.push_back(*nIt);
                bag.insert(G[*nIt].id);
            }
        }

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(c, G); nIt != nEnd; nIt++){
            if(*nIt != x){
                Nc.push_back(*nIt);
                bag.insert(G[*nIt].id);
            }
        }

        if(bag.size() != 4){
            bag.clear();
            continue;
        }

        typename boost::graph_traits<G_t>::vertex_descriptor u,v,w;

        if(Na[0] == Nb[0]){
            u = Na[0];
            v = Na[1];
            w = Nb[1];
        }
        else if(Na[0] == Nb[1]){
            u = Na[0];
            v = Na[1];
            w = Nb[0];
        }
        else if(Na[1] == Nb[0]){
            u = Na[1];
            v = Na[0];
            w = Nb[1];
        }
        else if(Na[1] == Nb[1]){
            u = Na[1];
            v = Na[0];
            w = Nb[0];
        }
        else
            continue;

        if((Nc[0] == v && Nc[1] == w) || (Nc[1] == v && Nc[0] == w)){
            bag.clear();
            bag.insert(G[u].id);
            bag.insert(G[v].id);
            bag.insert(G[x].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[a].id, bag));

            bag.clear();
            bag.insert(G[w].id);
            bag.insert(G[v].id);
            bag.insert(G[x].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[c].id, bag));

            bag.clear();
            bag.insert(G[w].id);
            bag.insert(G[u].id);
            bag.insert(G[x].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[b].id, bag));

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

//checks if there exists a vertex, such that its neighbours induce a clique (Simplicial)
template <typename G_t>
bool Simplicial(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 0)
            continue;

        //N is a clique, if no "edge miss" occures
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
            std::set<unsigned int> bag;
            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++)
                bag.insert(G[*nIt1].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vIt].id, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > (int)bag.size())? low : (int)bag.size();
            return true;
        }
    }
    return false;
}


//checks if there exists a almost simplicial vertex in G (includes the case of a simplicial vertex)
template <typename G_t>
bool AlmostSimplicial(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 0)
            continue;

        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        N.push_back(*vIt);

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++)
            N.push_back(*nIt);

        //N except one vertex now potentially is a clique
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
                if(specialNeighbourFound && (*nIt1 == specialNeighbour || *nIt2 == specialNeighbour))
                    continue;

                if(!boost::edge(*nIt1, *nIt2, G).second){
                    if(specialNeighbourFound){
                        //#special neighbours > 1
                        isAlmostSimplicial = false;
                        goto DOUBLE_BREAK;
                    }
                    //*nIt1 or *nIt2 is a special neighbour
                    cand1 = *nIt1;
                    cand2 = *nIt2;
                    missingEdgesCount++;
                }
            }

            if(missingEdgesCount > 0){
                if(missingEdgesCount == 1){
                    //cand2 has to be the special neighbour
                    specialNeighbour = cand2;
                }
                else{
                    //cand1 has to be the special neighbour
                    specialNeighbour = cand1;
                }
                specialNeighbourFound = true;
            }
        }

        DOUBLE_BREAK:

        if(isAlmostSimplicial){
            //adding the edges, if specialNeighbourFound is true, N is a clique and vertexIt is a simplicial vertex
            if(specialNeighbourFound){
                for(unsigned int i = 0; i < N.size(); i++){
                    if(N.at(i) != specialNeighbour)
                        boost::add_edge(specialNeighbour, N.at(i), G);
                }
            }

            std::set<unsigned int> bag;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++)
                bag.insert(G[*nIt].id);

            bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(G[*vIt].id, bag));

            boost::clear_vertex(*vIt, G);

            low = (low > (int)bag.size())? low : (int)bag.size();
            return true;
        }
    }
    return false;
}

//glues a single bag with the current tree decomposition
template<typename T_t>
void glue_bag_preprocessing(std::set<unsigned int> &bag, unsigned int preprocessed_node, T_t &T){
    if(boost::num_vertices(T) == 0){
        bag.insert(preprocessed_node);
        typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
        T[t_dec_node].bag = bag;

        return;
    }

    typename boost::graph_traits<T_t>::vertex_iterator vertexIt, vertexEnd;

    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(T); vertexIt != vertexEnd; vertexIt++){
        if(std::includes(T[*vertexIt].bag.begin(), T[*vertexIt].bag.end(), bag.begin(), bag.end())){
            bag.insert(preprocessed_node);
            typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
            T[t_dec_node].bag = bag;

            boost::add_edge(*vertexIt, t_dec_node, T);
            return;
        }
    }

    //case for a disconnected graph
    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
    bag.insert(preprocessed_node);
    T[t_dec_node].bag = bag;
    boost::tie(vertexIt, vertexEnd) = boost::vertices(T);
    boost::add_edge(*vertexIt, t_dec_node, T);

}

//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version stores the resulting bags in a vector and does not call further algorithms
template <typename G_t>
void _preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    if(Twig(G, bags, low) || Series(G, bags, low) || Triangle(G, bags, low) || Buddy(G, bags, low) || Cube(G, bags, low) ||
       Simplicial(G, bags, low) || AlmostSimplicial(G, bags, low)){

        _preprocessing(G, bags, low);
    }
    else if (boost::num_edges(G) != 0)
        low = (low > 4)? low : 4;
}

template <typename G_t>
void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, int &low){
    Islet(G, bags, low);
    _preprocessing(G, bags, low);
}

template <typename G_t>
void preprocessing(G_t &G, std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags){
    int low = -1;
    preprocessing(G, bags, low);
}

//glues bags with the current tree decomposition
template<typename T_t>
void preprocessing_glue_bags(std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > &bags, T_t &T){
    for(unsigned int i = bags.size(); i > 0; i--)
        glue_bag_preprocessing(bags[i-1].get<1>(), bags[i-1].get<0>(), T);
}

}

#endif
