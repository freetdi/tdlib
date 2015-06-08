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

#ifndef TD_PARTIAL_TW4
#define TD_PARTIAL_TW4

#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
#include "TD_preprocessing.hpp"

#define DEBUG

namespace treedec{

//YO reduction
template <typename G_t>
bool YO(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,y,z;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all (at most 3) combinations of choosing neighbours y and w of x with y,z of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                y = Nx[i];
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[j], G) != 3)
                        continue;
                    z = Nx[j];

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Ny, Nz;

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Ny.push_back(*nIt);

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(z, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Nz.push_back(*nIt);

                    //y and z must have a vertex 'a' in common, call the neighbour of z, that is not a or x, 'c'
                    //call the neighbour of y that is not a or x, 'b'. There exists 4 possibilities that can occure.
                    typename boost::graph_traits<G_t>::vertex_descriptor a,b,c,d = NULL;

                    if(Ny[0] == Nz[0]){
                        a = Ny[0];
                        b = Ny[1];
                        c = Nz[1];
                    }
                    else if(Ny[0] == Nz[1]){
                        a = Ny[0];
                        b = Ny[1];
                        c = Nz[0];
                    }
                    else if(Ny[1] == Nz[0]){
                        a = Ny[1];
                        b = Ny[0];
                        c = Nz[1];
                    }
                    else if(Ny[1] == Nz[1]){
                        a = Ny[1];
                        b = Ny[0];
                        c = Nz[0];
                    }
                    else
                        continue;

                    //call the neighbour of x, that is not y or z,'d'
                    if(i == 0){
                        if(j == 1)
                            d = Nx[2];
                        else
                            d = Nx[1];
                    }
                    else
                        d = Nx[0];

                    //check for all edges, that must/must not exist
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ab = boost::edge(a,b,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ac = boost::edge(a,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ad = boost::edge(a,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bc = boost::edge(b,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bd = boost::edge(b,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_cd = boost::edge(c,d,G);

                    if(!existsEdge_ab.second && !existsEdge_ac.second && !existsEdge_ad.second && !existsEdge_bc.second && existsEdge_bd.second && existsEdge_cd.second){
                        preprocessed_node = G[x].id;
                        for(unsigned int l = 0; l < Nx.size(); l++)
                            bag.insert(G[Nx[l]].id);

                        boost::add_edge(y,z,G);
                        boost::add_edge(y,d,G);
                        boost::add_edge(z,d,G);
                        boost::clear_vertex(x, G);
                        boost::remove_vertex(x, G);
#ifdef DEBUG
                        std::cout << "YO [x=" << G[x].id << ", y=" << G[y].id << ", z=" << G[z].id << ", a=" << G[a].id << ", b=" << G[b].id << ", c=" << G[c].id << ", d=" << G[d].id << "]" << std::endl;
#endif
                        return true;
                    }
                }
            }   
        }
    }
    return false;
}

template <typename G_t>
bool H7(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,y;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all 3 combinations of choosing y
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                y = Nx[i];

                std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Ny;

                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++){
                    if(*nIt != x)
                        Ny.push_back(*nIt);
                }

                //for both possibilities of choosing c and both possibilities of choosing b: check if the edge {c,b} exists
                //if this edge exists, it is clear, which neighbours of x or y have to be a,b,c,d. 
                typename boost::graph_traits<G_t>::vertex_descriptor a,b,c,d = NULL;
                std::pair<typename G_t::edge_descriptor, bool> existsEdge;
                
                std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx_;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++){
                    if(*nIt != y)
                        Nx_.push_back(*nIt);
                }
 
                if(boost::edge(Ny[0],Nx_[0],G).second){
                    a = Ny[1];
                    b = Nx_[0];                         
                    c = Ny[0];
                    d = Nx_[1];
                }  
                else if(boost::edge(Ny[0],Nx_[1],G).second){
                    a = Ny[1];
                    b = Nx_[1];                         
                    c = Ny[0];
                    d = Nx_[0];
                }
                else if(boost::edge(Ny[1],Nx_[0],G).second){
                    a = Ny[0];
                    b = Nx_[0];                         
                    c = Ny[1];
                    d = Nx_[1];
                }
                else if(boost::edge(Ny[1],Nx_[1],G).second){
                    a = Ny[0];
                    b = Nx_[1];                         
                    c = Ny[1];
                    d = Nx_[0];
                }
                else
                    continue;
                

                //check for all edges, that must/must not exist, {b,c} already exists
                std::pair<typename G_t::edge_descriptor, bool> existsEdge_ab = boost::edge(a,b,G);
                std::pair<typename G_t::edge_descriptor, bool> existsEdge_ac = boost::edge(a,c,G);
                std::pair<typename G_t::edge_descriptor, bool> existsEdge_ad = boost::edge(a,d,G);
                std::pair<typename G_t::edge_descriptor, bool> existsEdge_bd = boost::edge(b,d,G);
                std::pair<typename G_t::edge_descriptor, bool> existsEdge_cd = boost::edge(c,d,G);
                
                if(existsEdge_ab.second && !existsEdge_ac.second && !existsEdge_ad.second && !existsEdge_bd.second && !existsEdge_cd.second){
                    preprocessed_node = G[x].id;
                    for(unsigned int l = 0; l < Nx.size(); l++)
                        bag.insert(G[Nx[l]].id);

                    boost::add_edge(y,b,G);	
                    boost::add_edge(y,d,G);
                    boost::add_edge(b,d,G);

                    boost::clear_vertex(x, G);
                    boost::remove_vertex(x, G);
#ifdef DEBUG
                    std::cout << "H7 [x=" << G[x].id << ", y=" << G[y].id << ", a=" << G[a].id << ", b=" << G[b].id << ", c=" << G[c].id << ", d=" << G[d].id << "]" << std::endl;
#endif
                    return true;
                }
            }
        }
    }
    return false;
}

template <typename G_t>
bool TO(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //all neighbours of x must have degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    goto RTN_LOOP;
            }

            typename boost::graph_traits<G_t>::vertex_descriptor y, w, z;

            //for all (at most 3) combinations of choosing neighbours y and w of x with y, w of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                y = Nx[i];
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[j], G) != 3)
                        continue;
                    w = Nx[j];

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Ny, Nw;

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Ny.push_back(*nIt);

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Nw.push_back(*nIt);


                    //y and w must have a vertex 'c' in common, call the neighbour of y, that is not c or x, 'a'
                    //call the neighbour of w that is not c or x, 'd'. There exists 4 possibilities that can occure.
                    typename boost::graph_traits<G_t>::vertex_descriptor a, b, c, d;

                    if(Ny[0] == Nw[0]){
                        a = Ny[1];
                        c = Ny[0];
                        d = Nw[1];
                    }
                    else if(Ny[0] == Nw[1]){
                        a = Ny[1];
                        c = Ny[0];
                        d = Nw[0];
                    }
                    else if(Ny[1] == Nw[0]){
                        a = Ny[0];
                        c = Ny[1];
                        d = Nw[1];
                    }
                    else if(Ny[1] == Nw[1]){
                        a = Ny[0];
                        c = Ny[1];
                        d = Nw[0];
                    }
                    else
                        continue;

                    //call the neighbour of x, that is not y or w, 'z'
                    if(i == 0){
                        if(j == 1)
                            z = Nx[2];
                        else
                            z = Nx[1];
                    }
                    else
                        z = Nx[0];

                    //call the neighbour of z, that is not x or d, 'b'
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(z, G); nIt != nEnd; nIt++){
                        if(*nIt != x && *nIt != d){
                            b = *nIt;
                            break;
                        }
                    }

                    //check for all edges, that must/must not exist
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ab = boost::edge(a,b,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ac = boost::edge(a,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ad = boost::edge(a,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bc = boost::edge(b,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bd = boost::edge(b,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_cd = boost::edge(c,d,G);
                
                    if(existsEdge_ab.second && !existsEdge_ac.second && !existsEdge_ad.second && !existsEdge_bc.second && !existsEdge_bd.second && !existsEdge_cd.second){
                        preprocessed_node = G[x].id;
                        for(unsigned int l = 0; l < Nx.size(); l++)
                            bag.insert(G[Nx[l]].id);

                        boost::add_edge(y,w,G);	
                        boost::add_edge(y,z,G);
                        boost::add_edge(w,z,G);

                        boost::clear_vertex(x, G);
                        boost::remove_vertex(x, G);
#ifdef DEBUG
                        std::cout << "TO [x=" << G[x].id << ", y=" << G[y].id << ", z=" << G[z].id << ", w=" << G[w].id << ", a=" << G[a].id << ", b=" << G[b].id << ", c=" << G[c].id << ", d=" << G[d].id << "]" << std::endl;
#endif
                        return true;
                    }
                }
            }
        }
    RTN_LOOP:
        continue;
    }
    return false;
}


template <typename G_t>
bool YI(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,y,z;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all (at most 3) combinations of choosing neighbours y and w of x with y,z of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                y = Nx[i];
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[j], G) != 3)
                        continue;
                    z = Nx[j];

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Ny, Nz;

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Ny.push_back(*nIt);

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(z, G); nIt != nEnd; nIt++)
                        if(*nIt != x)
                            Nz.push_back(*nIt);

                    //y and z must have a vertex 'a' in common, call the neighbour of z, that is not a or x, 'c'
                    //call the neighbour of y that is not a or x, 'b'. There exists 4 possibilities that can occure.
                    typename boost::graph_traits<G_t>::vertex_descriptor a,b,c,d = NULL;

                    if(Ny[0] == Nz[0]){
                        a = Ny[0];
                        b = Ny[1];
                        c = Nz[1];
                    }
                    else if(Ny[0] == Nz[1]){
                        a = Ny[0];
                        b = Ny[1];
                        c = Nz[0];
                    }
                    else if(Ny[1] == Nz[0]){
                        a = Ny[1];
                        b = Ny[0];
                        c = Nz[1];
                    }
                    else if(Ny[1] == Nz[1]){
                        a = Ny[1];
                        b = Ny[0];
                        c = Nz[0];
                    }
                    else
                        continue;

                    //call the neighbour of x, that is not y or z,'d'
                    if(i == 0){
                        if(j == 1)
                            d = Nx[2];
                        else
                            d = Nx[1];
                    }
                    else
                        d = Nx[0];

                    //check for all edges, that must/must not exist
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ab = boost::edge(a,b,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ac = boost::edge(a,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_ad = boost::edge(a,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bc = boost::edge(b,c,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_bd = boost::edge(b,d,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_cd = boost::edge(c,d,G);

                    if(!existsEdge_ab.second && !existsEdge_ac.second && !existsEdge_ad.second && existsEdge_bc.second && !existsEdge_bd.second && !existsEdge_cd.second){
                        preprocessed_node = G[x].id;
                        for(unsigned int l = 0; l < Nx.size(); l++)
                            bag.insert(G[Nx[l]].id);

                        boost::add_edge(y,z,G);
                        boost::add_edge(y,d,G);
                        boost::add_edge(z,d,G);
                        boost::clear_vertex(x, G);
                        boost::remove_vertex(x, G);
#ifdef DEBUG
                        std::cout << "YI [x=" << G[x].id << ", y=" << G[y].id << ", z=" << G[z].id << ", a=" << G[a].id << ", b=" << G[b].id << ", c=" << G[c].id << ", d=" << G[d].id << "]" << std::endl;
#endif
                        return true;
                    }
                }
            }   
        }
    }
    return false;
}

template <typename G_t>
bool L1(G_t &G, std::set<unsigned int> &bag1, std::set<unsigned int> &bag2, unsigned int &preprocessed_node1, unsigned int &preprocessed_node2){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,v3,v4;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all (at most 3) combinations of choosing neighbours v3 and v4 of x with v3, v4 of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                v3 = Nx[i];
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[j], G) != 3)
                        continue;
                    v4 = Nx[j];

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv3, Nv4;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v3, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv3.push_back(*nIt);
                    }
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v4, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv4.push_back(*nIt);
                    }

                    typename boost::graph_traits<G_t>::vertex_descriptor v1,v2,v5,v6,y = NULL;

                    if(Nv3[0] == Nv4[0]){
                        y = Nv3[0];
                        v6 = Nv4[1];
                        v1 = Nv3[1];
                    }
                    else if(Nv3[0] == Nv4[1]){
                        y = Nv3[0];
                        v6 = Nv4[0];
                        v1 = Nv3[1];
                    }
                    else if(Nv3[1] == Nv4[0]){
                        y = Nv3[1];
                        v6 = Nv4[1];
                        v1 = Nv3[0];
                    }
                    else if(Nv3[1] == Nv4[1]){
                        y = Nv3[1];
                        v6 = Nv4[0];
                        v1 = Nv3[0];
                    }
                    else
                        continue;

                    if(boost::out_degree(y, G) != 3)
                        continue;

                    //call the neighbour of x, that is not v3 or v4,'v2'
                    if(i == 0){
                        if(j == 1)
                            v2 = Nx[2];
                        else
                            v2 = Nx[1];
                    }
                    else
                        v2 = Nx[0];

                    //call the neighbour of y that is not v3 or v4, 'v5' 
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++){
                        if(*nIt != v3 && *nIt != v4){
                            v5 = *nIt;
                            break;
                        }
                    }

                    //check for all edges, that must/must not exist
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_12 = boost::edge(v1,v2,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_15 = boost::edge(v1,v5,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_16 = boost::edge(v1,v6,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_25 = boost::edge(v2,v5,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_26 = boost::edge(v2,v6,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_56 = boost::edge(v5,v6,G);

                    if(existsEdge_12.second && !existsEdge_15.second && !existsEdge_16.second && !existsEdge_25.second && !existsEdge_26.second && existsEdge_56.second){
                        preprocessed_node1 = G[x].id;
                        preprocessed_node2 = G[y].id;
                        for(unsigned int l = 0; l < Nx.size(); l++)
                            bag1.insert(G[Nx[l]].id);
                        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(y, G); nIt != nEnd; nIt++)
                            bag2.insert(G[*nIt].id);

                        boost::add_edge(v2,v3,G);
                        boost::add_edge(v2,v4,G);
                        boost::add_edge(v3,v4,G);
                        boost::clear_vertex(x, G);
                        boost::remove_vertex(x, G);

                        boost::add_edge(v3,v5,G);
                        boost::add_edge(v4,v5,G);
                        boost::clear_vertex(y, G);
                        boost::remove_vertex(y, G);
#ifdef DEBUG
                        std::cout << "L1 [x=" << G[x].id << ", y=" << G[y].id << ", v1=" << G[v1].id << ", v2=" << G[v2].id << ", v3=" << G[v3].id << ", v4=" << G[v4].id << ", v5=" << G[v5].id << ", v6=" << G[v6].id << "]" << std::endl;
#endif
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

template <typename G_t>
bool L2(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,v3,v4;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all (at most 3) combinations of choosing neighbours v3 and v4 of x with v3 of degree 4, v3 of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[i], G) == 3 && boost::out_degree(Nx[j], G) == 4){
                        v3 = Nx[j];
                        v4 = Nx[i];
                    }
                    else if(boost::out_degree(Nx[j], G) == 3 && boost::out_degree(Nx[i], G) == 4){
                        v3 = Nx[i];
                        v4 = Nx[j];
                    }
                    else
                        continue;

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv3, Nv4;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v3, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv3.push_back(*nIt);
                    }
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v4, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv4.push_back(*nIt);
                    }

                    typename boost::graph_traits<G_t>::vertex_descriptor v0,v1,v2,v5,v6 = NULL;

                    if(Nv3[0] == Nv4[0]){
                        v0 = Nv3[1];
                        v1 = Nv3[2];
                        v2 = Nv3[0];
                        v5 = Nv4[1];
                    }
                    else if(Nv3[0] == Nv4[1]){
                        v0 = Nv3[1];
                        v1 = Nv3[2];
                        v2 = Nv3[0];
                        v5 = Nv4[0];
                    }
                    else if(Nv3[1] == Nv4[0]){
                        v0 = Nv3[0];
                        v1 = Nv3[2];
                        v2 = Nv3[1];
                        v5 = Nv4[1];
                    }
                    else if(Nv3[1] == Nv4[1]){
                        v0 = Nv3[0];
                        v1 = Nv3[2];
                        v2 = Nv3[1];
                        v5 = Nv4[0];
                    }
                    else if(Nv3[2] == Nv4[0]){
                        v0 = Nv3[0];
                        v1 = Nv3[1];
                        v2 = Nv3[2];
                        v5 = Nv4[1];
                    }
                    else if(Nv3[2] == Nv4[1]){
                        v0 = Nv3[0];
                        v1 = Nv3[1];
                        v2 = Nv3[2];
                        v5 = Nv4[0];
                    }
                    else
                        continue;

                    if(boost::out_degree(v2, G) != 4)
                        continue;

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv3_, Nv2;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v3, G); nIt != nEnd; nIt++){
                        if(*nIt != x && *nIt != v2)
                            Nv3_.push_back(*nIt);
                    }
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v2, G); nIt != nEnd; nIt++){
                        if(*nIt != v3 && *nIt != v4)
                            Nv2.push_back(*nIt);
                    }

                    //check, if (the neighbours of v3 without x and v2) == (the neighbours of v2 without v3 and v4) 
                    if(!((Nv3_[0] == Nv2[0] && Nv3_[1] == Nv2[1]) || (Nv3_[0] == Nv2[1] && Nv3_[1] == Nv2[0])))
                        continue;

                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++){
                        if(*nIt != v3 && *nIt != v4){
                            v6 = *nIt;
                            break;
                        }
                    }

                    //check for all edges, that must/must not exist
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_01 = boost::edge(v0,v1,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_05 = boost::edge(v0,v5,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_06 = boost::edge(v0,v6,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_15 = boost::edge(v1,v5,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_16 = boost::edge(v1,v6,G);
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_56 = boost::edge(v5,v6,G);

                    if(!existsEdge_01.second && !existsEdge_05.second && !existsEdge_06.second && !existsEdge_15.second && !existsEdge_16.second && existsEdge_56.second){
                        preprocessed_node = G[x].id;
                        for(unsigned int l = 0; l < Nx.size(); l++)
                            bag.insert(G[Nx[l]].id);

                        boost::add_edge(v3,v4,G);
                        boost::add_edge(v3,v6,G);
                        boost::add_edge(v4,v6,G);
                        boost::clear_vertex(x, G);
                        boost::remove_vertex(x, G);
#ifdef DEBUG
                        std::cout << "L2 [x=" << G[x].id << "v0=" << G[v0].id << ", v1=" << G[v1].id << ", v2=" << G[v2].id << ", v3=" << G[v3].id << ", v4=" << G[v4].id << ", v5=" << G[v5].id << ", v6=" << G[v6].id << "]" << std::endl;
#endif
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


template <typename G_t>
bool L3(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex x of degree 3
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 3){
            typename boost::graph_traits<G_t>::vertex_descriptor x,v3,v4;
            x = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nx;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++)
                Nx.push_back(*nIt);

            //for all (at most 3) combinations of choosing neighbours v3 and v4 of x with v3, v4 of degree 3
            for(unsigned int i = 0; i < Nx.size(); i++){
                if(boost::out_degree(Nx[i], G) != 3)
                    continue;
                v3 = Nx[i];
                for(unsigned int j = i+1; j < Nx.size(); j++){
                    if(boost::out_degree(Nx[j], G) != 3)
                        continue;
                    v4 = Nx[j];

                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv3, Nv4;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v3, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv3.push_back(*nIt);
                    }
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v4, G); nIt != nEnd; nIt++){
                        if(*nIt != x)
                            Nv4.push_back(*nIt);
                    }

                    typename boost::graph_traits<G_t>::vertex_descriptor v0,v1,v2,v5,v6;

                    for(unsigned int k = 0; k < Nv4.size(); k++){
                        if(boost::out_degree(Nv4[k], G) != 3)
                            continue;
                        v2 = Nv4[k];

                        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv2;
                        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v2, G); nIt != nEnd; nIt++){
                            if(*nIt != v4)
                                Nv2.push_back(*nIt);
                        }

                        //check, if (the neighbours of v3 without x) == (the neighbours of v2 without v4) 
                        if(!((Nv3[0] == Nv2[0] && Nv3[1] == Nv2[1]) || (Nv3[0] == Nv2[1] && Nv3[1] == Nv2[0])))
                            continue;

                        v0 = Nv3[0];
                        v1 = Nv3[1];

                        if(k == 0)
                            v5 = Nv4[1];
                        else
                            v5 = Nv4[0];

                        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++){
                            if(*nIt != v3 && *nIt != v4){
                                v6 = *nIt;
                                break;
                            }
                        }

                        //check for all edges, that must/must not exist
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_01 = boost::edge(v0,v1,G);
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_05 = boost::edge(v0,v5,G);
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_06 = boost::edge(v0,v6,G);
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_15 = boost::edge(v1,v5,G);
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_16 = boost::edge(v1,v6,G);
                        std::pair<typename G_t::edge_descriptor, bool> existsEdge_56 = boost::edge(v5,v6,G);

                        if(!existsEdge_01.second && !existsEdge_05.second && !existsEdge_06.second && !existsEdge_15.second && !existsEdge_16.second && existsEdge_56.second){
                            preprocessed_node = G[x].id;
                            for(unsigned int l = 0; l < Nx.size(); l++)
                                bag.insert(G[Nx[l]].id);

                            boost::add_edge(v3,v4,G);
                            boost::add_edge(v3,v6,G);
                            boost::add_edge(v4,v6,G);
                            boost::clear_vertex(x, G);
                            boost::remove_vertex(x, G);
#ifdef DEBUG
                            std::cout << "L3 [x=" << G[x].id << "v0=" << G[v0].id << ", v1=" << G[v1].id << ", v2=" << G[v2].id << ", v3=" << G[v3].id << ", v4=" << G[v4].id << ", v5=" << G[v5].id << ", v6=" << G[v6].id << "]" << std::endl;
#endif
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

template <typename G_t>
bool L4(G_t &G, std::set<unsigned int> &bag, unsigned int &preprocessed_node){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd; 
    //find a vertex v4 of degree 4
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 4){
            typename boost::graph_traits<G_t>::vertex_descriptor v2,v3,v4,x;
            v4 = *vIt;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv4;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v4, G); nIt != nEnd; nIt++)
                Nv4.push_back(*nIt);

            //for all (at most 6) combinations of choosing neighbours v3 and v4 of x with v3, v4 of degree 4
            for(unsigned int i = 0; i < Nv4.size(); i++){
                if(boost::out_degree(Nv4[i], G) != 4)
                    continue;
                v2 = Nv4[i];

                for(unsigned int j = i+1; j < Nv4.size(); j++){
                    if(boost::out_degree(Nv4[j], G) != 4)
                        continue;
                    v3 = Nv4[j];

                    std::pair<typename G_t::edge_descriptor, bool> existsEdge_34 = boost::edge(v3,v4,G);

                    if(!existsEdge_34.second)
                        continue;

                    //for the both choices of choosing x with x of degree 3
                    for(unsigned int k = 0; k < Nv4.size(); k++){
                        if(boost::out_degree(Nv4[k], G) == 3 && k != i && k != j){
                            x = Nv4[k]; 

                            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> Nv3, Nx;
                            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v3, G); nIt != nEnd; nIt++){
                                if(*nIt != v2 && *nIt != v4)
                                    Nv3.push_back(*nIt);
                            }
                            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(x, G); nIt != nEnd; nIt++){
                                if(*nIt != v4)
                                    Nx.push_back(*nIt);
                            }

                            typename boost::graph_traits<G_t>::vertex_descriptor v0,v1,v5,v6 = NULL;


                            //v3 and x must have a vertex v5 in common. With v5, we know v1 and v6 
                            if(Nv3[0] == Nx[0]){
                                v5 = Nv3[0];
                                v1 = Nv3[1];
                                v6 = Nx[1];
                            }
                            else if(Nv3[0] == Nx[1]){
                                v5 = Nv3[0];
                                v1 = Nv3[1];
                                v6 = Nx[0];
                            }
                            else if(Nv3[1] == Nx[0]){
                                v5 = Nv3[1];
                                v1 = Nv3[0];
                                v6 = Nx[1];
                            }
                            else if(Nv3[1] == Nx[1]){
                                v5 = Nv3[1];
                                v1 = Nv3[0];
                                v6 = Nx[0];
                            }
                            else
                                continue;

                            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v2, G); nIt != nEnd; nIt++){
                                if(*nIt != v1 && *nIt != v3 && *nIt != v4){
                                    v0 = *nIt;
                                    break;
                                }
                            }

                            //check for the remaining edges
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_01 = boost::edge(v0,v1,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_02 = boost::edge(v0,v2,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_04 = boost::edge(v0,v4,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_05 = boost::edge(v0,v5,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_06 = boost::edge(v0,v6,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_12 = boost::edge(v1,v2,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_15 = boost::edge(v1,v5,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_16 = boost::edge(v1,v6,G);
                            std::pair<typename G_t::edge_descriptor, bool> existsEdge_56 = boost::edge(v5,v6,G);

                            if(!existsEdge_01.second && existsEdge_02.second && existsEdge_04.second && !existsEdge_05.second && !existsEdge_06.second
                            && existsEdge_12.second && !existsEdge_15.second && !existsEdge_16.second && !existsEdge_56.second){
                                preprocessed_node = G[x].id;
                                for(unsigned int l = 0; l < Nx.size(); l++)
                                    bag.insert(G[Nx[l]].id);

                                boost::add_edge(v4,v5,G);
                                boost::add_edge(v4,v6,G);
                                boost::add_edge(v5,v6,G);
                                boost::clear_vertex(x, G);
                                boost::remove_vertex(x, G);
#ifdef DEBUG
                                std::cout << "L4 [" << "x=" << G[x].id << "v0=" << G[v0].id << ", v1=" << G[v1].id << ", v2=" << G[v2].id << ", v3=" << G[v3].id << ", v4=" << G[v4].id << ", v5=" << G[v5].id << ", v6=" << G[v6].id << "]" << std::endl;
#endif
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}



template <typename G_t, typename T_t>
void partial_tw4_preprocessing(G_t &G, T_t &T){
    std::set<unsigned int> bag, bag2;
    unsigned int preprocessed_node, preprocessing_node2;
    int low = -1;
    std::cout << "|V|: " << boost::num_vertices(G) << std::endl;

    if(Islet(G, bag, preprocessed_node, low) || Twig(G, bag, preprocessed_node, low) || Series(G, bag, preprocessed_node, low)
    || Triangle(G, bag, preprocessed_node, low) || Buddy(G, bag, preprocessed_node, low) || Cube(G, bag, preprocessed_node, low)
    || YO(G, bag, preprocessed_node) || H7(G, bag, preprocessed_node) || TO(G, bag, preprocessed_node) || YI(G, bag, preprocessed_node) 
    || L1(G, bag, bag2, preprocessed_node, preprocessing_node2) || L2(G, bag, preprocessed_node) || L3(G, bag, preprocessed_node) || L4(G, bag, preprocessed_node)){

        partial_tw4_preprocessing(G, T);
    }
    else{
        if(boost::num_vertices(G) == 0)
            std::cout << "graph has been fully tw4-preprocessed" << std::endl;
        else
            std::cout << "no tw4-rule applicable" << std::endl;
    }
}

}

#endif
