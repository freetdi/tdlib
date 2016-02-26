// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
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
// Offers functionality to compute tree decompositions of graphs
// according to various heuristic and functions, that convert
// tree decompositions to elimination orderings and vice versa.
// Also the LEX-M algorithm is included in this header
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
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;
//
//
// These functions are most likely to be interesting for outside use:
//
// void minDegree_decomp(G_t G&, T_t &T)
// void fillIn_decomp(G_t G&, T_t &T)
// void minDegree_ordering(G_t G, std::vector<unsigned int> &elim_ordering)
// void fillIn_ordering(G_t G, std::vector<unsigned int> &elim_ordering)
// void ordering_to_treedec(G_t G, std::vector<unsigned int> &elimination_ordering, T_t &T)
// void treedec_to_ordering(T_t T, std::vector<unsigned int> &elimination_ordering)
// void LEX_M_minimal_ordering(G_t &G, std::vector<unsigned int> &alpha)
//

#ifndef TD_ELIMINATION_ORDERING
#define TD_ELIMINATION_ORDERING

#include <cmath>
#include <climits>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "TD_preprocessing.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_misc.hpp"

namespace treedec{

namespace detail{

    // DRAFT. no useful interface
template<typename G>
struct degree_mod : public noboost::vertex_callback<G>{
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    degree_mod(misc::DEGS<G>* d, G* g) : _degs(d), _g(g){}

    // reinsert with degree-1
    void operator()(vertex_descriptor v){ untested();
        unsigned deg=boost::degree(v,*_g);
        assert(deg);
        bool done=(*_degs)[deg-1].insert(v).second;
        assert(done);
    }
    //private: not yet.
//    BUG:: hardcoded type
        misc::DEGS<G>* _degs;
    private:
        degree_mod(const degree_mod&){}
        G* _g;
};

}

//Constructs a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic. Ignores isolated vertices.
template <typename G_t, typename T_t>
void _minDegree_decomp(G_t &G, T_t &T){
    std::vector<std::set<unsigned int> > bags(boost::num_vertices(G));
    std::vector<unsigned int> elim_vertices(boost::num_vertices(G));
    misc::DEGS<G_t> degs(G);
    detail::degree_mod<G_t> cb(&degs, &G);

    unsigned int i = 0;
    unsigned min_ntd = 1; // minimum nontrivial vertex degree
    while(boost::num_edges(G) > 0){
        //Search a minimum degree vertex.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        boost::tie(vIt, vEnd) = boost::vertices(G);
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;
        unsigned min_degree = boost::num_vertices(G);
        unsigned num_vert = boost::num_vertices(G);
        auto mdvi = degs[0].begin();

        // recompute ntd can only increase from here
        while(degs[min_ntd].empty()){
            ++min_ntd;
            // min_ntd==num_vert contradicts the outer loop condition
            // (this loop should be safe)
            assert(min_ntd != num_vert);
        }
        assert(!degs[min_ntd].empty());
        assert(min_ntd);
        assert(min_ntd==1 || degs[min_ntd-1].empty());

#ifndef NDEBUG // remove later.
        for(vIt=boost::vertices(G).first; vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree < min_degree && degree > 0){
                min_degree = degree;
            }
        }
        assert(min_degree == min_ntd);
#endif
        mdvi = degs[min_ntd].begin(); // min degree vertex iterator

        std::set<unsigned int> bag;

        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*mdvi, G);
                nIt1 != nEnd; nIt1++){ untested();
            unsigned id=noboost::get_id(G, *nIt1);
            bag.insert(id); // collect neighborhood, inefficient.
        }
        degs.check();
        misc::make_clique(boost::adjacent_vertices(*mdvi, G), G, &cb);
//         degs.check(); not here.
//         degrees of neighbors are wrong

        bags[i] = bag; // std::move?!
        unsigned id=noboost::get_id(G, *mdvi);
        elim_vertices[i++] = id;

        unsigned n=degs[min_ntd].erase(*mdvi);
//        degs[min_ntd].erase(mdvi);
        boost::clear_vertex(*mdvi, G);
        if(min_ntd>1){
            --min_ntd;
        }
        assert(boost::degree(*mdvi, G)==0);
        
        assert(n==1);

#ifndef NDEBUG //redundant, for checking only
        bool done=degs[0].insert(*mdvi).second;
        assert(done);
#endif
        degs.check(); // g has a node of deg 445 that degs doesn't know about.
    }

    for(; i > 0; i--){
        glue_bag(bags[i-1], elim_vertices[i-1], T);
    }
}

//Constructs a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic.
template <typename G_t, typename T_t>
void minDegree_decomp(G_t &G, T_t &T){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    int low = 0;
    treedec::Islet(G, bags, low);
    treedec::_minDegree_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}


//Constructs a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic. Ignores isolated vertices.
template <typename G_t, typename T_t>
void _fillIn_decomp(G_t &G, T_t &T){
    std::vector<std::set<unsigned int> > bags(boost::num_vertices(G));
    std::vector<unsigned int> elim_vertices(boost::num_vertices(G));

    unsigned int i = 0;
    while(boost::num_edges(G) > 0){
        //Search a vertex v such that least edges are missing for making the neighbourhood of v a clique.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;

        boost::tie(vIt, vEnd) = boost::vertices(G);

        unsigned int min_fill = UINT_MAX;
        for(; vIt != vEnd; vIt++){
            if(boost::out_degree(*vIt, G) == 0){
                continue;
            }

            unsigned int current_fill = 0;

            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
                nIt2 = nIt1;
                nIt2++;
                for(; nIt2 != nEnd; nIt2++){
                    if(!boost::edge(*nIt1, *nIt2, G).second){
                        current_fill++;
                    }
                }
            }

            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0){
                    break;
                }
            }
        }

        //Make the neighbourhood of 'min_vertex' a clique.
        std::set<unsigned int> bag;

        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd; nIt1++){
            unsigned id=noboost::get_id(G, *nIt1);
            bag.insert(id);
            nIt2 = nIt1;
            nIt2++;
            for(; nIt2 != nEnd; nIt2++){
                boost::add_edge(*nIt1, *nIt2, G);
            }
        }

        bags[i] = bag;
        unsigned id=noboost::get_id(G, min_vertex);
        elim_vertices[i++] = id;

        boost::clear_vertex(min_vertex, G);
    }

    for(; i > 0; i--){
        glue_bag(bags[i-1], elim_vertices[i-1], T);
    }
}

//Constructs a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic.
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    int low = 0;
    treedec::Islet(G, bags, low);
    treedec::_fillIn_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}


//Computes an elimination ordering according to the minDegree heuristic (version used for postprocessing algorithms).
template<typename G_t>
void _minDegree_ordering(G_t G, std::vector<unsigned int> &elim_ordering, std::vector<bool> &visited){
    while(true){
        //Search a minimum degree vertex.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_degree = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned id=noboost::get_id(G,*vIt);
            if(visited[id]){
                continue;
            }

            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree < min_degree){
                min_degree  = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree == UINT_MAX){
            return;
        }

        //Make the neighbourhood of 'min_vertex' a clique.
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd; nIt1++){
            nIt2 = nIt1;
            nIt2++;
            for(; nIt2 != nEnd; nIt2++){
                boost::add_edge(*nIt1, *nIt2, G);
            }
        }

        unsigned id=noboost::get_id(G, min_vertex);
        elim_ordering.push_back(id);
        visited[id] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

//Computes an elimination ordering according to minDegree heuristic (version used for postprocessing algorithms).
template<typename G_t>
void minDegree_ordering(G_t G, std::vector<unsigned int> &elim_ordering){
    std::vector<bool> visited(boost::num_vertices(G), false);
    _minDegree_ordering(G, elim_ordering, visited);
}


//Computes an elimination ordering according to fillIn heuristic (version used for postprocessing algorithms).
template<typename G_t>
void _fillIn_ordering(G_t G, std::vector<unsigned int> &elim_ordering, std::vector<bool> &visited){
    while(true){
        //Search a vertex v such that least edges are missing for making the neighbourhood of v a clique.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_fill = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned id=noboost::get_id(G, *vIt);
            if(visited[id]){
                continue;
            }

            if(boost::out_degree(*vIt, G) == 0){
                min_vertex = *vIt;
                min_fill = 0;
                break;
            }

            unsigned int current_fill = 0;

            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
                nIt2 = nIt1;
                nIt2++;
                for(; nIt2 != nEnd; nIt2++){
                    if(!boost::edge(*nIt1, *nIt2, G).second){
                        current_fill++;
                    }
                }
            }

            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0){
                    break;
                }
            }
        }

        if(min_fill == UINT_MAX){
            return;
        }

        //Make the neighbourhood of 'min_vertex' a clique.
        for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd; nIt1++){
            nIt2 = nIt1;
            nIt2++;
            for(; nIt2 != nEnd; nIt2++){
                boost::add_edge(*nIt1, *nIt2, G);
            }
        }

        unsigned id=noboost::get_id(G, min_vertex);
        elim_ordering.push_back(id);
        visited[id] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

//Computes an elimination ordering according to fillIn heuristic (version used for postprocessing algorithms).
template<typename G_t>
void fillIn_ordering(G_t &G, std::vector<unsigned int> &elim_ordering){
    std::vector<bool> visited(boost::num_vertices(G)+1, false);
    _fillIn_ordering(G, elim_ordering, visited);
}


//Make G a filled graph according to the provided elimination_ordering. Stores the cliques in C and the additional
//edges in F.
template <typename G_t>
void make_filled_graph(G_t &G, std::vector<unsigned int> &elim_ordering, std::vector<std::set<unsigned int> > &C, std::vector<std::vector<std::vector<unsigned int> > > &F){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    C.resize(elim_ordering.size());
    F.resize(elim_ordering.size());

    std::vector<bool> visited(boost::num_vertices(G), false);

    for(unsigned int i = 0; i < elim_ordering.size(); i++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        std::set<unsigned int> N_i, E_i;
        C[i].insert(elim_ordering[i]);
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[elim_ordering[i]], G); nIt != nEnd; nIt++){
            unsigned id=noboost::get_id(G, *nIt);
            if(!visited[id]){
                C[i].insert(id);
            }
        }

        for(std::set<unsigned int>::iterator sIt1 = C[i].begin(); sIt1 != C[i].end(); sIt1++){
            std::set<unsigned int>::iterator sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != C[i].end(); sIt2++){
                if(!boost::edge(idxMap[*sIt1], idxMap[*sIt2], G).second){
                    std::vector<unsigned int> edge;
                    edge.push_back(*sIt1);
                    edge.push_back(*sIt2);
                    F[i].push_back(edge);
                    boost::add_edge(idxMap[*sIt1], idxMap[*sIt2], G);
                }
            }
        }
        visited[elim_ordering[i]] = true;
    }
}

template <typename G_t, typename T_t>
void _ordering_to_treedec(G_t &G, std::vector<unsigned int> &elimination_ordering, T_t &T, unsigned int idx){
    if(idx == elimination_ordering.size())
        return;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);
    typename boost::graph_traits<G_t>::vertex_descriptor elim_vertex = *vIt++;
    for(; vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        if(id == elimination_ordering[idx]){
            elim_vertex = *vIt;
            break;
        }
    }

    //Collect the neighbours of 'elim_vertex'.
    std::set<unsigned int> bag;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(elim_vertex, G); nIt != nEnd; nIt++){
        unsigned id=noboost::get_id(G, *nIt);
        bag.insert(id);
        neighbours.push_back(*nIt);
    }

    //Make the neighbourhood of 'min_vertex' a clique.
    for(unsigned int i = 0; i < neighbours.size(); i++){
        for(unsigned int j = i+1; j < neighbours.size(); j++){
            boost::add_edge(neighbours[i], neighbours[j], G);
        }
    }

    unsigned id=noboost::get_id(G, elim_vertex);
    unsigned int elim_vertex_id = id;

    boost::clear_vertex(elim_vertex, G);

    _ordering_to_treedec(G, elimination_ordering, T, idx+1);
    glue_bag(bag, elim_vertex_id, T);
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t G, std::vector<unsigned int> &elimination_ordering, T_t &T){
    _ordering_to_treedec(G, elimination_ordering, T, 0);
}

template <typename G_t, typename T_t>
void _ordering_to_treedec(G_t &G, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, T_t &T, unsigned int idx){
    std::vector<std::set<unsigned int> > bags;
    std::vector<unsigned int> elim_vertices;

    for(unsigned int i = 0; i < elimination_ordering.size(); i++){
        //Collect the neighbours of elimination vertex i.
        std::set<unsigned int> bag;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(elimination_ordering[i], G); nIt != nEnd; nIt++){
            unsigned id=noboost::get_id(G, *nIt);
            bag.insert(id);
            neighbours.push_back(*nIt);
        }

        //Make the neighbourhood of elimination vertex i a clique.
        for(unsigned int j = 0; j < neighbours.size(); j++){
            for(unsigned int k = j+1; k < neighbours.size(); k++){
                boost::add_edge(neighbours[j], neighbours[k], G);
            }
        }

        boost::clear_vertex(elimination_ordering[i], G);

        bags.push_back(bag);
        unsigned id=noboost::get_id(G, elimination_ordering[i]);
        elim_vertices.push_back(id);

    }

    for(unsigned int i = bags.size(); i > 0; i--){
        glue_bag(bags[i-1], elim_vertices[i-1], T);
    }
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t G, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, T_t &T){
    _ordering_to_treedec(G, elimination_ordering, T, 0);
}

template <typename T_t>
void _treedec_to_ordering(T_t &T, std::vector<unsigned int> &elimination_ordering){
    bool leaf_found = false;

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor leaf, parent;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::out_degree(*tIt, T) <= 1 && !T[*tIt].bag.empty()){
            leaf = *tIt;
            leaf_found = true;
            break;
        }
    }

    if(leaf_found){
        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(leaf, T);
        parent = *nIt;

        std::set<unsigned int> difference;

        if(boost::out_degree(leaf, T) == 1){
            if(!std::includes(T[parent].bag.begin(), T[parent].bag.end(),
                              T[leaf].bag.begin(), T[leaf].bag.end())){
                std::set_difference(T[leaf].bag.begin(), T[leaf].bag.end(), T[parent].bag.begin(),
                                    T[parent].bag.end(), std::inserter(difference, difference.begin()));
            }
            boost::clear_vertex(leaf, T);
        }
        else{
            difference = T[leaf].bag;
        }

        for(std::set<unsigned int>::iterator sIt = difference.begin(); sIt != difference.end(); sIt++){
            elimination_ordering.push_back(*sIt);
        }

        T[leaf].bag.clear();

        _treedec_to_ordering(T, elimination_ordering);
    }
}

template <typename T_t>
void treedec_to_ordering(T_t T, std::vector<unsigned int> &elimination_ordering){
    if(boost::num_vertices(T) == 0){
        return;
    }

    _treedec_to_ordering(T, elimination_ordering);
}

template <typename G_t>
void LEX_M_fill_in(G_t &G, std::vector<std::vector<unsigned int> > &fill_in_edges){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    unsigned int max = idxMap.size();
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<unsigned int> alpha_inv(max+1);
    std::vector<std::vector<unsigned int> > reached_i(max+1);

    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        label[id] = 1.0;
        alpha_inv[i++] = 0;
        visited[id] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v=*vEnd;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned id=noboost::get_id(G, *vIt);
            if(alpha_inv[id] == 0){
                if(label[id] > max){
                    max = (unsigned int) label[id];
                    v = *vIt;
                }
            }
        }
        unsigned id=noboost::get_id(G, v);
        visited[id] = true;
        alpha_inv[id] = i+1;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned idn=noboost::get_id(G, *nIt);
            if(alpha_inv[idn] == 0){
                reached_i[(int)label[idn]-1].push_back(idn);
                visited[idn] = true;
                label[idn] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                unsigned int w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[w], G); nIt != nEnd; nIt++){
                    unsigned idn=noboost::get_id(G, *nIt);
                    if(visited[idn]){
                        continue;
                    }

                    visited[idn] = true;
                    if((unsigned int)label[idn]-1 > j){
                        reached_i[(int)label[idn]].push_back(idn);
                        label[idn] += 0.5;
                        std::vector<unsigned int> edge;
                        unsigned idv=noboost::get_id(G, v);
                        edge.push_back(idv);
                        edge.push_back(idn);
                        fill_in_edges.push_back(edge);
                    }
                    else{
                        reached_i[j].push_back(idn);
                    }
                }
            }
        }

        for(unsigned int j = 0; j < label.size(); j++){
            label[j] = (float)roundf(label[j]);
            k = (k > (unsigned int)label[j])? k : (unsigned int)label[j];
        }
    }
}

template <typename G_t>
void LEX_M_minimal_ordering(G_t &G, std::vector<unsigned int> &alpha){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    unsigned int max = idxMap.size();
    alpha.resize(boost::num_vertices(G));
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<unsigned int> alpha_inv(max+1);
    std::vector<std::vector<unsigned int> > reached_i(max+1);

    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        label[id] = 1.0;
        alpha_inv[i++] = 0;
        visited[id] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned id=noboost::get_id(G, *vIt);
            if(alpha_inv[id] == 0){
                if((unsigned int)label[id] > max){
                    max = (unsigned int) label[id];
                    v = *vIt;
                }
            }
        }
        unsigned idv=noboost::get_id(G, v);
        visited[idv] = true;
        alpha[i] = idv;
        alpha_inv[idv] = i+1;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned idn=noboost::get_id(G, *nIt);
            if(alpha_inv[idn] == 0){
                reached_i[(int)label[idn]-1].push_back(idn);
                visited[idn] = true;
                label[idn] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                unsigned int w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[w], G); nIt != nEnd; nIt++){
                    unsigned idn=noboost::get_id(G, *nIt);
                    if(visited[idn]){
                        continue;
                    }

                    visited[idn] = true;
                    if((unsigned int)label[idn]-1 > j){
                        reached_i[(int)label[idn]].push_back(idn);
                        label[idn] += 0.5;
                    }
                    else{
                        reached_i[j].push_back(idn);
                    }
                }
            }
        }

        for(unsigned int j = 0; j < label.size(); j++){
            label[j] = (float)roundf(label[j]);
            k = (k > (unsigned int)label[j])? k : (unsigned int)label[j];
        }
    }
}

template <typename G_t>
void LEX_M_minimal_ordering(G_t &G, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &alpha){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    unsigned int max = idxMap.size();
    alpha.resize(boost::num_vertices(G));
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<unsigned int> alpha_inv(max+1);
    std::vector<std::vector<unsigned int> > reached_i(max+1);

    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        label[id] = 1.0;
        alpha_inv[i++] = 0;
        visited[id] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v=*vEnd;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned id=noboost::get_id(G, *vIt);
            if(alpha_inv[id] == 0){
                if((unsigned int)label[id] > max){
                    max = (unsigned int) label[id];
                    v = *vIt;
                }
            }
        }
        unsigned idv=noboost::get_id(G, v);
        visited[idv] = true;
        alpha[i] = v;
        alpha_inv[idv] = i+1;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned idn=noboost::get_id(G, *nIt);
            if(alpha_inv[idn] == 0){
                reached_i[(int)label[idn]-1].push_back(idn);
                visited[idn] = true;
                label[idn] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                unsigned int w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[w], G); nIt != nEnd; nIt++){
                    unsigned idn=noboost::get_id(G, *nIt);
                    if(visited[idn])
                        continue;

                    visited[idn] = true;
                    if((unsigned int)label[idn]-1 > j){
                        reached_i[(int)label[idn]].push_back(idn);
                        label[idn] += 0.5;
                    }
                    else{
                        reached_i[j].push_back(idn);
                    }
                }
            }
        }

        for(unsigned int j = 0; j < label.size(); j++){
            label[j] = (float)roundf(label[j]);
            k = (k > (unsigned int)label[j])? k : (unsigned int)label[j];
        }
    }
}

} //namespace treedec

#endif //ifdef TD_ELIMINATION_ORDERING

// vim:ts=8:sw=4:et
