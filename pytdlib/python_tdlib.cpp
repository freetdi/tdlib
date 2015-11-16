#include <boost/tuple/tuple.hpp>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include "TD_preprocessing.hpp"
#include "TD_combinations.hpp"
#include "TD_lower_bounds.hpp"
#include "TD_elimination_orderings.hpp"
#include "TD_misc.hpp"


#ifndef TD_STRUCT_VERTEX
#define TD_STRUCT_VERTEX

struct Vertex{
    unsigned int id;
};

#endif

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;

struct bag{
    std::set<unsigned int> bag;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag> TD_tree_dec_t;

#include "python_tdlib.hpp"


void make_tdlib_graph(TD_graph_t &G, std::vector<unsigned int> &V, std::vector<unsigned int> &E){
    unsigned int max = 0;
    for(unsigned int i = 0; i < V.size(); i++)
        max = (V[i]>max)? V[i] : max;

    std::vector<TD_graph_t::vertex_descriptor> idxMap(max+1);

    for(unsigned int i = 0; i < V.size(); i++){
        idxMap[V[i]] = boost::add_vertex(G);
        G[idxMap[V[i]]].id = V[i];
    }

    if(E.size() != 0){
        for(unsigned int j = 0; j < E.size()-1; j++){
            boost::add_edge(idxMap[E[j]], idxMap[E[j+1]], G);
            j++;
        }
    }
}

void make_tdlib_decomp(TD_tree_dec_t &T, std::vector<std::vector<int> > &V, std::vector<unsigned int> &E){
    std::vector<TD_tree_dec_t::vertex_descriptor> idxMap(V.size()+1);

    for(unsigned int i = 0; i < V.size(); i++){
        idxMap[i] = boost::add_vertex(T);
        std::set<unsigned int> bag;
        for(unsigned int j = 0; j < V[i].size(); j++)
            bag.insert((unsigned int) V[i][j]);
        T[idxMap[i]].bag = bag;
    }

    if(E.size() != 0){
        for(unsigned int j = 0; j < E.size()-1; j++){
            boost::add_edge(idxMap[E[j]], idxMap[E[j+1]], T);
            j++;
        }
    }

}

void make_python_graph(TD_graph_t &G, std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    boost::graph_traits<TD_graph_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        V_G.push_back(G[*vIt].id);

    boost::graph_traits<TD_graph_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        E_G.push_back(G[boost::source(*eIt, G)].id);
        E_G.push_back(G[boost::target(*eIt, G)].id);
    }
}

void make_python_decomp(TD_tree_dec_t &T, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int> vertex_map;
    boost::graph_traits<TD_tree_dec_t>::vertex_iterator tIt, tEnd;
    unsigned int id = 0;
    
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        vertex_map.insert(std::pair<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>(*tIt, id++));
        std::vector<int> bag;
        for(std::set<unsigned int>::iterator sIt = T[*tIt].bag.begin(); sIt != T[*tIt].bag.end(); sIt++)
            bag.push_back((int)*sIt);
        V_T.push_back(bag);
    }
    
    boost::graph_traits<TD_tree_dec_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(T); eIt != eEnd; eIt++){
        std::map<boost::graph_traits<TD_tree_dec_t>::vertex_descriptor, unsigned int>::iterator v, w;
        v = vertex_map.find(boost::source(*eIt, T));
        w = vertex_map.find(boost::target(*eIt, T));
        E_T.push_back(v->second);
        E_T.push_back(w->second);
    }
}


/* PREPROCESSING */

int gc_preprocessing(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &bags, int lb){
    TD_graph_t G, H;
    make_tdlib_graph(G, V_G, E_G);

    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > td_bags;
    treedec::preprocessing(G, td_bags, lb);

    V_G.clear();
    E_G.clear();

    treedec::remove_isolated_vertices(H, G);
    G = H;

    make_python_graph(G, V_G, E_G);

    for(unsigned int i = 0; i < td_bags.size(); i++){
        std::vector<int> bag;
        bag.push_back(td_bags[i].get<0>());
        for(std::set<unsigned int>::iterator sIt = td_bags[i].get<1>().begin(); sIt != td_bags[i].get<1>().end(); sIt++)
            bag.push_back((int)*sIt);
        bags.push_back(bag);
    }   

    return lb;
}


int gc_PP_MD(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::PP_MD(G, T, lb);

    treedec::make_small(T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}


int gc_PP_FI_TM(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::PP_FI_TM(G, T, lb);

    treedec::make_small(T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

int gc_preprocessing_glue_bags(std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, std::vector<std::vector<int> > &bags){
    TD_tree_dec_t T;

    make_tdlib_decomp(T, V_T, E_T);

    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > td_bags;
    for(unsigned int i = bags.size(); i > 0; i--){
        unsigned int v = (unsigned int)bags[i-1][0];
        std::set<unsigned int> bag;
        for(unsigned int j = 1; j < bags[i-1].size(); j++)
            bag.insert(bags[i-1][j]);

        td_bags.push_back(boost::tuple<unsigned int, std::set<unsigned int> >(v, bag));
    }

    treedec::preprocessing_glue_bags(td_bags, T);

    V_T.clear();
    E_T.clear();

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

/* LOWER BOUNDS */


int gc_deltaC_min_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::deltaC_min_d(G);
}

int gc_deltaC_max_d(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::deltaC_max_d(G);
}

int gc_deltaC_least_c(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::deltaC_least_c(G);
}

int gc_LBN_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBN_deltaC(G);
}

int gc_LBNC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBNC_deltaC(G);
}

int gc_LBP_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBP_deltaC(G);
}

int gc_LBPC_deltaC(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    return treedec::lb::LBPC_deltaC(G);
}


/* EXACT TREE DECOMPOSITIONS */

int gc_exact_decomposition_cutset(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::exact_decomposition_cutset(G, T, lb);

    treedec::make_small(T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}


int gc_exact_decomposition_dynamic(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, int lb){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::exact_decomposition_dynamic(G, T, lb);

    treedec::make_small(T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

/* APPOXIMATIVE TREE DECOMPOSITIONS */


int gc_seperator_algorithm(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::seperator_algorithm(G, T);

    treedec::make_small(T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

void gc_minDegree_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E, std::vector<unsigned int> &elim_ordering){
    TD_graph_t G;
    make_tdlib_graph(G, V, E);

    treedec::minDegree_ordering(G, elim_ordering);
}

void gc_fillIn_ordering(std::vector<unsigned int> &V, std::vector<unsigned int> &E, std::vector<unsigned int> &elim_ordering){
    TD_graph_t G;
    make_tdlib_graph(G, V, E);

    treedec::fillIn_ordering(G, elim_ordering);
}

int gc_ordering_to_treedec(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T, std::vector<unsigned int> &elim_ordering){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;
    treedec::ordering_to_treedec(G, elim_ordering, T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

void gc_treedec_to_ordering(std::vector<std::vector<int> > &V, std::vector<unsigned int> &E, std::vector<unsigned int> &elim_ordering){
    TD_tree_dec_t T;
    make_tdlib_decomp(T, V, E);

    return;

    treedec::treedec_to_ordering(T, elim_ordering);
}


/* POSTPROCESSING */

int gc_MSVS(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;	
    make_tdlib_decomp(T, V_T, E_T);

    treedec::MSVS(G, T);

    V_T.clear();
    E_T.clear();

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

void gc_minimalChordal(std::vector<unsigned int> &V, std::vector<unsigned int> &E, std::vector<unsigned int> &old_elimination_ordering, std::vector<unsigned int> &new_elimination_ordering){
    TD_graph_t G;
    make_tdlib_graph(G, V, E);

    typename std::vector<TD_graph_t::vertex_descriptor> new_elimination_ordering_tmp;

    treedec::minimalChordal(G, old_elimination_ordering, new_elimination_ordering_tmp);

    TD_graph_t::vertex_iterator vIt, vEnd;
    for(unsigned int i = 0; i < new_elimination_ordering_tmp.size(); i++)
        new_elimination_ordering.push_back(G[new_elimination_ordering_tmp[i]].id);
}


/* MISC */


int gc_is_valid_treedecomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;
    make_tdlib_decomp(T, V_T, E_T);

    return treedec::is_valid_treedecomposition(G, T);
}

int gc_trivial_decomposition(std::vector<unsigned int> &V_G, std::vector<unsigned int> &E_G, std::vector<std::vector<int> > &V_T, std::vector<unsigned int> &E_T){
    TD_graph_t G;
    make_tdlib_graph(G, V_G, E_G);

    TD_tree_dec_t T;

    treedec::trivial_decomposition(G, T);

    make_python_decomp(T, V_T, E_T);

    return treedec::get_width(T);
}

int gc_get_width(std::vector<std::vector<int> > &V_T){
    int width = 0;
    for(unsigned int i = 0; i< V_T.size(); i++)
        width = ((int)V_T[i].size() > width)? (int)V_T[i].size() : width;

    return width-1;
}
