#include <tdlib/graph.hpp>

template <typename G_t>
bool is_valid_clique(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    for(typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt1 = C.begin(); sIt1 != C.end(); sIt1++){
        typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != C.end(); sIt2++){
            if(!boost::edge(*sIt1, *sIt2, G).second){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_vertex_cover(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() || C.find(boost::target(*eIt, G)) != C.end()){
        }
        else{
            return false;
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_independent_set(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() && C.find(boost::target(*eIt, G)) != C.end()){
            return false;
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_dominating_set(G_t &G, typename std::set<unsigned int> &set)
{
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(set.find(*vIt) == set.end()){
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            bool hit = false;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                if(set.find(*nIt) != set.end()){
                    hit = true;
                    break;
                }
            }
            if(!hit){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_coloring(G_t &G, std::vector<typename treedec::treedec_traits<TD_tree_dec_t>::bag_type> &vec)
{
    std::vector<bool> visited(boost::num_vertices(G), false);
    for(unsigned i = 0; i < vec.size(); i++){
        for(typename treedec::treedec_traits<TD_tree_dec_t>::bag_type::iterator sIt = vec[i].begin(); sIt != vec[i].end(); sIt++){
            unsigned pos = treedec::get_pos(*sIt, G);
            visited[pos] = true;
        }
    }
    for(unsigned i = 0; i < visited.size(); i++){
        if(!visited[i]){
            return false;
        }
    }

    std::vector<unsigned> col(boost::num_vertices(G));
    for(unsigned i = 0; i < vec.size(); i++){
        for(typename treedec::treedec_traits<TD_tree_dec_t>::bag_type::iterator sIt = vec[i].begin(); sIt != vec[i].end(); sIt++){
            col[treedec::get_pos(*sIt, G)] = i;
        }
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(col[treedec::get_pos(boost::source(*eIt, G), G)] == col[treedec::get_pos(boost::target(*eIt, G), G)]){
            return false;
        }
    }
    return true;
}
