template <typename G_t>
bool is_valid_coloring(G_t &G, std::vector<int> &coloring){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned int spos = noboost::get_pos(boost::source(*eIt, G), G);
        unsigned int tpos = noboost::get_pos(boost::target(*eIt, G), G);
        if(coloring[spos] == coloring[tpos]){
            return false;
        }
    }
    return true;
}
