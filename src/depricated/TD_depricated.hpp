
template <typename G_t>
void all_k_colorings(unsigned int n, unsigned int k, unsigned int max,
        std::set<unsigned int> &M, std::vector<std::vector<int> > &colorings, G_t &G)
{
    std::vector<int> coloring(n, -1);
    std::set<unsigned int>::iterator iM = M.begin();
    while(iM != M.end()){
        coloring[*(iM++)]++;
    }

    iM = M.begin();
    unsigned int c = 0;

    if(treedec::app::is_valid_coloring(G, M, coloring)){
        colorings[c++] = coloring;
    }

    while(iM != M.end() && c < colorings.size()){
        if(coloring[*iM] < k-1){
            if(iM == M.end()){ break; }
            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }
        }
        else{
            while(coloring[*iM] == k-1 && iM != M.end()){
                coloring[*iM] = 0;
                iM++;
            }
            if(iM == M.end()){ break; }

            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }

            iM = M.begin();
        }
    }

    colorings.resize(c);
}


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

template <typename G_t>
bool is_valid_vertex_cover(G_t &G, typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() || C.find(boost::target(*eIt, G)) != C.end()){
        }
        else{ return false;
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_independent_set(G_t &G, typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() && C.find(boost::target(*eIt, G)) != C.end()){
            return false;
        }
    }
    return true;
}
