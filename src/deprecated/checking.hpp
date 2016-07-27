template <typename G_t, typename O_t>
bool is_minimum_degree_ordering(G_t &G, O_t &O){
    for(unsigned int i = 0; i < O.size(); i++){
        unsigned deg = boost::degree(O[i], G);
        unsigned min_deg = min_degree(G);
        if(deg != min_deg && min_deg > 0){
            std::cout << "deg: " << deg << ", min_deg: " << min_deg << std::endl;
            return false;
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(O[i], G);
        treedec::make_clique(nIt, nEnd, G);
        boost::clear_vertex(O[i], G);
    }
    return true;
}
