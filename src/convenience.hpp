namespace treedec{

/* TODO: cleanup?!
template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T,
                   unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
{
    return minDegree_decomp(
        G, T,
        (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL,
        ub, ignore_isolated_vertices
       );
}

template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   minDegree_decomp(G_t &G)
{
    return minDegree_decomp(G, (typename treedec_chooser<G_t>::type*)NULL);
}


template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  fillIn_decomp(G_t &G, T_t &T,
                unsigned ub = UINT_MAX, bool ignore_isolated_vertices=false)
{
    return fillIn_decomp
      (G, T,
       (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL,
       ub, ignore_isolated_vertices);
}

template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   fillIn_decomp(G_t &G)
{
    return fillIn_decomp(G, (typename treedec_chooser<G_t>::type*)NULL);
}
*/


namespace gen_search{

template <typename G_t>
void generic_elimination_search_CFG1(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<bool> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, G_t, configs::CFG_DFS_1<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
               cur_ordering,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,                       //nodes generated
               0                        //orderings generated
       );

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

/*
    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    G_t H(G);

    int w1_check = treedec::get_width_of_elimination_ordering(H, ordering)+1;

    std::cout << "width of elimination ordering (check): " << w1_check << std::endl;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering2(boost::num_vertices(G));

    std::cout << "postrefiner..";

    int w2 = treedec::minimalChordal(G, ordering, ordering2)+1;

    std::cout << " done." << std::endl;

    int w2_check = treedec::get_width_of_elimination_ordering(G, ordering2)+1;

    std::cout << "width of elimination refined ordering: " << w2_check << std::endl;

    if(w1_check != generic_elim_DFS_test.global_upper_bound_bagsize() || w2 != w2_check){
        std::cout << "widthcheck error!!!!!" << std::endl;
        std::cout << "w1: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
        std::cout << "w1_check: " << w1_check << std::endl;
        std::cout << "w2: " << w2 << std::endl;
        std::cout << "w2_check: " << w2_check << std::endl;
        exit(-666);
    }
*/
}




template <typename G_t>
void generic_elimination_search_CFG2(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<bool> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, G_t, configs::CFG_DFS_2<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
               cur_ordering,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,                       //nodes generated
               0                        //orderings generated
       );

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

/*
    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    G_t H(G);

    int w1_check = treedec::get_width_of_elimination_ordering(H, ordering)+1;

    std::cout << "width of elimination ordering (check): " << w1_check << std::endl;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering2(boost::num_vertices(G));

    std::cout << "postrefiner..";

    int w2 = treedec::minimalChordal(G, ordering, ordering2)+1;

    std::cout << " done." << std::endl;

    int w2_check = treedec::get_width_of_elimination_ordering(G, ordering2)+1;

    std::cout << "width of elimination refined ordering: " << w2_check << std::endl;

    if(w1_check != generic_elim_DFS_test.global_upper_bound_bagsize() || w2 != w2_check){
        std::cout << "widthcheck error!!!!!" << std::endl;
        std::cout << "w1: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
        std::cout << "w1_check: " << w1_check << std::endl;
        std::cout << "w2: " << w2 << std::endl;
        std::cout << "w2_check: " << w2_check << std::endl;

        exit(-666);
    }
*/
}


template <typename G_t>
void generic_elimination_search_CFG3(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<bool> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, G_t, configs::CFG_DFS_3<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
               cur_ordering,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,                       //nodes generated
               0                        //orderings generated
       );

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

/*
    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    G_t H(G);

    int w1_check = treedec::get_width_of_elimination_ordering(H, ordering)+1;

    std::cout << "width of elimination ordering (check): " << w1_check << std::endl;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering2(boost::num_vertices(G));

    std::cout << "postrefiner..";

    int w2 = treedec::minimalChordal(G, ordering, ordering2)+1;

    std::cout << " done." << std::endl;

    int w2_check = treedec::get_width_of_elimination_ordering(G, ordering2)+1;

    std::cout << "width of elimination refined ordering: " << w2_check << std::endl;

    if(w1_check != generic_elim_DFS_test.global_upper_bound_bagsize() || w2 != w2_check){
        std::cout << "widthcheck error!!!!!" << std::endl;
        std::cout << "w1: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
        std::cout << "w1_check: " << w1_check << std::endl;
        std::cout << "w2: " << w2 << std::endl;
        std::cout << "w2_check: " << w2_check << std::endl;
        exit(-666);
    }
*/
}

} //namespace gen_search

} //namespace treedec



