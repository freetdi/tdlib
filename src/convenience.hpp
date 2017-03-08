//TODO:guard

#include <gala/graph.h>
#include <gala/boost.h>

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
    std::vector<bool> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, G_t, configs::CFG_DFS_1<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
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

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

/*
    assert(treedec::get_width_of_elimination_ordering(G, ordering)
            == generic_elim_DFS_test.global_upper_bound_bagsize());
*/
    std::cout << "width of elimination ordering (check): " << treedec::get_width_of_elimination_ordering(G, ordering) << std::endl;

}




template <typename G_t>
void generic_elimination_search_CFG2(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<bool> active(boost::num_vertices(G), true);

    typedef gala::graph<std::vector, std::vector, uint32_t> ssg_vec_vec32i;

    typedef G_t             Underlying_t;
    typedef ssg_vec_vec32i  Overlay_t;


//    overlay<G_t, G_t> olay(G, active);

    overlay<Underlying_t, Overlay_t> olay(G, active);

    generic_elimination_search_DFS<Underlying_t, Overlay_t, configs::CFG_DFS_2<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
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

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    std::cout << "width of elimination ordering (check): " << treedec::get_width_of_elimination_ordering(G, ordering) << std::endl;
/*
    assert(treedec::get_width_of_elimination_ordering(G, ordering)
            == generic_elim_DFS_test.global_upper_bound_bagsize());
*/
}


template <typename G_t>
void generic_elimination_search_CFG3(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<bool> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, G_t, configs::CFG_DFS_3<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               ordering,
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

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    std::cout << "width of elimination ordering (check): " << treedec::get_width_of_elimination_ordering(G, ordering) << std::endl;

    assert(treedec::get_width_of_elimination_ordering(G, ordering)
            == generic_elim_DFS_test.global_upper_bound_bagsize());
}

} //namespace gen_search

} //namespace treedec



