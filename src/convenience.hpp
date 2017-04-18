#ifndef CONVENIENCE
#define CONVENIENCE

// BUG: inclusion is upside down.
#include "generic_base.hpp"
#include "generic_elimination_search.hpp"

namespace treedec{

namespace gen_search{

template <typename G_t>
void generic_elimination_search_CFG1(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_1> //TODO: constructor...
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
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl; //bug: is not propagated

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize()
			 ==treedec::get_bagsize_of_elimination_ordering(H, ordering));
}




template <typename G_t>
void generic_elimination_search_CFG2(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_2 > //TODO: constructor...
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

    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
	 size_t B=treedec::get_bagsize_of_elimination_ordering(H, ordering);
    if(A != B){
		 unreachable();
		 std::cerr << A << " vs " << B << "\n";
	 }else{
	//	 std::cerr << "ok: " << A << " vs " << B << "\n";
	 }
}


template <typename G_t>
void generic_elimination_search_CFG3(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

    overlay<G_t, G_t> olay(G, active);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_3 > //TODO: constructor...
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

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize()
           ==treedec::get_bagsize_of_elimination_ordering(H, ordering)+1);
}


template <typename G_t>
void generic_elimination_search_CFG4(G_t &G, unsigned max_nodes, unsigned max_orderings){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

#ifdef HAVE_GALA_NOTYET
    typedef gala::graph<std::vector, std::vector, uint32_t> ssg_vec_vec32i;

    typedef G_t Underlying_t;
    typedef ssg_vec_vec32i  Overlay_t;

    overlay<Underlying_t, Overlay_t> olay(G, active);
#else
    typedef G_t Underlying_t;
    typedef G_t Overlay_t;

    overlay<Underlying_t, Overlay_t> olay(G, active);
#endif

    generic_elimination_search_DFS<Underlying_t, configs::CFG_DFS_2 > //TODO: constructor...
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

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize() == treedec::get_bagsize_of_elimination_ordering(H, ordering));
}


} //namespace gen_search

} //namespace treedec



#endif //guard
