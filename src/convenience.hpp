#ifndef CONVENIENCE
#define CONVENIENCE

// TODO: inclusion is upside down.
#include "generic_base.hpp"
#include "generic_elimination_search_configs.hpp"
#include "preprocessing.hpp"

#ifdef HAVE_GALA_NOTYET
#include <gala/graph.h>
#include <gala/boost.h>
#endif

namespace treedec{

template <typename G_t>
void generic_elimination_search_CFG1(G_t const &G, unsigned max_nodes, unsigned max_orderings)
{
    gen_search::configs::CFG_DFS_1<G_t, algo::default_config>
       generic_elim_DFS_test (G);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
    size_t B=treedec::get_bagsize_of_elimination_ordering(H, generic_elim_DFS_test.ordering());
    assert(A == B);
}

template <typename G_t>
void generic_elimination_search_CFG2(G_t const &G, unsigned max_nodes, unsigned max_orderings)
{

    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G));
    ord_type cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true); // BUG

    gen_search::configs::CFG_DFS_2<G_t, algo::default_config>
      generic_elim_DFS_test (G);

//    gen_search::generic_elimination_search_DFS<G_t,
//        gen_search::configs::CFG_DFS_2 > //TODO: constructor...
//       generic_elim_DFS_test
//              (olay,
//               active,
//               ordering,
//               cur_ordering);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
    size_t B=treedec::get_bagsize_of_elimination_ordering(H, generic_elim_DFS_test.ordering());
    assert(A == B);
}


template <typename G_t>
void generic_elimination_search_CFG3(G_t const &G, unsigned max_nodes, unsigned max_orderings)
{
//    gen_search::overlay<G_t, G_t> olay(G);

    gen_search::configs::CFG_DFS_3<G_t, algo::default_config>
       generic_elim_DFS_test (G);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
    size_t B=treedec::get_bagsize_of_elimination_ordering(H, generic_elim_DFS_test.ordering());
    assert(A == B);
}


template <typename G_t>
void generic_elimination_search_p17(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    std::cout << "edges before PP: " << boost::num_edges(G) << std::endl;

    impl::preprocessing<G_t> PP(G);
    PP.do_it();
    PP.get_graph(G);

    std::cout << "edges after PP: " << boost::num_edges(G) << std::endl;

    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G));
    ord_type cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

/*
#ifdef HAVE_GALA_NOTYET
    typedef gala::graph<std::vector, std::vector, uint32_t> ssg_vec_vec32i;

    typedef G_t Underlying_t;
    typedef ssg_vec_vec32i  Overlay_t;

    overlay<Underlying_t, Overlay_t> olay(G, active);
#else

//    gen_search::overlay<Underlying_t, Overlay_t> olay(G);
#endif
*/

    gen_search::configs::CFG_DFS_p17<G_t, algo::default_config>
       generic_elim_DFS_test (G /* ... more? */);

    //gen_search::generic_elimination_search_DFS<G_t, gen_search::configs::CFG_DFS_2 >
    //   generic_elim_DFS_test
    //          (olay, // FIXME.
    //           active,
    //           ordering,
    //           cur_ordering);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
    size_t B=treedec::get_bagsize_of_elimination_ordering(H, generic_elim_DFS_test.ordering());
    assert(A == B);
}

template <typename G_t>
void generic_elimination_search_p17_jumper(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    std::cout << "vertices before PP: " << boost::num_vertices(G) << std::endl;
    std::cout << "edges before PP: " << boost::num_edges(G) << std::endl;

    impl::preprocessing<G_t> PP(G);
    PP.do_it();

    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags;

    PP.get_bags(bags);

    G_t G2;
    PP.get_graph(G2);

    std::cout << "vertices after PP: " << boost::num_vertices(G2) << std::endl;
    std::cout << "edges after PP: " << boost::num_edges(G2) << std::endl;


    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G2));
    ord_type cur_ordering(boost::num_vertices(G2));

    std::vector<BOOL> active(boost::num_vertices(G2), true);

/*
#ifdef HAVE_GALA_NOTYET
    typedef gala::graph<std::vector, std::vector, uint32_t> ssg_vec_vec32i;

    typedef G_t Underlying_t;
    typedef ssg_vec_vec32i  Overlay_t;

    overlay<Underlying_t, Overlay_t> olay(G2, active);
#else

//    gen_search::overlay<Underlying_t, Overlay_t> olay(G);
#endif
*/

    gen_search::configs::CFG_DFS_p17<G_t, algo::default_config>
       generic_elim_DFS_test (G2 /* ... more? */);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    unsigned lb=generic_elim_DFS_test.global_lower_bound_bagsize();
    unsigned ub=generic_elim_DFS_test.global_upper_bound_bagsize();
    ord_type best=generic_elim_DFS_test.ordering();

    while(true){
        gen_search::configs::CFG_DFS_p17_2<G_t, algo::default_config>
           generic_elim_DFS_test2 (G2);

        generic_elim_DFS_test2.set_max_nodes_generated(max_nodes);
        generic_elim_DFS_test2.set_max_orderings_generated(max_orderings);

        generic_elim_DFS_test2.set_lb(lb);
        generic_elim_DFS_test2.set_ub(ub);
        generic_elim_DFS_test2.set_best_ordering(best);

        generic_elim_DFS_test2.do_it();

        lb=generic_elim_DFS_test2.global_lower_bound_bagsize();
        if(generic_elim_DFS_test2.global_upper_bound_bagsize() == ub){
            return;
        }
        ub=generic_elim_DFS_test2.global_upper_bound_bagsize();
        best=generic_elim_DFS_test2.ordering();
    }

    G_t H(G2);
    size_t A=ub;
    size_t B=treedec::get_bagsize_of_elimination_ordering(H, best);
    assert(A == B);
}

} //namespace treedec



#endif //guard

// vim:ts=8:sw=4:et
