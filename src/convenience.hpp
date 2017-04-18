#ifndef CONVENIENCE
#define CONVENIENCE

// BUG: inclusion is upside down.
#include "generic_base.hpp"
#include "generic_elimination_search.hpp"

// TODO: what's this? move to test?

namespace treedec{

// namespace draft?

template <typename G_t>
void generic_elimination_search_CFG1(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G));
    ord_type cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true); // BUG.

    gen_search::overlay<G_t, G_t> olay(G);

    //TODO: constructor...
    gen_search::generic_elimination_search_DFS<G_t, gen_search::configs::CFG_DFS_1>
       generic_elim_DFS_test
              (olay,
               active, // BUG: optional.
               ordering,
               cur_ordering);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize()
			 ==treedec::get_bagsize_of_elimination_ordering(H, ordering));
}

template <typename G_t>
void generic_elimination_search_CFG2(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G));
    ord_type cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true); // BUG

    gen_search::overlay<G_t, G_t> olay(G);

    gen_search::generic_elimination_search_DFS<G_t,
        gen_search::configs::CFG_DFS_2 > //TODO: constructor...
       generic_elim_DFS_test
              (olay,
               active,
               ordering,
               cur_ordering);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();


    G_t H(G);
    size_t A=generic_elim_DFS_test.global_upper_bound_bagsize();
	 size_t B=treedec::get_bagsize_of_elimination_ordering(H, ordering);
	 if(A != B){
                unreachable();
                std::cerr << A << " vs " << B << "\n";
        }else{
       //       std::cerr << "ok: " << A << " vs " << B << "\n";
        }

}


template <typename G_t>
void generic_elimination_search_CFG3(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    gen_search::overlay<G_t, G_t> olay(G);

    //TODO: constructor...
    gen_search::generic_elimination_search_DFS<G_t, gen_search::configs::CFG_DFS_3 >
       generic_elim_DFS_test (olay);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize()
           ==treedec::get_bagsize_of_elimination_ordering(H, generic_elim_DFS_test.ordering())+1);
}


template <typename G_t>
void generic_elimination_search_CFG4(G_t &G, unsigned max_nodes, unsigned max_orderings)
{
    typedef std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ord_type;
    ord_type ordering(boost::num_vertices(G));
    ord_type cur_ordering(boost::num_vertices(G));

    std::vector<BOOL> active(boost::num_vertices(G), true);

#ifdef HAVE_GALA_NOTYET
    typedef gala::graph<std::vector, std::vector, uint32_t> ssg_vec_vec32i;

    typedef G_t Underlying_t;
    typedef ssg_vec_vec32i  Overlay_t;

    overlay<Underlying_t, Overlay_t> olay(G, active);
#else
    typedef G_t Underlying_t;
    typedef G_t Overlay_t;

    gen_search::overlay<Underlying_t, Overlay_t> olay(G);
#endif

    //TODO: constructor...
    gen_search::generic_elimination_search_DFS<G_t, gen_search::configs::CFG_DFS_2 >
       generic_elim_DFS_test
              (olay, // FIXME.
               active,
               ordering,
               cur_ordering);

    generic_elim_DFS_test.set_max_nodes_generated(max_nodes);
    generic_elim_DFS_test.set_max_orderings_generated(max_orderings);

    generic_elim_DFS_test.do_it();

    G_t H(G);
    assert(generic_elim_DFS_test.global_upper_bound_bagsize() == treedec::get_bagsize_of_elimination_ordering(H, ordering));
}

} //namespace treedec



#endif //guard
// vim:ts=8:sw=4:et
