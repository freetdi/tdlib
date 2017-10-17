
#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <tdlib/preprocessing.hpp>
#include <tdlib/graph.hpp>
#include <boost/graph/copy.hpp>
#ifdef USE_GALA
#include <tdlib/exact_ta.hpp>
#endif


typedef int var_t; // relevant here?
typedef bool cfg_alive_t;
typedef bool cfg_dying_t;

// import sdcc types
struct tree_dec_node
{
  std::set<unsigned int> bag;
  std::set<var_t> alive;
  // assignment_list_t assignments; // not relevant
  unsigned weight;
};
struct cfg_node
{
  // iCode *ic; // bnot relevant
  // operand_map_t operands;
  cfg_alive_t alive;
  cfg_dying_t dying;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, cfg_node> cfg_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, tree_dec_node> tree_dec_t;

// some algorithms
#ifdef USE_GALA
namespace choice{
  template<class G, template<class G_, class ...> class C>
  using exact_ta_=treedec::exact_ta<G, C>;
}
typedef treedec::draft::exact_decomposition<cfg_t,
             treedec::algo::default_config,
             choice::exact_ta_> alg_A;
#endif

int main()
{
	cfg_t g(3);
	treedec::add_edge(0, 1, g);
	treedec::add_edge(1, 2, g);
	treedec::add_edge(2, 0, g);
	assert(boost::num_edges(g)==6);
	assert(treedec::num_edges(g)==3);

#ifdef USE_GALA
	alg_A A(g);
	A.do_it(1);
#endif

//	std::cout << A.get_treewidth();

}
