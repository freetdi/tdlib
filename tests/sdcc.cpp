#include <set>
#include <tdlib/trace.hpp>

namespace treedec{
	struct bag_t;
}

typedef int var_t; // relevant here?
typedef std::set<unsigned int> sdcc_bagtype;
struct tree_dec_node
{
  std::set<unsigned int> bag;
  std::set<var_t> alive;
  unsigned weight;
};

#include <tdlib/treedec_traits.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_utility.hpp>

typedef bool cfg_alive_t;
typedef bool cfg_dying_t;

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

REGISTER_GRAPH_WITH_BUNDLED_BAGS(tree_dec_t, bag);

#include <tdlib/graph.hpp>
#include <tdlib/preprocessing.hpp>
#ifdef USE_GALA
#include <tdlib/exact_ta.hpp>
#endif
#include <boost/graph/copy.hpp>
#include <tdlib/thorup.hpp>

#ifdef USE_GALA
// some algorithms, need to move to comb::p17...
namespace choice{
  template<class G, template<class G_, class ...> class C>
  using exact_ta_=treedec::exact_ta<G, C>;
}
typedef treedec::draft::exact_decomposition<cfg_t,
             treedec::algo::default_config,
             choice::exact_ta_> alg_A;
#endif

typedef treedec::thorup<cfg_t> thorup;

template<class G, class A>
void do_it(G& g){
	tree_dec_t t;
	A B(g);
	B.do_it();
	B.get_tree_decomposition(t);
	assert(treedec::get_bagsize(t) == 3);
	auto nt=boost::num_vertices(t);

	boost::print_graph(t);

	for(unsigned i=0; i<nt; ++i){
		std::cout << i << ":";
		for(auto x: t[i].bag){
			std::cout << " " << x;
		}
		std::cout << "\n";
	}
}

int main()
{

	cfg_t g(4);
	treedec::add_edge(0, 1, g);
	treedec::add_edge(1, 2, g);
	boost::add_edge(2, 3, g);
	boost::add_edge(3, 0, g);
	assert(boost::num_edges(g)==4); // oops?
	assert(treedec::num_edges(g)==4);
	auto p=boost::edges(g);
	unsigned ii=0;
	for(; p.first!=p.second; ++p.first){
		++ii;
	}
	assert(ii==4);

	cfg_t h(g);
#if USE_GALA
	std::cout << "p17 test\n";
	//test<cfg_t, alg_A>(h);
	alg_A A(g);
	A.do_it(1);
	// A.get_tree_decomposition(t); almost
#endif

	std::cout << "thorup\n";
	do_it<cfg_t, thorup>(h);
}
