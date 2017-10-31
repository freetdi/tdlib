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
#include <boost/graph/copy.hpp>


#include <tdlib/thorup.hpp>
#include <tdlib/combinations.hpp>


typedef treedec::comb::PP_MD<cfg_t> PP_MD;
typedef treedec::comb::PP_FI<cfg_t> PP_FI;
typedef treedec::comb::PP_FI_TM<cfg_t> PP_FI_TM;
typedef treedec::thorup<cfg_t> thorup;

#ifdef USE_GALA
typedef treedec::comb::ex17<cfg_t> ppta;
#endif


#include <tdlib/nice_decomposition.hpp>

template<class G, class A>
static void do_it(G& g){
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

	treedec::nice::nicify(t);
}

#ifndef INLINE_CC
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

	cfg_t h1(g);
	cfg_t h2(g);
	cfg_t h3(g);
	cfg_t h4(g);
	cfg_t h5(g);

	std::cout << "thorup\n";
	do_it<cfg_t, thorup>(h1);

	std::cout << "PP+MD\n";
	do_it<cfg_t, PP_MD>(h2);

	std::cout << "PP+FI\n";
	do_it<cfg_t, PP_FI>(h3);

	std::cout << "PP+FI\n";
	do_it<cfg_t, PP_FI_TM>(h4);


#if USE_GALA
	std::cout << "ppta\n";
	//test<cfg_t, alg_A>(h);
	ppta A(h5);
	A.do_it(1);
	// A.get_tree_decomposition(t); almost
#endif


}
#endif
