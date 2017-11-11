#include "config.h"
#include <set>
#include <treedec/trace.hpp>

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

#include <treedec/treedec_traits.hpp>
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
                              GDIR, cfg_node> cfg_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              TDIR, tree_dec_node> tree_dec_t;

TREEDEC_TREEDEC_BAG_TRAITS(tree_dec_t, bag);

#include <treedec/graph.hpp>
#include <treedec/preprocessing.hpp>
#include <boost/graph/copy.hpp>


#include <treedec/thorup.hpp>
#include <treedec/combinations.hpp>


typedef treedec::he::fill_in<cfg_t> FI;
typedef treedec::he::thorup<cfg_t> thorup;
typedef treedec::comb::PP_MD<cfg_t> PP_MD;
typedef treedec::pending::PP_FI<cfg_t> PP_FI;
typedef treedec::pending::PP_FI_TM<cfg_t> PP_FI_TM;

#ifdef USE_GALA
typedef treedec::comb::ex17<cfg_t> ppta;
#endif


#include <treedec/nice_decomposition.hpp>

template<class G, class A>
static void do_it(G const& g_){
	auto g=g_;
	tree_dec_t t;
	A B(g);
	B.do_it();
	B.get_tree_decomposition(t);
	boost::print_graph(t);
	//assert(treedec::get_bagsize(t) == 3);
	std::cout << treedec::get_bagsize(t) << "\n";
	auto nt=boost::num_vertices(t);

	assert(treedec::is_valid_treedecomposition(g, t));

	for(unsigned i=0; i<nt; ++i){
		std::cout << i << ":";
		for(auto x: t[i].bag){
			std::cout << " " << x;
		}
		std::cout << "\n";
	}

#ifndef NO_NICIFY
	treedec::nice::nicify(t);
	assert(treedec::is_valid_treedecomposition(g, t));
#else
	incomplete();
#endif
}

#ifndef INLINE_CC
int main()
{
	const unsigned n=16;
	cfg_t h(n);
#include "g.h"
	cfg_t g;
	boost::copy_graph(h, g);
//	boost::copy_graph(h, g);
	boost::add_vertex(g);

	std::cout << "==test graph==\n";
	boost::print_graph(g);
	std::cout << "====\n";

	std::cout << "PP+FI+TM\n";
	do_it<cfg_t, PP_FI_TM>(g);

	std::cout << "FI\n";
	do_it<cfg_t, FI>(g);

	std::cout << "thorup\n";
	do_it<cfg_t, thorup>(g);

	std::cout << "FI\n";
	do_it<cfg_t, FI>(g);

	std::cout << "PP+MD\n";
	do_it<cfg_t, PP_MD>(g);
//
#ifdef USE_GALA
	std::cout << "ppta\n";
	do_it<cfg_t, ppta>(g);
#endif

	std::cout << "PP+FI\n";
	do_it<cfg_t, PP_FI>(g);




}
#endif
