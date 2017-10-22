#include <set>
#include <tdlib/trace.hpp>
#include <boost/graph/properties.hpp>

namespace treedec{
	struct bag_t;
}

typedef int var_t; // relevant here?
typedef std::set<unsigned int> sdcc_bagtype;
struct tree_dec_node
{
//  typedef boost::vertex_property_tag kind; // do we need it?
  std::set<unsigned int> bag;
  std::set<var_t> alive;
  // assignment_list_t assignments; // not relevant
  unsigned weight;

  //tree_dec_node& operator=(const boost::property<treedec::bag_t, std::set<unsigned int> >&);
};


#include <tdlib/graph_traits.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

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

typedef treedec::thorup<cfg_t> alg_B;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			  boost::property<treedec::bag_t, std::vector<unsigned> > > prop_tdt;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, tree_dec_node> sbib_tdt;

void dryrun()
{
	typedef boost::adjacency_list<boost::vecS, boost::vecS,
			  boost::undirectedS, treedec::bag_t> bag_t_tdt;

	prop_tdt p;
	bag_t_tdt b;
	tree_dec_t t;
	boost::copy_graph(p, t);
	boost::copy_graph(b, t);

}

#if 0
namespace boost{

  inline bagstuff::gtob<tree_dec_t>::type \
  get(treedec::bag_t, tree_dec_t const&t, unsigned k)\
  { untested();\
	  return t[k].bag;\
  }\

}
#endif

int main()
{
	dryrun();

	cfg_t g(4);
	treedec::add_edge(0, 1, g);
	treedec::add_edge(1, 2, g);
	treedec::add_edge(2, 3, g);
	treedec::add_edge(3, 0, g);
	assert(boost::num_edges(g)==8);
	assert(treedec::num_edges(g)==4);
	auto p=boost::edges(g);
	unsigned ii=0;
	for(; p.first!=p.second; ++p.first){
		++ii;
	}
	assert(ii==8);


	cfg_t h(g);
	tree_dec_t t;
#if USE_GALA
	std::cout << "p17 test\n";
	alg_A A(g);
	A.do_it(1);
	// A.get_tree_decomposition(t); almost
#endif

	std::cout << "thorup test\n";
	alg_B B(h);
	B.do_it();
	B.get_tree_decomposition(t);
	std::cout << "thorup test done\n";

	boost::get(treedec::bag_t(), t); // test
	auto const&ct(t);
	boost::get(treedec::bag_t(), ct); // test

	assert(treedec::get_bagsize(t) == 3);

//	std::cout << A.get_treewidth();

}
