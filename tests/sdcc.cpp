#include <set>
#include <boost/graph/properties.hpp>
#include <tdlib/trace.hpp>

namespace treedec{
	struct bag_t;
}

typedef int var_t; // relevant here?
typedef std::set<unsigned int> sdcc_bagtype;
struct tree_dec_node
{
  typedef boost::vertex_property_tag kind; // do we need it?
  std::set<unsigned int> bag;
  std::set<var_t> alive;
  // assignment_list_t assignments; // not relevant
  unsigned weight;

  // HACK HACK. looking for a sulution...
  tree_dec_node& operator=(std::set<unsigned int> const& b){ untested();
	  bag = b;
	  return *this;
  }
  tree_dec_node& operator=(const boost::property<treedec::bag_t, std::set<unsigned int> >&);
};


#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <tdlib/preprocessing.hpp>
#include <tdlib/graph.hpp>
#include <boost/graph/copy.hpp>
#include <tdlib/thorup.hpp>
#ifdef USE_GALA
#include <tdlib/exact_ta.hpp>
#endif

tree_dec_node& tree_dec_node::operator=(
		const boost::property<treedec::bag_t, std::set<unsigned int> >& p){
	bag=p.m_value;
	return *this;
}


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

namespace boost{
sdcc_bagtype& get(treedec::bag_t, tree_dec_t& t, unsigned u)
{ untested();
	return t[u].bag;
}
}

#ifdef USE_GALA
// some algorithms, need to move to some header, ex::p17...
namespace choice{
  template<class G, template<class G_, class ...> class C>
  using exact_ta_=treedec::exact_ta<G, C>;
}
typedef treedec::draft::exact_decomposition<cfg_t,
             treedec::algo::default_config,
             choice::exact_ta_> alg_A;
#endif

typedef treedec::thorup<cfg_t> alg_B;

int main()
{

	cfg_t g(3);
	treedec::add_edge(0, 1, g);
	treedec::add_edge(1, 2, g);
	treedec::add_edge(2, 0, g);
	assert(boost::num_edges(g)==6);
	assert(treedec::num_edges(g)==3);

	cfg_t h(g);
	tree_dec_t t;
	
#ifdef USE_GALA
	alg_A A(g);
	A.do_it(1);
#endif

	alg_B B(h);
	B.do_it();
	// auto t1=B.get_tree_decomposition();
	B.get_tree_decomposition(t);

//	std::cout << A.get_treewidth();

}
