#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <tdlib/preprocessing.hpp>
#include <tdlib/graph.hpp>
#include <boost/graph/copy.hpp>
#include <tdlib/exact_ta.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;

namespace choice{
  template<class G,
		template<class G_, class ...> class C>
  using exact_ta_=treedec::exact_ta<G, C>;
}

typedef treedec::draft::exact_decomposition<G,
             treedec::algo::default_config,
             choice::exact_ta_> alg_A;

int main()
{
	G g(3);
	treedec::add_edge(0, 1, g);
	treedec::add_edge(1, 2, g);
	treedec::add_edge(2, 0, g);
	assert(boost::num_edges(g)==3);
	assert(treedec::num_edges(g)==3);

	alg_A A(g);

	G alsvu1;
	boost::copy_graph(g, alsvu1);
	assert(boost::num_edges(alsvu1)==3);

	A.do_it(1);

//	std::cout << A.get_treewidth();
}
