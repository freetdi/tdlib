#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <tdlib/preprocessing.hpp>
#include <tdlib/graph.hpp>
#include <boost/graph/copy.hpp>
#include <tdlib/combinations.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
typedef treedec::graph_traits<G>::treedec_type T;

typedef treedec::comb::ex17<G, treedec::algo::default_config> alg_A;

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

	A.do_it();

	T t;
	A.get_tree_decomposition(t);

	std::cout << "tw is " << treedec::get_bagsize(t)-1 << "\n";
	assert(3==treedec::get_bagsize(t));
}
