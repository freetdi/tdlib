#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <tdlib/preprocessing.hpp>
#include <tdlib/graph.hpp>
#include <boost/graph/copy.hpp>
#include <tdlib/exact_ta.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> ALSVU;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> ALVVU;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> ALSVD;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> ALVVD;

int main()
{
	return 0; //for now
	ALVVD alvvd1(3);
	treedec::add_edge(0, 1, alvvd1);
	treedec::add_edge(1, 2, alvvd1);
	treedec::add_edge(2, 0, alvvd1);
	assert(boost::num_edges(alvvd1)==6);
	assert(treedec::num_edges(alvvd1)==3);

	ALSVU alsvu1;
	boost::copy_graph(alvvd1, alsvu1);
	assert(boost::num_edges(alsvu1)==3);

	treedec::exact_ta<ALSVU> A(alsvu1);
	A.do_it(1);

//	std::cout << A.get_treewidth();
}
