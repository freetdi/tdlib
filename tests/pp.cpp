#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <treedec/preprocessing.hpp>
#include <treedec/graph.hpp>
#include <boost/graph/copy.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> ALSVU;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> ALVVU;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> ALSVD;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> ALVVD;

#ifdef HAVE_GALA_GRAPH_H
// ... later
#endif

typedef boost::tuple< unsigned, std::set<unsigned> > mypair_type;
typedef std::vector< mypair_type  > bags_type;

int main()
{
	ALVVD alvvd1(3);
	treedec::add_edge(0, 1, alvvd1);
	treedec::add_edge(1, 2, alvvd1);
	treedec::add_edge(2, 0, alvvd1);
	assert(boost::num_edges(alvvd1)==6);
	assert(treedec::num_edges(alvvd1)==3);

	ALSVU alsvu1;
	boost::copy_graph(alvvd1, alsvu1);
	assert(boost::num_edges(alsvu1)==3);

	// legacy interface (does it make sense?)
	bags_type bags;
	treedec::preprocessing(alsvu1, bags);
	assert(bags.size()==2);
	bags.clear();
	treedec::preprocessing(alvvd1, bags);
	assert(bags.size()==2);
}
