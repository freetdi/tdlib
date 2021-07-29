#include <treedec/graph.hpp>
#include <iostream>

#ifndef TESTGRAPH_T
#include <boost/graph/adjacency_list.hpp>
#define TESTGRAPH_T boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
#endif

int main()
{
	TESTGRAPH_T G(5);

	boost::add_edge(0,1,G);
	boost::add_edge(0,2,G);
	boost::add_edge(0,3,G);
	boost::add_edge(0,4,G);

	boost::add_edge(1,0,G);
	boost::add_edge(1,2,G);
	boost::add_edge(1,4,G);

	BOOST_AUTO(i, treedec::common_out_edges(0,1,G).first);
	BOOST_AUTO(e, treedec::common_out_edges(0,1,G).second);

	unsigned count=0;
	for(;i!=e;++i){
		++count;
		std::cout << *i << "\n";
	}

	assert(count==2);
}
