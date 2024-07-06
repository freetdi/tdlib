
#include <gala/cbset.h>
#include <gala/graph.h>
#include <gala/boost.h>

// typedef unsigned __int128 uint128_t; // GCC

typedef cbset::BSET_DYNAMIC<2, uint64_t, cbset::nohowmany_t, cbset::nooffset_t, cbset::nosize_t> myset;
template<class A, class...>
using myset_=myset;
typedef gala::graph<myset_, std::vector, unsigned> graph_t;

#include "tr_iodec.h"
#include "tr_bag.h"

int main(int argc, char** argv)
{
	int n = 5;

	if(argc>1){
		n = atoi(argv[1]);
	}

	graph_t g(n*n);

	std::vector<std::pair<int, int>> edges;
	// making edges of graph representing grid n by n
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < n; y++) {
			if (x < n - 1){
				boost::add_edge(x + y * n, x + 1 + y * n, g);
			}else{
			}
			if (y < n - 1){
				boost::add_edge(x + y * n, x + (y + 1) * n, g);
			}else{
			}
		}
	}

	auto i = boost::edges(g);
	auto k=0;
	for(; i.first!=i.second; ++i.first){
		++k;
		trace2("", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == 2 * (n-1) * n);

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){
//		trace1("grid", *j.first);
//	}

	auto vi = boost::vertices(g);
	for(; vi.first!=vi.second; ++vi.first){
		unsigned d = boost::degree(*vi.first, g);
		trace2("grid", *vi.first, d);
	}

	TR<graph_t> a(g);

	a.do_it();

//	Tree t;
	// Graph g(edges.begin(), edges.end(), n * n);
	// treedec::exact_decomposition_cutset(g, t);

	// std::cout << treedec::get_bagsize(t) << "\n";
	return 0;
}
