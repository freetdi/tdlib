
#define DEBUG_FILL

#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include "graph_traits.hpp"
// #include "treedec_traits.hpp"

#include <gala/cbset.h>
#include <gala/td.h>
#include <gala/boost.h>
#include <treedec/directed_view.hpp>


#include "induced_supergraph.hpp"
#include "elimination_orderings.hpp"
#include <gala/graph.h>

#include <stdarg.h>
#include <stdio.h>

#define assert_symmetric(g) { \
    unsigned i=0; \
    auto E=boost::edges(g); \
    for(;E.first!=E.second; ++E.first){ \
        ++i;\
        assert( boost::edge(boost::source(*E.first, g), \
                    boost::target(*E.first, g),g).second); \
        assert( boost::edge(boost::target(*E.first, g), \
                    boost::source(*E.first, g),g).second); \
    }\
    trace1("symmetric", i); \
}

template<class G>
struct uvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=false; // TODO: add "false" test.
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    // static constexpr bool force_symmetric=true; // meaninngless (undirected)
    // typedef tdDEGS<G> degs_type; // obsolete.
};
//typedef gala::graph<std::vector, std::vector, unsigned, uvv_config> graph_t;
typedef gala::graph<std::vector, std::vector, unsigned> graph_t;

template<class X, class ... rest>
struct algo_config : treedec::algo::default_config<X, rest...>{
    static void message(int badness, const char* fmt, ...) {
        if (badness >= bLOG){
            char buffer[2048] = "c ";
            va_list arg_ptr;
            va_start(arg_ptr,fmt);
            vsprintf(buffer+2,fmt,arg_ptr);
            va_end(arg_ptr);
            std::cout << buffer;
        }else{
        }
    }
};
typedef treedec::pending::impl::fillIn<graph_t, algo_config> algo_type;

template<class A, class G>
void biedge(A a, A b, G& g)
{
	boost::add_edge(a,b,g);
	boost::add_edge(b,a,g);
}

void grid(int n)
{
	std::cout << "grid " << n << " fi...";
	std::cout << std::flush;
	graph_t g(n*n);

	std::vector<std::pair<int, int>> edges;
	// making edges of graph representing grid n by n
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < n; y++) {
			if (x < n - 1){
				biedge(x + y * n, x + 1 + y * n, g);
			}else{
			}
			if (y < n - 1){
				biedge(x + y * n, x + (y + 1) * n, g);
			}else{
			}
		}
	}
	assert_symmetric(g);

	auto i = boost::edges(g);
	auto k=0;
	for(; i.first!=i.second; ++i.first){
		++k;
		trace2("", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == 2 * (n-1) * n);

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){ untested();
//		trace1("grid", *j.first);
//	}

	auto vi = boost::vertices(g);
	for(; vi.first!=vi.second; ++vi.first){
		unsigned d = boost::degree(*vi.first, g);
		trace2("grid", *vi.first, d);
	}

	algo_type a(g);

	std::cout << ".";
	a.do_it();
	std::cout << " " << a.get_bagsize() << "\n";

//	Tree t;
	// Graph g(edges.begin(), edges.end(), n * n);
	// treedec::exact_decomposition_cutset(g, t);

	// std::cout << treedec::get_bagsize(t) << "\n";
}

void cycle(int n)
{
	graph_t g(n*n);

	std::vector<std::pair<int, int>> edges;
	// making edges of graph representing grid n by n
	for (int x = 1; x <= n; x++) {
		biedge(x%n, x-1, g);
	}
	assert_symmetric(g);

	auto i = boost::edges(g);
	auto k=0;
	for(; i.first!=i.second; ++i.first){
		++k;
		trace2("", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == n);

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){ untested();
//		trace1("grid", *j.first);
//	}

	auto vi = boost::vertices(g);
	for(; vi.first!=vi.second; ++vi.first){
		unsigned d = boost::degree(*vi.first, g);
		trace2("grid", *vi.first, d);
	}

	algo_type a(g);

	std::cout << "cycle " << n << " fi...";
	a.do_it();
	std::cout << " " << a.get_bagsize() << "\n";

//	Tree t;
	// Graph g(edges.begin(), edges.end(), n * n);
	// treedec::exact_decomposition_cutset(g, t);

	// std::cout << treedec::get_bagsize(t) << "\n";
}

void clique(int n)
{
	std::cout << "clique " << n << " fi...";
	std::cout << std::flush;
	graph_t g(n);

	std::vector<std::pair<int, int>> edges;

	for (int x = 0; x < n; x++) {
		for (int y = 0; y < x; y++) {
			biedge(x, y, g);
		}
	}
	assert_symmetric(g);

	auto i = boost::edges(g);
	auto k=0;
	for(; i.first!=i.second; ++i.first){
		++k;
		trace2("", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == (n-1) * n / 2);

	std::cout << "." << std::flush;

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){ untested();
//		trace1("grid", *j.first);
//	}

	algo_type a(g);

	std::cout << ".";
	a.do_it();
	std::cout << a.get_bagsize() << "\n";

//	Tree t;
	// Graph g(edges.begin(), edges.end(), n * n);
	// treedec::exact_decomposition_cutset(g, t);

	// std::cout << treedec::get_bagsize(t) << "\n";
}

void rect(int n, int /*m*/)
{
	graph_t g(n*(n-1));

	std::vector<std::pair<int, int>> edges;
	// making edges of graph representing rect n by n
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < (n-1); y++) {
			if (x < n - 1){
				biedge(x + y * n, x + 1 + y * n, g);
			}else{
			}
			if (y < n - 2){
				biedge(x + y * n, x + (y + 1) * n, g);
			}else{
			}
		}
	}
	assert_symmetric(g);

	auto i = boost::edges(g);
	auto k=0;
	for(; i.first!=i.second; ++i.first){
		++k;
		trace2("grid rect", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == (n) * (n-2) + (n-1) * (n-1));

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){
//		trace1("grid", *j.first);
//	}

	auto vi = boost::vertices(g);
	for(; vi.first!=vi.second; ++vi.first){
		unsigned d = boost::degree(*vi.first, g);
		trace2("grid", *vi.first, d);
	}

	algo_type a(g);

	a.do_it();
	std::cout << "grid rect " << n << " fi " << a.get_bagsize() << "\n";

//	Tree t;
	// Graph g(edges.begin(), edges.end(), n * n);
	// treedec::exact_decomposition_cutset(g, t);

	// std::cout << treedec::get_bagsize(t) << "\n";
}

int main(int, char**)
{
	grid(14);
	grid(10);
	grid(20);
	grid(50);
//	grid(100);
	rect(2, 1);
	rect(3, 2);
	rect(4, 3);
	rect(5, 4);
	rect(6, 5);
	rect(7, 6);
	rect(8, 7);
	rect(10, 9);
	clique(3);
	clique(6);
	clique(10);
	clique(15);
	clique(20);
	grid(3);
	grid(5);
	cycle(5);
	cycle(100);
//
	return 0;
//#ifndef DEBUG_FILL
//	grid(100);
//	clique(100);
//#endif
	//test0(50);
	//return 0;
}
