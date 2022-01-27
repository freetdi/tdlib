
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

template<class X, class ... rest>
struct algo_config : treedec::algo::default_config<X, rest...>{
    static void message(int badness, const char* fmt, ...) { itested();
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
typedef typename treedec::graph_traits<graph_t>::treedec_type T;

typedef treedec::pending::impl::fillIn<graph_t, algo_config> algo_type;

template<class A, class G>
void biedge(A a, A b, G& g)
{
	boost::add_edge(a,b,g);
	boost::add_edge(b,a,g);
}

void clique(int n)
{
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
		trace2("clique edg", boost::source(*i.first, g), boost::target(*i.first, g));
	}
	trace2("dbg", k, n);
	assert(k == (n-1) * n / 2);

//	auto j = boost::adjacent_vertices(4, g);
//	for(; j.first!=j.second; ++j.first){
//		trace1("grid", *j.first);
//	}

	auto vi = boost::vertices(g);
	for(; vi.first!=vi.second; ++vi.first){
		unsigned d = boost::degree(*vi.first, g);
		trace2("clique:", *vi.first, d);
	}

	algo_type a(g);

	a.do_it();

	T t;
	a.get_tree_decomposition(t);

	if(treedec::get_bagsize(t) == a.get_bagsize()){
	}else{
		incomplete();
	}
	assert(treedec::is_valid_treedecomposition(g, t));
	std::cout << "done clique " << n << " fi " << a.get_bagsize() << "\n";
}

int main(int, char**)
{
	clique(1);
	clique(2);
	clique(3);
	clique(6);
	clique(10);
	clique(15);
	clique(20);
}
