#include <boost/graph/adjacency_list.hpp>
#include <treedec/preprocessing.hpp>
#include <treedec/treedec.hpp>
#include <treedec/misc.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G; // need undirected?
typedef typename treedec::graph_traits<G>::treedec_type T;


int main()
{
	int n = 10;

	G g(n*n);
	for (int x = 0; x < n; x++) {
		for (int y = 0; y < n; y++) {
			if (x < n - 1){
				boost::add_edge(x + y * n, x + 1 + y * n, g);
			}else{
			}

			if (y < n - 1) {
				boost::add_edge(x + y * n, x + (y + 1) * n, g);
			}else{
			}
		}
	}
	treedec::impl::preprocessing<G> A(g);
	A.do_it();
	T t;
	A.get_tree_decomposition(t);
	trace1("dbg", boost::num_vertices(t));
	auto tb = boost::vertices(t);
	for(; tb.first!=tb.second; ++tb.first){
		trace1("bag", *tb.first);
		auto bb = boost::get(treedec::bag_t(), t, *tb.first);
		for(auto x : bb){
			trace2("bag", *tb.first, x);
		}
	}

	auto ee = boost::edges(t);
	for(; ee.first!=ee.second; ++ee.first){
		auto src = boost::source(*ee.first, t);
		auto tgt = boost::target(*ee.first, t);
		trace2("t edg", src, tgt);
	}

	trace2("treedec?", boost::num_vertices(t), boost::num_edges(t));
	int chk = treedec::check_treedec(g, t);
	std::cout << "check " << chk << "\n";

	return chk;
}
