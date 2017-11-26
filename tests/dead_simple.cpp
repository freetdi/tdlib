#include <treedec/thorup.hpp>
#include <treedec/nice_decomposition.hpp> // nicify
#include <treedec/misc.hpp> // is_valid_td

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::undirectedS> G;
typedef typename treedec::graph_traits<G>::treedec_type T;
typedef treedec::thorup<G> thorupAlgorithm;

int main()
{
	const unsigned n=16;
	G g(n);
	for(unsigned i=0; i<n; ++i){
		boost::add_edge(i, (i+1)%n, g);
	}

	auto const& h(g);
	thorupAlgorithm a(h);
	a.do_it();

	T t;
	a.get_tree_decomposition(t);

	incomplete();
	// treedec::nice::nicify(t); ouch, segfault.

	assert(treedec::is_valid_treedecomposition(g,t));
	auto w=treedec::get_width(t);
	assert(w==2);
	std::cout << "cycle has treewidth <= " << w << "\n";
}
