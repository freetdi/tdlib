#include <boost/graph/adjacency_list.hpp>
#include <treedec/preprocessing.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> G;
typedef typename treedec::graph_traits<G>::treedec_type T;

int main()
{

	G g(10);
	treedec::impl::preprocessing<G> A(g);
	A.do_it();
	T t;
	A.get_tree_decomposition(t);

	int ret = treedec::check_treedec(g, t);
	assert(ret==0);
	return ret;
}
