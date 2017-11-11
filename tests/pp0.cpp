#include <boost/graph/adjacency_list.hpp>
#include <treedec/preprocessing.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> G;


int main()
{

	G g;
	treedec::impl::preprocessing<G> A(g);
	A.do_it();

	return 0;
}
