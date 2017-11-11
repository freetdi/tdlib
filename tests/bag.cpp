#include <treedec/container.hpp>
#include <treedec/graph_traits.hpp>

typedef
 boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
typedef typename treedec::graph_traits<G>::treedec_type T;

#include <treedec/treedec.hpp>
int main(){

	T t;
	boost::get(treedec::bag_t(), t);

}

