#include <tdlib/graph_traits.hpp>
#include <tdlib/elimination_orderings.hpp>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;

#include <tdlib/treedec.hpp>

int main()
{
	using namespace boost;
	G g(10);

	boost::add_edge(0,1,g);
	boost::add_edge(0,2,g);
	boost::add_edge(1,2,g);
	boost::add_edge(3,4,g);
	boost::add_edge(3,7,g);

	boost::add_edge(8,9,g);
	boost::add_edge(8,5,g);

	unsigned j=0;
	std::vector<size_t> o(num_vertices(g));

	BOOST_AUTO(i, o.begin());
	for(;i!=o.end(); ++i){
		*i = j++;
	}
	std::vector<size_t> io(o.size());

	typename treedec::graph_traits<G>::treedec_type t;

	treedec::draft::vec_ordering_to_tree(g, o, t, &io);
	assert(!treedec::check_treedec(g, t));

	std::cerr << "has bagsize " << treedec::get_width(t) << "\n";

	// check if the inverse permutation is indeed one
	// needed??
	j = 0;
	BOOST_AUTO(ii, io.begin());
	for(; ii!=io.end(); ++ii){
		if(*ii!=j++){
			std::cerr << "not an inverse perm" << *ii << " " << j << "\n";
		}
	}
}
