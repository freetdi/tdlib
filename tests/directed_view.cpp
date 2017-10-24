#include <tdlib/directed_view.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <tdlib/graph.hpp> // num_edges
#include <boost/graph/graph_utility.hpp>

struct cfg_node
{
  // iCode *ic; // bnot relevant
  // operand_map_t operands;
  bool alive;
  bool dying;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, cfg_node> cfg_t;
typedef treedec::draft::directed_view<cfg_t> cfgd_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::directedS, cfg_node> w_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::undirectedS, cfg_node> u_t;

template<class G>
void test()
{
	G g;
	auto a=boost::add_vertex(g);
	auto b=boost::add_vertex(g);
	boost::add_edge(a, b, g);
	boost::print_graph(g);

	auto q=boost::adjacent_vertices(a, g);
	for(;q.first!=q.second;++q.first){
		std::cout << a << ";" << *q.first << "\n";
	}
	{
	auto q=boost::adjacent_vertices(b, g);
	for(;q.first!=q.second;++q.first){
		std::cout << b << ";" << *q.first << "\n";
	}
	}

	auto p=boost::edges(g);
	for(;p.first!=p.second;++p.first){
		auto e=*p.first;
		std::cout << boost::source(e, g) << " : " << boost::target(e, g) << "\n";
	}

	{
	typedef treedec::draft::directed_view<G> D;
	D d(g);

	auto p=boost::edges(d);
	for(;p.first!=p.second;++p.first){
		auto e=*p.first;
		std::cout << boost::source(e, g) << " : " << boost::target(e, g) << "\n";
	}
	assert(boost::num_edges(d)==2); // backend sees two edges
	assert(treedec::num_edges(d)==1); // logical edges.

	std::cerr << "outedges\n";
	{
	auto q=boost::adjacent_vertices(a, g);
	for(;q.first!=q.second;++q.first){
		std::cout << a << ":" << *q.first << "\n";
	}
	}
	{
	auto q=boost::adjacent_vertices(b, g);
	for(;q.first!=q.second;++q.first){
		std::cout << b << ":" << *q.first << "\n";
	}
	std::cerr << "/outedges\n";
	}
	}
}

int main()
{
	BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<cfg_t> ));
	BOOST_CONCEPT_ASSERT(( boost::BidirectionalGraphConcept<cfg_t> ));

	// it's not.
	// BOOST_CONCEPT_ASSERT(( boost::BidirectionalGraphConcept<w_t> ));

	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<cfg_t> ));
	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<cfgd_t> ));
	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<u_t> ));
	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<w_t> ));

	BOOST_CONCEPT_ASSERT(( boost::PropertyGraphConcept<cfg_t, unsigned, bool cfg_node::*> ));

	BOOST_STATIC_ASSERT( std::is_same<typename cfgd_t::wrapped_type, w_t>::value);

	BOOST_CONCEPT_ASSERT(( boost::PropertyGraphConcept<cfgd_t, unsigned, bool cfg_node::*> ));

	test<cfg_t>();

}
