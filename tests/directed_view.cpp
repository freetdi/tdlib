#include <tdlib/directed_view.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>

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

int main()
{
	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<cfg_t> ));
	BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<cfgd_t> ));

	BOOST_CONCEPT_ASSERT(( boost::PropertyGraphConcept<cfg_t, unsigned, bool cfg_node::*> ));

	BOOST_STATIC_ASSERT( std::is_same<typename cfgd_t::wrapped_type, w_t>::value);

	BOOST_CONCEPT_ASSERT(( boost::PropertyGraphConcept<cfgd_t, unsigned, bool cfg_node::*> ));

}
