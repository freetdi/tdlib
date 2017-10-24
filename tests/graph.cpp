#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>

struct foo{
	int a;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, foo> Gb;

BOOST_STATIC_ASSERT(std::is_convertible<
		typename boost::graph_traits<Gb>::directed_category*,
		boost::directed_tag* >::value);
BOOST_STATIC_ASSERT(std::is_convertible<
		typename boost::graph_traits<Gb>::traversal_category*,
		boost::bidirectional_graph_tag* >::value);

int main(){
}
