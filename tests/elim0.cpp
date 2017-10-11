

#include <tdlib/combinations.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_t;

int main(){

    graph_t G;
	 typename treedec::graph_traits<graph_t>::treedec_type T;
    treedec::PP_MD(G, T);
}

