#include "induced_supergraph.hpp"
#include "numbering.hpp"
#include <boost/graph/adjacency_list.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;

void test0()
{
	bald_t g(10);
	treedec::draft::NUMBERING_1<bald_t> num(g);
	std::vector<unsigned> degrees(10);

	for(unsigned i=0; i<9; ++i){
		boost::add_edge(i, i+1, g);
		boost::add_edge(i+1, i, g);
		num.put(i);
		num.increment();
		degrees[i] = i;
	}

	treedec::Supergraph<bald_t, treedec::draft::NUMBERING_1<bald_t>, std::vector<unsigned> > s(g, num, degrees);

	assert(s.bagsize(5) == 6);

	for(unsigned i=0; i<5; ++i){
		auto R = boost::adjacent_vertices(i, g);
		std::cout << "v " << i << " pos " << num.get_position(i) << "\n";
		for(; R.first!=R.second; ++R.first){
			std::cout << "g " << i << " -> " << *R.first << "\n";
		}
	}
	
	for(unsigned i=0; i<5; ++i){
		auto R = s.bag_vertices(i);
		std::cout << "bag " << i << "\n";
		for(; R.first!=R.second; ++R.first){
			std::cout << " -- " << *R.first << "\n";
		}
	}
}

void biedge(int i, int j, bald_t& g)
{
		boost::add_edge(j, i, g);
		boost::add_edge(i, j, g);
}

void test1()
{
	std::cout << "======== test 1\n";
	int N = 6;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);
	std::vector<unsigned> degrees(6);

	/*       2--5    after deletion of 0
	 *      /        N(2) = bag(0) + unmarked neighbours of 2
	 *  1--0-3       after deletion of 2.
	 *      \|       ->  2 -> 0 needed to see N(2)
	 *       4       ->  0 -> 2 needed to see bag(0)
	 */

	for(unsigned i=1; i<5; ++i) {
		boost::add_edge(0, i, g);
		boost::add_edge(i, 0, g);
	}

	boost::add_edge(3, 4, g);
	boost::add_edge(4, 3, g);
	boost::add_edge(2, 5, g);
	boost::add_edge(5, 2, g);

	degrees[5] = 1; // TODO
	
	treedec::Supergraph<bald_t, treedec::draft::NUMBERING_1<bald_t>, std::vector<unsigned> > s(g, num, degrees);
	assert(boost::out_degree(5, s) == 1);

	std::cout << "== eliminate 1\n";
	num.put(1);
	num.increment();

//	bag 0: 2, 3, 4
// bag 1: 0
// bag 2: 0, 5
	for(unsigned k=0; k<3; ++k){
		if(k==1) continue;
		std::cout << "adj " << k << "\n";
		auto R = s.adjacent_vertices(k);
		for(; R.first!=R.second; ++R.first){
			std::cout << " " << *R.first;
		}
		std::cout << "\n";
	}

	std::cout << "== eliminate 0\n";
	num.put(0);
	num.increment();
	
//	adj 0: 2, 3, 4
// adj 1: 0
// adj 2: 5
	for(unsigned k=0; k<2; ++k){
		auto R = s.bag_vertices(k);
		std::cout << "adj " << k << "\n";
		for(; R.first!=R.second; ++R.first){
			std::cout << " " << *R.first;
		}
		std::cout << "\n";
	}
	for(unsigned k=2; k<3; ++k){
		auto R = boost::adjacent_vertices(k, s);
		std::cout << "adj " << k << "\n";
		for(; R.first!=R.second; ++R.first){
			std::cout << " " << *R.first;
		}
		std::cout << "\n";
	}
}

void test2()
{
	std::cout << "======== test 2\n";
	int N = 6;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);
	std::vector<unsigned> degrees(6);

	/*   0--3       after deletion of 0. bag(0) = {2,3}
	 *    \         N(2) = bag(0) + {1,4}
	 *     2---4    after deletion of 2.
	 *    /         bag(2) = {0,1,3,4,5}
	 *   1--5       N(4) = bag(2)
	 *              bag(2) = bag(0) + bag{1}
	 */

	biedge(0, 3, g);
	biedge(0, 2, g);
	biedge(2, 4, g);
	biedge(2, 1, g);
	biedge(5, 1, g);

	typedef treedec::Supergraph<bald_t, treedec::draft::NUMBERING_1<bald_t>, std::vector<unsigned> > SG;
	SG s(g, num, degrees);
	typedef typename boost::graph_traits<SG>::adjacency_range adjacency_range;

	num.put(0);
	num.increment();

	std::cout << "B0\n"; // 2,3
	auto R = s.bag_vertices(0);
	int k = 0;
	for(; R.first!=R.second; ++R.first){
		++k;
		assert(*R.first == 2
			||  *R.first == 3);
		std::cout << " -- " << *R.first << "\n";
	}
	assert(k==2);

	k = 0;
	std::cout << "N2\n"; // 3,1,4,2
	{
		adjacency_range R2 = boost::adjacent_vertices(2, s);
		std::cout << *R2.first << "\n";
	}

	adjacency_range R2 = boost::adjacent_vertices(2, s);
	for(; R2.first!=R2.second; ++R2.first){
		++k;
		std::cout << " " << *R2.first << "";
	}
	std::cout << "\n";
	assert(k==4);

	num.put(1);
	num.increment();

	k = 0;
	std::cout << "N2\n"; // 4,3,5,2,2
	auto R3 = boost::adjacent_vertices(2, s);
	for(; R3.first!=R3.second; ++R3.first){
		++k;
		std::cout << " " << *R3.first << "";
	}
	std::cout << "\n";
	assert(k==5);

}

int main()
{
	test0();
	test1();
	test2();
}
