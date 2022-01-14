

#include <boost/property_map/property_map.hpp> 
#include "numbering.hpp"
#include "induced_supergraph.hpp"
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/properties.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
typedef treedec::draft::NUMBERING_1<bald_t> Numbering;

struct visitor1{
	template<class V>
	bool operator()(V v){
		std::cout << v << "\n";
		return true;
	}
};

void test1()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(0);
	num.increment();

	num.indistinguishable(3, 4);

	num.put(1);
	num.increment();

	num.indistinguishable(6, 5);
	num.indistinguishable(9, 5);

	num.put(7);
	num.increment();


	num.indistinguishable(5, 2);
	num.indistinguishable(4, 2);
	num.put(2);

	num.increment(6); /// 2 3 4 5 6 9


	num.put(8);
	num.increment(1);

	sns[0] = 1;
	sns[1] = 1;
	sns[2] = 5;
	sns[3] = 0;
	sns[4] = -1;
	sns[5] = -2;
	sns[6] = 0;
	sns[7] = 1;
	sns[8] = 1;
	sns[9] = 0;
	std::vector<long> ord(10);

	num.get_ordering(ord, sns, Numbering::oTree);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << " " << sns[ord[i]] << "\n";
	}

	visitor1 v1;
	std::cout << "--1\n";
	treedec::draft::visit_supernode(1, num, ord, sns, v1);
	std::cout << "--2\n";
	treedec::draft::visit_supernode(2, num, ord, sns, v1);
	std::cout << "--3\n";
	treedec::draft::visit_supernode(3, num, ord, sns, v1);
	std::cout << "--4\n";
	treedec::draft::visit_supernode(4, num, ord, sns, v1);
	std::cout << "--5\n";
	treedec::draft::visit_supernode(5, num, ord, sns, v1);
}

int main()
{
	std::cout << "----1\n";
	test1();
	return 0;
}
