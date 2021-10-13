

#include <boost/property_map/property_map.hpp> 
#include "numbering.hpp"
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/properties.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;

void test0()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	for(unsigned i=0; i<N; ++i){
		num.put((i+3)%10);
		num.increment();
	}

	std::vector<long> ord(10);
	std::vector<long> sns(10, 1);

	num.get_ordering(ord, sns);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << "\n";
	}
}

void test1()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	for(unsigned i=0; i<4; ++i){
		num.indistinguishable(i, i+5);
	}
	num.indistinguishable(9, 4);

	std::vector<long> sns(10, 0);

	num.put(5);
	num.increment(2);
	sns[5] = 2;

	num.put(6);
	sns[6] = 2;
	num.increment(2);

	num.put(7);
	sns[7] = 2;
	num.increment(2);

	num.put(8);
	sns[8] = 2;
	num.increment(2);

	num.put(4);
	sns[4] = 2;
	num.increment(2);

	std::vector<long> ord(10);

	num.get_ordering(ord, sns);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << "\n";
	}
}

// chain 6 - 2 - 3 - 9 - 4
void test2()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(6);
	num.increment();
	num.put(3);
	num.indistinguishable(2, 6);
	num.indistinguishable(4, 9);
	num.indistinguishable(3, 2);
	num.indistinguishable(9, 3);
	sns[2] = 0;
	sns[3] = 0;
	sns[4] = 0;
	sns[9] = 0;
	sns[6] = 5;
	num.increment(4);

	num.put(8);
	num.indistinguishable(7, 8);
	sns[8] = 2;
	sns[7] = 0;
	num.increment(2);

	num.put(1);
	num.increment();
	num.put(5);
	num.increment();
	num.put(0);
	num.increment();

	std::vector<long> ord(10);

	num.get_ordering(ord, sns);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << "\n";
	}

	assert(ord[0] == 6);
	assert(ord[1] == 2);
	assert(ord[2] == 3);
	assert(ord[3] == 4);
	assert(ord[4] == 9);

}

void test3()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(0);
	num.indistinguishable(1, 0);
	num.indistinguishable(2, 1);
	num.increment(3);

	sns[0] = 2;
	sns[1] = 0;
	sns[2] = 0;

	num.put(3);
	num.increment();
	num.put(4);
	num.increment();
	num.put(6);
	num.increment();
	num.put(5);
	num.increment();
	num.put(7);
	num.increment();
	num.put(8);
	num.increment();
	num.put(9);
	num.increment();

	std::vector<long> ord(10);

	num.get_ordering(ord, sns);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << " sns " << sns[ord[i]] << "\n";
	}
}

int main()
{
	test0();
	std::cout << "----1\n";
	test1();
	std::cout << "----2\n";
	test2();
	std::cout << "----3\n";
	test3();
	return 0;
}
