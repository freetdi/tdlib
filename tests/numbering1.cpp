

#include <boost/property_map/property_map.hpp> 
#include "numbering.hpp"
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/properties.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
typedef treedec::draft::NUMBERING_1<bald_t> Numbering;

void test1()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(8);
	num.increment();

	num.indistinguishable(2, 6);
	num.indistinguishable(3, 6);
	num.indistinguishable(4, 6);
	num.indistinguishable(9, 6);
	num.indistinguishable(7, 6);

	num.put(6);
	sns[6] = 6;
	sns[2] = 0;
	sns[3] = 0;
	sns[4] = 0;
	sns[7] = 0;
	sns[9] = 0;
	num.increment(sns[6]);

	num.put(0);
	num.increment();
	num.put(1);
	num.increment();
	num.put(5);
	num.increment();

	std::vector<long> ord(10);

	num.get_ordering(ord, sns, Numbering::oTree);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << " " << sns[ord[i]] << "\n";
	}

	assert(ord[0]==8);
	assert(ord[1]==6);   assert(sns[6]==6);
	  assert(ord[2]==2);
	  assert(ord[3]==3);
	  assert(ord[4]==4);
	  assert(ord[5]==7);
	  assert(ord[6]==9);
	
	assert(ord[7]==0);
	assert(ord[8]==1);
	assert(ord[9]==5);
}

// a chain. 8 - 2 - 4 - 6
void test2()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.indistinguishable(6, 2);
	num.indistinguishable(4, 2);
	num.indistinguishable(2, 8);
	sns[8] = 4;
	sns[2] = -2;
	sns[6] = 0;
	sns[4] = 0;
	num.put(8);
	num.increment(4);

	num.put(0);
	num.increment();
	num.put(1);
	num.increment();
	num.put(3);
	num.increment();
	num.put(7);
	num.increment();
	num.put(5);
	num.increment();
	num.put(9);
	num.increment();

	std::vector<long> ord(10);

	num.get_ordering(ord, sns, Numbering::oTree);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": " << ord[i] << " " << sns[ord[i]] <<"\n";
	}

	assert(ord[0]==8);
	assert(ord[1]==2); assert(sns[2]==-2);
	assert(ord[4]==0);
	assert(sns[8]==4); // ??
}
// 2
 /// // 4
 /// //  3
void test3()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(1);
	num.increment();

	num.indistinguishable(4, 2);
	num.indistinguishable(3, 4);
	sns[4] = -1;
	sns[2] = 3;
	sns[3] = 0;
	num.put(2);
	num.increment(sns[2]);

	num.put(0);
	num.increment();
	num.put(5);
	num.increment();
	num.put(6);
	num.increment();
	num.put(7);
	num.increment();
	num.put(8);
	num.increment();
	num.put(9);
	num.increment();

	std::vector<long> ord(10);

	for(unsigned i=0; i<N; ++i){
		std::cout << "input " << i << ": sns " << sns[i] << " pos " << num.get_position(i) << "\n";
	}

	num.get_ordering(ord, sns, Numbering::oTree);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": ord " << ord[i] << " " << sns[ord[i]] << "\n";
	}
	assert(ord[0]==1); assert(sns[1]==1);
	assert(ord[1]==2); assert(sns[2]==3);
}
// 2
 /// // 3
 /// //  4
 /// //  5
// 1
// 0
// 7->3.
// 8->7
void test4()
{
	const unsigned N=10;
	bald_t g(N);
	treedec::draft::NUMBERING_1<bald_t> num(g);

	std::vector<long> sns(10, 1);

	num.put(2);
	num.increment();

	//num.put(3);
	num.indistinguishable(4, 3); // 4->3
	num.indistinguishable(5, 3); // 5->3
	sns[4] = 0;
	sns[5] = 0;
	// num.increment(sns[3]);

	num.put(1);
	num.increment();
	num.put(0);
	num.increment();

	num.indistinguishable(3, 7); // 3->7
	sns[3] = -2; // 3 has subtree of size 2, but is not root.
	sns[7] = 4;  // 7 has them all.
//	num.put(8);

	num.indistinguishable(8, 7); // 8->7
	++sns[7];
	sns[8] = 0;
	num.put(7);
	num.increment(sns[7]);

	num.put(6);
	num.increment();


	num.put(9);
	num.increment();

	std::vector<long> ord(10);

	for(unsigned i=0; i<N; ++i){
		std::cout << "input " << i << ": sns " << sns[i] << " pos " << num.get_position(i) << "\n";
	}

	num.get_ordering(ord, sns, Numbering::oTree);

	for(unsigned i=0; i<N; ++i){
		std::cout << i << ": ord " << ord[i] << " " << sns[ord[i]] << "\n";
	}

	assert(ord[8]==6);
}

int main()
{
	std::cout << "----1\n";
	test1();
	std::cout << "----2\n";
	test2();
	std::cout << "----3\n";
   test3();
	std::cout << "----4\n";
   test4();
	return 0;
}
