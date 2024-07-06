
#include <gala/cbset.h>
#include <gala/graph.h>
#include <gala/boost.h>

#define NC 1

typedef cbset::BSET_DYNAMIC<NC, uint64_t, cbset::nohowmany_t, cbset::nooffset_t, cbset::nosize_t> myset;
template<class A, class...>
using myset_=myset;
typedef gala::graph<myset_, std::vector, unsigned> graph_t;

#include <treedec/tr_bs.hpp>

struct Mysetless {
	bool operator()(myset const& a, myset const& b) const {
		return a.compare_int(b) < 0;
	}
}setless;

int main(int, char**)
{
	unsigned M = 1000;
	unsigned N = 64*NC;

	typedef BlockSieve<myset, myset, 128> bs_t;
	bs_t bs;
	bs.set_n(N);

	std::vector<myset> U(M);
	std::set<myset, Mysetless> S;

	int x=0;
	for(unsigned i=0; i<100000; ++i){
		x = (x*319901 + 101) % 2623333;
		U[(i)%M].insert( x % N );
	}

	for(auto key : U){
		std::cout << "ins " << key << "\n";
		bs[key] = key;
		S.insert(key);
	}
	std::cout << "bs size " << bs.size() << "\n";

	assert(S.size() == bs.size());


	myset key;
	myset s;
	for(unsigned i=0; i<25; ++i){
		bs_t::set_out_hack buf;
		bs.collectSuperblocks( 999, 999, key, s, buf);
		std::cout << "superblocks " << i << ": " << key.size() << " " << bs.size() << "\n";
		key.insert((i*17) % N);
	}

	size_t k=0;
	for(auto x : bs){
		(void)x;
		++k;
	}
	assert(bs.size()==k);
	return 0;
#if 0

	trace1("key", key);
	for(auto i:buf){
		trace1("val ", *i);
	}
	assert(buf.size()==4);

	return 0;
#endif
}
