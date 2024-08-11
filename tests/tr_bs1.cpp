
#include <gala/cbset.h>
#include <gala/graph.h>
#include <gala/boost.h>

typedef cbset::BSET_DYNAMIC<2, uint64_t, cbset::nohowmany_t, cbset::nooffset_t, cbset::nosize_t> myset;
template<class A, class...>
using myset_=myset;
typedef gala::graph<myset_, std::vector, unsigned> graph_t;

#include <treedec/tr_bs.hpp>

int main(int argc, char** argv)
{

	unsigned a = 60;
	unsigned b = 70;

	BlockSieve<myset, myset, 4> bs;
	bs.set_n(90);

	myset key;
	for(unsigned i = a; i<b; ++i){
		//trace1("put", key);
		bs.put(key, key);
		key.insert(i);
	}

	return 0;
	key.clear();
	for(unsigned i = 0; i<b; ++i){
		//trace1("put2", key);
		bs.put(key, key);
		key.insert(i);
	}

	key.clear();
	for(unsigned i = 3; i<25; ++i){
		key.insert(i);
	}
	key.insert(45);
	myset s;

	std::vector<myset const*> buf;
	bs.collectSuperblocks( 99, 99, key, s, buf);
	trace1("key", key);
	for(auto i:buf){
		trace1("val ", *i);
	}
	assert(buf.size()==4);

	return 0;
}
