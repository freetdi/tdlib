#include <boost/graph/adjacency_list.hpp>
//#include <treedec/treedec_traits.hpp>
#include <treedec/iter.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_t;

int main(int, char**)
{
	unsigned n = 3;
	graph_t g(n*n);

	for (unsigned x = 0; x < n; x++) {
		for (unsigned y = 0; y < n; y++) {
			if (x < n-1){
				boost::add_edge(x + y*n, x + y*n + 1, g);
			}else{
			}
		}
	}

	{
		auto R = treedec::make_components_range(g);
		assert(R.first!=R.second);

		unsigned k = 0;
		for(; R.first != R.second; ++R.first){
			unsigned kk = 0;
			auto P = *R.first;
			for(; P.first!=P.second; ++P.first){
				/// 0 1 2 ; 3 4 5 ; 6 7 8
				assert(*P.first == kk + n*k);
				++kk;
			}
			assert(kk==n);
			++k;
		}
		trace2("===========", k, n);
		assert(k==n);
	}
	if(0)
	{
		auto R = treedec::make_components_range(g);
		unsigned kk = 0;
		for(; R.first != R.second; ){ untested();
			trace1("partial", kk);
			auto P = *R.first;

			for(unsigned m=0; m<kk; ++m){
				trace3("triangle visit", *P.first, kk, m);
//				assert(*P.first == n*kk + m);
				++P.first;
			}

			//++P.first;
			//trace2("visit", *P.first, kk);
			//assert(*P.first == n*kk + 1);

			// ++P.first;
			// trace2("visit", *P.first, kk);
			// assert(*P.first == n*kk + 2);

			++kk;
			++R.first;
		}
		trace2("cmp", kk, n);
		assert(R.first==R.second);
		assert(kk==n);
	}

	trace0("=====empty graph======");
	{
		graph_t g(n*n);
		auto R = treedec::make_components_range(g);
		int k = 0;
		for(;R.first != R.second; ++R.first){
			++k;
		}
		assert(k == n*n);

	}
	trace0("=====omit contents======");

	if(1){
		auto R = treedec::make_components_range(g);
		for(unsigned i=0; i<n; ++i){
			auto range = *R.first;
			trace2("next comp", i, *range.first);
			assert(*range.first == n*i);

			++R.first;
			trace0("done seek");
		}
		assert(R.first==R.second);
	}
}
