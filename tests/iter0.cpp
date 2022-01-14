#include <set>
#include <iostream>
#include <stdlib.h> // malloc
#include <string.h> // memcpy
#include <vector>
#include <assert.h>
#include <tr1/utility>

#include <boost/typeof/typeof.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edge_list.hpp>
#include <iostream>

#include <treedec/iter.hpp>
#include <treedec/misc.hpp>
#include <treedec/graph_util.hpp>


int main()
{
	boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> g;

	auto v1 = boost::add_vertex(g);
	auto v2 = boost::add_vertex(g);
	auto v3 = boost::add_vertex(g);
	auto v4 = boost::add_vertex(g);

	std::cout << v1 << ": " <<boost::get(boost::vertex_index, g, v1) << "\n";
	std::cout << v2 << ": " <<boost::get(boost::vertex_index, g, v2) << "\n";
	std::cout << v3 << ": " <<boost::get(boost::vertex_index, g, v3) << "\n";

	{
		auto R = make_components_range(g);
		auto& cmpi = R.first;
		auto ii = (*cmpi).first;
		std::cout << *ii << "\n";
		std::cout << "\n";
	}
	std::cout << "========\n";

	boost::add_edge(v1,v2,g);
	{
		std::cout << "\n";

		auto vir=boost::vertices(g);
		std::vector<BOOL> visited(boost::num_vertices(g), false);

		auto M = treedec::util::make_incidence_mask(visited);

		auto R = make_components_range(vir.first, vir.second, g, M);
		auto& cmpi = R.first;

		unsigned k=0;
		for(; cmpi!=R.second;){
			auto P=(*cmpi);
			auto ii=P.first;
			auto ee=P.second;
			for(; ii!=ee; ++ii){
				assert(k==*ii);
				++k;
			}

			++cmpi;
		}
		std::cout << k << "\n";
		std::cout << "========\n";

		{
			auto R = make_components_range(g);
			auto& cmpi = R.first;

			unsigned k=0;
			for(; cmpi!=R.second;){
				auto P = (*cmpi);
				auto ii = P.first;
				auto ee = P.second;
				for(; ii != ee; ++ii){
					assert(k==*ii);
					++k;
				}

				++cmpi;
			}
		}
		std::cout << "========\n";

	}
	{
		boost::add_edge(v3,v4,g);
		// 0 -- 1
		// 2 -- 3

		assert(boost::num_edges(g) == 2);
		assert(boost::num_vertices(g) == 4);
		auto vi=boost::vertices(g).first;
		auto ve=boost::vertices(g).second;
		++vi;
		++vi;
		++vi;
		std::vector<BOOL> visited(boost::num_vertices(g), false);
		auto vm = treedec::util::make_incidence_mask(visited);

		auto R = make_components_range(vi, ve, g, vm);
		auto& cmpi = R.first;
		int k = 0;

		for(; cmpi!=R.second;){
			auto P = (*cmpi);
			auto ii = P.first;
			auto ee = P.second;
			for(; ii != ee; ++ii){
				assert(int(*ii)==3-k);
				++k;
			}

			++cmpi;
		}
		assert(k == 2);
	}
	std::cout << "========\n";
	{
		auto vi = boost::vertices(g).first;
		auto ve = boost::vertices(g).first;
		++vi;
		++ve;
		++ve;
		++ve;
		std::vector<BOOL> visited(boost::num_vertices(g), false);
		auto vm = treedec::util::make_incidence_mask(visited);

		auto R = make_components_range(vi,ve,g,vm);
		//	auto R=make_components_range(g);
		auto& cmpi = R.first;

		for(; cmpi!=R.second;){
			auto P = (*cmpi);
			auto ii = P.first;
			auto ee = P.second;
			for(; ii != ee; ++ii){
				std::cout << " " << *ii;
			}
			std::cout << "\n";

			++cmpi;
		}
	}

	std::cout << "========\n";
	{
		auto v5 = boost::add_vertex(g);
		auto v6 = boost::add_vertex(g);
		auto v7 = boost::add_vertex(g);
		auto v8 = boost::add_vertex(g);

		auto vi=boost::vertices(g).first;
		auto ve=boost::vertices(g).first;
		++vi;
		++ve;
		++ve;
		++ve;
		std::vector<BOOL> visited(boost::num_vertices(g), false);
		visited[0]=true;
		visited[1]=true;
		boost::add_edge(v1,v5,g);
		boost::add_edge(v2,v6,g);
		boost::add_edge(v3,v6,g);
		boost::add_edge(v3,v7,g);
		boost::add_edge(v3,v8,g);

		auto vm=treedec::util::make_incidence_mask(visited);
		auto R=make_components_range(vi,ve,g,vm);
		//	auto R=make_components_range(g);
		auto& cmpi = R.first;

		for(; cmpi!=R.second;){
			auto P=(*cmpi);
			auto ii=P.first;
			auto ee=P.second;
			for(; ii!=ee; ++ii){
				std::cout << " " << *ii;
			}
			std::cout << "\n";

			++cmpi;
		}
	}
	std::cout << "========\n";
	{
		g.clear();
		auto v0=boost::add_vertex(g);
		(void) v0;
		auto v1 = boost::add_vertex(g);
		auto v2 = boost::add_vertex(g);
		auto v3 = boost::add_vertex(g);
		auto v4 = boost::add_vertex(g);
		auto v5 = boost::add_vertex(g);
		auto v6 = boost::add_vertex(g);
		auto v7 = boost::add_vertex(g);
		auto v8 = boost::add_vertex(g);
		(void) v8;

		auto vi=boost::vertices(g).first;
		auto ve=boost::vertices(g).first;
		++vi;
		++ve;
		++ve;
		++ve;
		std::vector<BOOL> visited(boost::num_vertices(g), false);
		boost::add_edge(v1,v3,g);
		boost::add_edge(v2,v4,g);
		boost::add_edge(v3,v5,g);
		boost::add_edge(v4,v6,g);
		boost::add_edge(v5,v7,g);

		auto vm = treedec::util::make_incidence_mask(visited);
		auto R = make_components_range(vi,ve,g,vm);
		auto& cmpi = R.first;

		unsigned c = 0;
		for(; cmpi!=R.second; ++cmpi){
			auto P = (*cmpi);
			auto ii = P.first;
			auto ee = P.second;
			unsigned m = 0;
			for(; ii!=ee; ++ii){
				assert (2*m+c+1 == *ii);
				++m;
			}

			++c;
		}
	}
	std::cout << "========\n";
	{
		g.clear();
		auto v0 = boost::add_vertex(g);
		(void) v0;
		auto v1 = boost::add_vertex(g);
		auto v2 = boost::add_vertex(g);
		auto v3 = boost::add_vertex(g);
		auto v4 = boost::add_vertex(g);
		auto v5 = boost::add_vertex(g);
		auto v6 = boost::add_vertex(g);
		auto v7 = boost::add_vertex(g);
		auto v8 = boost::add_vertex(g);
		(void) v8;

		auto vi = boost::vertices(g).first;
		auto ve = boost::vertices(g).first;
		++vi;
		++ve;
		++ve;
		++ve;
		std::vector<BOOL> visited(boost::num_vertices(g), false);

		// 0  1 -- 3 -- 5 -- 7
		//    2 -- 4 -- 6
		//
		// [vi, ve) = {1,2}
		//
		boost::add_edge(v1,v3,g);
		boost::add_edge(v2,v4,g);
		boost::add_edge(v3,v5,g);
		boost::add_edge(v4,v6,g);
		boost::add_edge(v5,v7,g);
		visited[0]=true;
		visited[2]=true;
		visited[7]=true;
		// visit 1 3 5...


		auto vm = treedec::util::make_incidence_mask(visited);
		auto R = make_components_range(vi, ve, g, vm);
		auto& cmpi = R.first;

		unsigned k = 0;
		for(; cmpi!=R.second;){
			unsigned kk = 0;
			auto P = (*cmpi);
			auto ii = P.first;
			auto ee = P.second;
			for(; ii!=ee; ++ii){
				assert(k==0);
				trace2("visit", *ii, k);
				assert(*ii == 2* kk + 1);
				++kk;
			}
			std::cout << "\n";

			++k;
			++cmpi;
		}
	}
}
