#include <random>
#include <algorithm>
#include <iterator>

#include "induced_supergraph.hpp"
#include "numbering.hpp"
#include "graph_util.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/graph_utility.hpp>
//#include <boost/random.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> baluvv_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS> baldss_t;

template<class G, class V>
void biedge(V i, V j, G& g)
{
	assert(j!=i);
	boost::add_edge(j, i, g);
	boost::add_edge(i, j, g);
}

void test0(unsigned size, int ne)
{

	std::random_device rd;
	std::mt19937 rng(rd());
	// if(argc>1){
	// 	//rng.seed(atoi(argv[2]));
	// 	size  = atoi(argv[1]);
	// }else{
	// }
	rng.seed(5);
	baluvv_t rg;
	boost::generate_random_graph(rg, size, ne, rng, false, false);

	std::vector<int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

	std::shuffle(v.begin(), v.end(), rng);

	bald_t g(size);
	treedec::draft::NUMBERING_1<bald_t> num(g);
	std::vector<unsigned> degrees(size);

	baldss_t h(size);
	typedef baldss_t::vertex_descriptor vd_ss;
	std::vector<vd_ss> vh(size);
	std::map<vd_ss,int> vm;
	std::vector<unsigned> seq(size);
	auto idx=0;
	for(auto V=boost::vertices(h); V.first!=V.second; ++V.first){
		trace2("setgraph", idx, *V.first);
		vh[idx] = *V.first;
		vm[*V.first] = idx;
		seq[idx] = idx;
		idx++;
	}

	std::shuffle(seq.begin(), seq.end(), rng);

	// boost::copy_graph(g, h); // does not work.
	auto x = boost::edges(rg);
	for(; x.first!=x.second; ++x.first){
		auto s = boost::source(*x.first, rg);
		auto t = boost::target(*x.first, rg);
		biedge(vh[s], vh[t], h);
		biedge(s, t, g);
		trace2("edge", s, t);
	}

	treedec::Supergraph<bald_t, treedec::draft::NUMBERING_1<bald_t>, std::vector<unsigned> > s(g, num, degrees);

	std::vector<std::set<vd_ss>> bag_sets(size/2);

	for(unsigned i=0; i<size/2; ++i) {
		auto t = seq[i];
		auto tvh = boost::adjacent_vertices(vh[t], h);
		int k = 0;
		for(; tvh.first!=tvh.second; ++tvh.first){
			bag_sets[i].insert(*tvh.first);
		}

		std::set<long> chk;
		auto tvs = boost::adjacent_vertices(t, s);
		for(; tvs.first!=tvs.second; ++tvs.first){
			if(t!=*tvs.first){
				trace1("deg", *tvs.first);
				chk.insert(*tvs.first);
				++k;
			}else{
			}
		}

		trace3("deg", k, boost::out_degree(vh[t], h), boost::out_degree(t, s));
		assert( boost::out_degree(vh[t], h) == boost::out_degree(t, s));
		assert( boost::out_degree(vh[t], h) == chk.size());

		trace3("elim", i, t, vh[t]);
//		boost::print_graph(g);
		treedec::eliminate_vertex(t, s);
		treedec::eliminate_vertex(vh[t], h);
		trace3("elimd", i, t, vh[t]);
//		boost::print_graph(g);
		trace2("dbg", boost::num_vertices(h), boost::num_vertices(s));
		assert( boost::num_vertices(h) == boost::num_vertices(s));

		auto E=boost::edges(g);
//		for(; E.first!=E.second; ++E.first){ untested();
//			auto s = boost::source(*E.first, rg);
//			auto t = boost::target(*E.first, rg);
//			trace2("EDG", s, t);
//		}

		{
			auto r = boost::vertices(s);
			for(;r.first!=r.second; ++r.first){
				auto v = *r.first;
				trace1("vertex s", v);
			}
		}

		{
			auto r = boost::vertices(s);
			for(;r.first!=r.second; ++r.first){
				auto v = *r.first;
				auto meh = treedec::count_missing_edges(vh[v], h);
				trace1("===== degcheck", v);
				int k = 0;

				std::set<long> K;
				auto tvs = boost::adjacent_vertices(v, s);
				for(; tvs.first!=tvs.second; ++tvs.first) {
					auto n = *tvs.first;
					trace1(".... degcheck", n);
					assert(!num.is_numbered(n));
					if(n == v){
					}else{
						++k;
						K.insert(n);
					}
				}
				trace4("degcheck", *r.first, k, K.size(), boost::out_degree(*r.first, s));
				trace4("degcheck", *r.first, k, K.size(), boost::out_degree(vh[*r.first], h));
				assert(boost::out_degree(*r.first, s) ==  boost::out_degree(vh[*r.first], h));
				assert(boost::out_degree(*r.first, s) == K.size());
#if 1
				auto mes = treedec::count_missing_edges(*r.first, s);
				trace3("mecheck", meh, mes, v);
				trace1("mecheck", boost::out_degree(v, *s));
				assert(meh==mes);
#endif
			}
		}

	}


	int k = 0;
	auto r = boost::vertices(h);
	for(;r.first!=r.second; ++r.first){
		trace1("check", *r.first);
		++k;
	}

	auto rs = boost::vertices(s);
	for(;rs.first!=rs.second; ++rs.first){
		--k;
	}
	assert(k==0);

	assert(size==seq.size());
	for(unsigned j=size/2; j<size; ++j){
		auto t = seq[j];
		assert(t<size);
		auto dh = boost::out_degree(vh[t], h);
		auto ds = boost::out_degree(t, s);

		trace3("check", t, dh, ds);

		auto tvh = boost::adjacent_vertices(vh[t], h);
		long unsigned k = 0;
		for(; tvh.first!=tvh.second; ++tvh.first){
			trace3("h ++adj", k, vm[*tvh.first], *tvh.first);
			++k;
		}
		assert(k==dh);

		std::set<long> chk;
		auto tvs = boost::adjacent_vertices(t, s);
		for(; tvs.first!=tvs.second; ++tvs.first){
			if(*tvs.first == t){
			}else{
				trace2("--adj", k, *tvs.first);
				chk.insert(*tvs.first);
			}
			assert(!num.is_numbered(*tvs.first));
		}
		trace5("deg", j, t, dh, ds, k);

		assert(k==chk.size());
		assert(dh==ds);
	}

	// bag contents.
	trace1("======== bags", size);
	for(unsigned i=0; i<size/2; ++i) {
		auto ii = seq[i];
		auto bb = s.bag_vertices(ii);
		size_t k = 0;
		for(; bb.first!=bb.second; ++bb.first){
			++k;
		}
		// trace5("bag compare", bag_sets[i].size(), i, ii, k, boost::out_degree(ii ,s));
		assert(k==bag_sets[i].size());
//		assert(k==boost::out_degree(ii ,s)); incomplete()
	}
}

int main(int, char**)
{
	std::cout << "1\n";
	test0(10, 5);
	std::cout << "2\n";
	test0(100, 1000);
	std::cout << "3\n";
	test0(130, 500);
}
