
#ifdef HAVE_GALA
#include <gala/boost.h>
#include <gala/graph.h>
#include <gala/td.h>
#include <gala/trace.h>
#endif

#include <boost/graph/copy.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <stdlib.h>
#include <tdlib/elimination_orderings.hpp>
#include <tdlib/graph.hpp>
#include <tdlib/minimum_degree_ordering.hpp>
#include <tdlib/trace.hpp>
// not in development yet
//#include <tdlib/printer.hpp>


#ifdef HAVE_GALA
template<class G>
struct dvv_config : public gala::graph_cfg_default<G>
{
	static constexpr bool is_directed = true;
};

typedef uint32_t unsignedType;

#if 1
typedef gala::graph<std::vector, std::vector, unsignedType, dvv_config> sg_dvv;
#else
typedef gala::graph<std::set, std::vector, uint16_t, dvv_config> sg_dvv;
#endif

#endif

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> balu_t;

// typedef bald_t parsegr_t;

int main(int argc, char** argv)
{
#ifndef HAVE_GALA
	(void) argc;
	(void) argv;
	return 77;
#else

#if 0
	PARSE* p;
	try{
		p = new PARSE(std::cin);
	}catch(...){
		std::cout << "uuh\n";
		exit(2);
	}
#endif

	size_t size=1<<5;
	if(argc>1) size=1<<atoi(argv[1]);
	size_t ne=size*size/8 + 1;

	boost::mt19937 rng;
	if(argc>2){
		rng.seed(atoi(argv[2]));
	}
	bool print=false;
	if(argc>3){
		print = true;
	}

	typedef sg_dvv G;
	sg_dvv g(size);
	g.reshape(0);
	boost::generate_random_graph(g, size, ne, rng);
	size_t e=boost::num_edges(g);
	size_t n=boost::num_vertices(g);
	std::cout << "generated " << e << " edges, " << n << " vertices\n";
	treedec::check(g);
	g.make_symmetric(false);
	treedec::check(g);
	e = boost::num_edges(g);

	// boost::add_edge(0,1,g);
	e=boost::num_edges(g);

	std::cout << "symmetric " << e << " edges\n";


	unsigned i=0;
	auto E=boost::edges(g);
	for(;E.first!=E.second; ++E.first){
		++i;

		std::cout << boost::source(*E.first, g) << " -- " <<
			boost::target(*E.first, g) << "\n";

		if(i==5) break;

	}

	std::deque<unsigned long > iso;
	auto V=boost::vertices(g);
	for(;V.first!=V.second; ++V.first){
			if(boost::out_degree(*V.first, g)){

			}else{
				iso.push_back(*V.first);
			}
	}

#if 0 // incomplete
	(vec_to_treedec now works anyway!)
	for(auto i: iso){
		boost::remove_vertex(i, g);
	}
#endif

	n = boost::num_vertices(g);

	V=boost::vertices(g);
	for(;V.first!=V.second; ++V.first){
		assert(*V.first<n);
		// assert(boost::out_degree(*V.first,g));
	}

	std::cout << "removed " << iso.size() <<"\n";

	i = 0;
	// boost md does not like cliques.
	trace3("clique check", n, e, size);
	if((n*(n-1u)) == boost::num_edges(g)){ untested();
		exit(0);
	}else{
		itested();
	}

	std::vector<int> inverse_perm(n, 0);
	std::vector<int> supernode_sizes(n, 1);
	auto id = boost::get(boost::vertex_index, g);
	std::vector<int> degree(n, 0);
	std::vector<int> io(n, 0);
	std::vector<int> o(n, 0);

	G h(g);
	assert(boost::num_edges(g)==boost::num_edges(g));

	// assert(n + iso.size() == size);

	/*
	 * (Graph& g,
	 *  DegreeMap degree,
	 *  InversePermutationMap inverse_perm,
	 *  PermutationMap perm,
	 *  SuperNodeMap supernode_size,
	 *  int delta,
	 *  VertexIndexMap vertex_index_map)
	 */

	unsigned ub=-1;

	unsigned w =
#ifndef HAVE_MINDEGREE_FORK
		0;
	itested();
#endif
	boost::minimum_degree_ordering
		(g,
		 boost::make_iterator_property_map(&degree[0], id, degree[0]),
		 &io[0],
		 &o[0],
		 boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
		 0,
		 id
#ifdef HAVE_MINDEGREE_FORK
		 , ub
#endif
		);
	typename treedec::graph_traits<G>::treedec_type t;
	g=h; // restore
	treedec::check(g);
	assert(boost::num_edges(g)==boost::num_edges(h));

	std::vector<bool> permcheck(n, false);
	for(auto x : o){
		permcheck[x]=true;
	}
	for(auto x : permcheck){
		assert(x);
	}

	int status = 0;
#ifdef CHECK_INEFFICIENT_VARIANT
	treedec::ordering_to_treedec(g, o, t ); // can kill g!
   h=g; // restore

	status = treedec::check_treedec(g,t);
	std::cout << "bagsize " << treedec::get_bagsize(t) << " status " << status <<"\n";
	std::cout << "same as bagsize! " << w <<"\n";
	assert(w == treedec::get_bagsize(t)); /// checks if BMD works!
#endif

	std::vector<unsigned>& ou=reinterpret_cast<std::vector<unsigned>&>(o);
	treedec::draft::vec_ordering_to_tree(g, ou, t );
	assert(boost::num_edges(t)+1==boost::num_vertices(t));
	assert(boost::num_edges(t)+1==boost::num_vertices(g));

	if(print){
		auto P = treedec::grtdprinter<G>(std::cout, g, "test");
		treedec::draft::vec_ordering_to_tree(g, ou, P );
	}


	// does not work on vectors...
	// status=treedec::check_treedec(g,t);

	if (!status) std::cout << "treedec is valid!!\n";
	std::cout << "bagsize " << treedec::get_bagsize(t) << " status " << status <<"\n";
	std::cout << "boost_minDegree_ordering said " << w << "\n";

	incomplete();
//	assert(!status);

	if(w != treedec::get_bagsize(t)){
	}
	assert(w == treedec::get_bagsize(t));

#endif // HAVE_GALA

}
