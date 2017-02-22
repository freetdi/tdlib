
#include <boost/graph/copy.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <stdlib.h>
#include <tdlib/elimination_orderings.hpp>
#include <tdlib/graph.hpp>
#include <tdlib/trace.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> balu_t;

// typedef bald_t parsegr_t;

int main(int argc, char** argv)
{
	size_t size=1<<5;
	if(argc>1) size=1<<atoi(argv[1]);
	size_t ne=size*size/8 + 1;

	boost::mt19937 rng;
	if(argc>2){
		rng.seed(atoi(argv[2]));
	}

	typedef bald_t G;
	bald_t g;
	boost::generate_random_graph(g, size, ne, rng);
	size_t e=boost::num_edges(g);
	size_t n=boost::num_vertices(g);
	std::cout << "generated " << e << " edges, " << n << " vertices\n";

	BOOST_AUTO(EE, boost::edges(g));
	for(;EE.first!=EE.second; ++EE.first){
		BOOST_AUTO(s, boost::source(*EE.first, g));
		BOOST_AUTO(t, boost::target(*EE.first, g));
		if(!boost::edge(t,s,g).second){
			boost::add_edge(t,s,g);
		}
	}

	// boost::add_edge(0,1,g);
	e=boost::num_edges(g);

	std::cout << "symmetric " << e << " edges\n";

	unsigned i=0;
	BOOST_AUTO(E, boost::edges(g));
	for(;E.first!=E.second; ++E.first){
		++i;

		std::cout << boost::source(*E.first, g) << " -- " <<
			boost::target(*E.first, g) << "\n";

		if(i==5) break;

	}

	std::deque<unsigned long > iso;
	BOOST_AUTO(V, boost::vertices(g));
	for(;V.first!=V.second; ++V.first){
			if(boost::out_degree(*V.first, g)){

			}else{
				iso.push_back(*V.first);
			}
	}

#ifdef REMOVE_ISOLATED
	// no longer necessary.
	for(auto i: iso){
		boost::remove_vertex(i, g);
	}
#endif

	n = boost::num_vertices(g);

	V=boost::vertices(g);
	for(;V.first!=V.second; ++V.first){
		assert(*V.first<n);
//		assert(boost::out_degree(*V.first,g));
	}

	std::cout << "tagged " << iso.size() <<"\n";

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
	BOOST_AUTO(id, boost::get(boost::vertex_index, g));
	std::vector<int> degree(n, 0);
	std::vector<int> io(n, 0);
	std::vector<int> o(n, 0);

	G h(g);
	assert(boost::num_edges(g)==boost::num_edges(g));

//	assert(n + iso.size() == size);

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

	int w =
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

	int status;
#ifdef TRY_INEFFICIENT_VARIANT
	treedec::ordering_to_treedec(g, o, t ); // can kill g!
   h=g; // restore

	status = treedec::check_treedec(g,t);
	std::cout << "bagsize " << treedec::get_bagsize(t) << " status " << status <<"\n";
	std::cout << "same as bagsize! " << w <<"\n";
	assert(w == treedec::get_bagsize(t)); /// checks if BMD works!
#endif

	treedec::draft::vec_ordering_to_tree(g, o, t );


	status=treedec::check_treedec(g,t);
	if (!status) std::cout << "treedec is valid!!\n";
	std::cout << "bagsize " << treedec::get_bagsize(t) << " status " << status <<"\n";
	std::cout << "boost_minDegree_ordering said " << w << "\n";
	assert(!status);

	if(w != treedec::get_bagsize(t)){
	}
	assert(w == treedec::get_bagsize(t));


}
