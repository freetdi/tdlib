#include <set>
#include <vector>

typedef int var_t;
typedef std::set<unsigned int> sdcc_bagtype;
struct tree_dec_node {
//  typedef boost::vertex_property_tag kind; // do we need it?
	sdcc_bagtype bag;
  std::set<var_t> alive;
  unsigned weight;
};
struct tree_dec_node2 {
//  typedef boost::vertex_property_tag kind; // do we need it?
	std::vector<unsigned> bag;
  std::set<var_t> alive;
  unsigned weight;
};

#include <treedec/treedec_traits.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, tree_dec_node> sbib_tdt;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              boost::bidirectionalS, tree_dec_node2> vbib_tdt;

TREEDEC_TREEDEC_BAG_TRAITS(sbib_tdt, bag);
TREEDEC_TREEDEC_BAG_TRAITS(vbib_tdt, bag);

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>

#include <iostream>
#include <vector>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <treedec/preprocessing.hpp>
#include <treedec/graph.hpp>
#include <treedec/thorup.hpp>
#ifdef HAVE_GALA_GRAPH_H
#include <treedec/exact_ta.hpp>
#endif


typedef boost::adjacency_list<boost::vecS, boost::vecS,
		  boost::undirectedS, treedec::bag_t> bag_t_tdt;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		  boost::property<treedec::bag_t, std::vector<unsigned> > > prop_tdt;
// this one is used in thorup for example.
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
		  boost::property<treedec::bag_t, std::set<unsigned int> > > example_tdt;


BOOST_STATIC_ASSERT(
				std::is_same< typename boost::vertex_bundle_type<sbib_tdt>::type,
								tree_dec_node >::value );

#include <boost/graph/copy.hpp>

#if 0
namespace boost{

  template<class U>
  inline void
  put(const put_get_helper<std::vector<U>,
		 property_map<prop_tdt , vertex_all_t >::type  >& pa, unsigned long k,
		 property<treedec::bag_t, std::vector<unsigned> > const& v)
  { untested();
	  auto& b=static_cast<bagstuff::treebagpmap<sbib_tdt> &>(pa)[k];
	  b.clear();
	  for(auto const& i : v.m_value){ untested();
		  treedec::push(b, i);
	  }
  }
  template<class U>
  inline void
  put(const put_get_helper<std::vector<U>,
		 property_map<prop_tdt , vertex_all_t >::type  >& pa, unsigned long k,
		 std::set<unsigned>const& v)
  { untested();
	  auto& b = pa[k];
	  b.clear();
	  for(auto const& i : v){ untested();
		  b.insert(i);
	  }
  }
  template<class U>
  inline void
  put(const put_get_helper<std::vector<U>,
		 property_map<prop_tdt , vertex_all_t >::type  >& pa, unsigned long k,
		 std::vector<unsigned>const& v)
  { untested();
	  auto& b = pa[k];
	  b.clear();
	  for(auto const& i : v){ untested();
		  b.insert(i);
	  }
  }
} // boost
#endif

// too late?
#include <treedec/thorup.hpp>

int main(int, char**)
{
	sbib_tdt st; // such as used in sdcc
	vbib_tdt vt;
	bag_t_tdt bt; // bundled property, using bag_t
	prop_tdt pt; // bag_t is a property, could be anything.
	example_tdt et;
	auto sv=boost::add_vertex(st);
	boost::add_vertex(pt);

	auto m=boost::get(treedec::bag_t(), st);
	auto xx=boost::get(m, 0);
	boost::get(treedec::bag_t(), st, unsigned(0));
	

	auto sb=boost::get(&tree_dec_node::bag, st, sv);
	auto sb2=boost::get(treedec::bag_t(), st, sv);

	sb.insert(0);

	boost::copy_graph(st, bt);
	boost::copy_graph(bt, st);
	boost::copy_graph(pt, bt);
	boost::copy_graph(et, bt);
	boost::copy_graph(et, st);
//	boost::copy_graph(bt, pt); not yet.
//	boost::copy_graph(st, pt); not yet.
//
	auto n=*boost::vertices(pt).first;
	auto mp=get(treedec::bag_t(), pt);
	auto& b=get(treedec::bag_t(), pt, n);
	treedec::push(b, 1);
	treedec::push(b, 2);
	treedec::push(b, 3);

	n=boost::add_vertex(pt);
	{
		auto& b=get(treedec::bag_t(), pt, n);
		treedec::push(b, 4);
		treedec::push(b, 5);
		treedec::push(b, 6);
	}

	boost::print_graph(pt);
	while(boost::num_vertices(st)){
		boost::remove_vertex(0, st);
	}
	std::cout<<"copy\n";
	boost::copy_graph(pt, st);
	boost::copy_graph(pt, vt);

	auto N=boost::num_vertices(st);
	assert(N==2);
	auto p=boost::vertices(st);

	for(;p.first!=p.second; ++p.first){
		std::cout << *p.first << ":";
		auto q=boost::get(&tree_dec_node::bag, st, *p.first);
		for(auto i : q){
			std::cout << " " << i;
		}
		std::cout <<"\n";
	}

	{
		typedef boost::adjacency_list<boost::vecS, boost::vecS,
				  boost::bidirectionalS > G;
		typedef treedec::thorup<G> alg_B;
		G h;
		alg_B B(h);
		B.do_it();
		auto X=B.get_tree_decomposition();
		boost::copy_graph(X, st);

		B.get_tree_decomposition(st);
	}
	auto const& cpt(pt);
	boost::copy_graph(cpt, st);
	boost::copy_graph(cpt, vt);
	auto const& cbt(bt);
	boost::copy_graph(cbt, st);
//	boost::copy_graph(cbt, vt); not yet
	boost::copy_graph(cbt, st);

}
