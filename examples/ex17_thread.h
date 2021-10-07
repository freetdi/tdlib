//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> bald_t;
//
//
#include <treedec/algo.hpp>
#include <treedec/exact_ta.hpp>
#include <treedec/tuple_td.hpp>
#include <gala/boost_copy.h>

//typedef cbset::BSET_DYNAMIC<K> BSET;

class ta_kernel_config{
    ///...
};

namespace ta{
  template<class GraphType>
  struct cfg_8{
    static constexpr unsigned max_vertex_index=255u;
  };
  template<class GraphType>
  struct cfg_7{
    static constexpr unsigned max_vertex_index=127u;
  };
  template<class GraphType>
  struct cfg_6{
    static constexpr unsigned max_vertex_index=63u;
  };

}

// TODO: move
namespace treedec{ //

    namespace detail{
        template<class G,
            template<class G_, class ...> class C>
        using exact_ta_=exact_ta<G, C>;
    }

    // to combinations?
template <typename G, typename T_t>
void exact_decomposition_ta(G &g, T_t &T, unsigned lb_bs=0)
{ untested();
    auto alg=draft::exact_decomposition<G, algo::default_config, detail::exact_ta_>(g);
    return alg.do_it(T, lb_bs);
}

}

template<class G, class O>
void print_matrix(G const& g, O& o)
{
	unsigned n=boost::num_vertices(g);

	for(unsigned i=0; i<n; ++i){ itested();
//		auto a=boost::adjacent_vertices(g, i);
//		for(;a.first!=a.second; ++a.first){
//			std::cout << "\n";
//
//		}
		for(unsigned j=0; j<n; ++j){ itested();
			if(boost::edge(i,j,g).second){ itested();
				o << "1";
			}else{ itested();
				o << "0";
			}
		}
		
		o << "\n";
	}
}

template<class G, class O, class S>
void print_matrix(G const& g, O const& o, S& s)
{ itested();
    typedef typename boost::graph_traits<G>::vertex_descriptor VD;
	unsigned n=boost::num_vertices(g);

	for(unsigned i=0; i<n; ++i){ itested();
//		auto a=boost::adjacent_vertices(g, i);
//		for(;a.first!=a.second; ++a.first){ itested();
//			std::cout << "\n";
//
//		}
		s << "c ";
		for(unsigned j=0; j<n; ++j){ itested();
                    VD I=boost::vertex(o[i], g);
                    VD J=boost::vertex(o[j], g);
			if(boost::edge(I, J, g).second){ itested();
				s << "1";
			}else{ itested();
				s << "0";
			}
		}
		
		s << "\n";
	}
}


#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#if 1 // slow preprocessing
    typedef ssg_16i EX17__GRAPH;
#else // broken preprocessing
typedef sg_dvv16 EX17__GRAPH;
#endif

template<class G, template<class H, class ... > class cfgt=treedec::algo::default_config>
class EX17_THREAD : public TWTHREAD<G, cfgt> { //
public:
    typedef TWTHREAD<G, cfgt> base;

    typedef EX17__GRAPH G_work; // BUG
    typedef decomp_t<G_work> T; // BUG.
private:
    EX17_THREAD(){ untested();
        incomplete();
    }
    EX17_THREAD(const EX17_THREAD&){ untested();
        incomplete();
    }
    EX17_THREAD(const EX17_THREAD&&){
        incomplete();
    }
public:

    EX17_THREAD(G& g, const std::string& name="EX17")
        : base(g, name, 0), _testg(g)
    {
        treedec::check(g);
        /// std::cerr << "orig\n";
        // rcmk(gg);
        base::go();
    }

    void do_print_results(std::ostream& o)
    {
        trace2("EXR", base::_result, treedec::get_bagsize(_t));
        base::print_results_tree(o, _t);
    }

    void run() {
# if 1 // HACK
        // PP old(?) interface (used in exact.hpp) refuses gala
        G_work _work;
        boost::copy_graph(base::_g, _work);
//        treedec::preprocessing(base::_g, bags, low);
#else
        typedef G G_work;
#endif


        auto alg=treedec::draft::exact_decomposition<G_work,
             treedec::algo::default_config, // FIXME
             treedec::detail::exact_ta_>(_work);
        //     treedec::detail::exact_ta_>(base::_g);

        // set_bagsize(lb,-1); do_it()?!
        // get_treedec() ... 
        alg.try_it(_t, 0);

        unsigned r = treedec::get_bagsize(_t); // inefficient
        trace1("excut done", r);
        base::commit_result(r);
        base::unlock_results();
//        kill(getpid(), SIGTERM);
    }
private:
//   result_t _result;
    T _t;
    G const& _testg;
};


// vim:ts=8:sw=4:et
