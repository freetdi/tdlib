#ifndef THTHREAD_H
#define THTHREAD_H

#include <treedec/thorup.hpp>

template<class G, template<class H, class ... > class cfgt=treedec::algo::default_config>
class TH_THREAD : public TWTHREAD<G, cfgt> {
public:
    typedef cfgt<G> CFG;
    typedef treedec::thorup<G, cfgt> algo_type;
    typedef TWTHREAD<G, cfgt> base;
	 using base::_g;
	 using base::_result;
    TH_THREAD( G const&g, const std::string& name )
        : base(g, name, 0) , _TH(NULL)
    {
		 //boost::print_graph(_work);
        _TH = new algo_type(_g);
        base::go();
    }

    void do_print_results(std::ostream& o)
    {
#if 1
		 // pkks ordering_to_treedec
		 // also seems to work on oriented directed graphs.
		  treedec::grtdprinter<G> P(o, _g); // index map?
		  size_t numbags = boost::num_vertices(_g) + 1; // why +1??
		  P.head(numbags, _result);
		  assert(_TH);
		  _TH->get_tree_decomposition(P);
#else
		  // use tdlib to compute decomposition
		  // this one needs a symmetric graph...
		  auto elimord=_TH->get_elimord();
        base::print_results_order(o, elimord);
#endif
    }

    void run() { untested();
		 _TH->do_it();

		 unsigned r = _TH->get_bagsize();
		 std::cerr << "thorup RESULT " << r << "\n";
		 base::commit_result(r);
		 base::unlock_results();
    }
private:
#ifdef HAVE_GALA_GRAPH_H
	algo_type* _TH;
	//(_work);
#else
   incomplete
   G _work;
#endif
   decomp_t<G> _t;
};

#endif
