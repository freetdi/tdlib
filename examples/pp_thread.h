#ifndef PP_THREAD_H
#define PP_THREAD_H

#define GWORKFI G // ?

template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class PP_THREAD : public TWTHREAD<G, cfgt> {
public:
    typedef cfgt<G> CFG;
    typedef treedec::impl::preprocessing<G> algo_type;
    typedef TWTHREAD<GWORKFI, cfgt> base;
	 using base::_g;
	 using base::_result;
    PP_THREAD( G const&g, const std::string& name )
        : base(g, name, 0), _work(g) // <= stored here
			 , _PP(NULL)
    { untested();
        _PP = new algo_type(_work);
        treedec::check(g);
        treedec::check(_work);
		trace2("PP_THREAD3", boost::num_vertices(_g), boost::num_edges(_g));
// #ifdef USE_GALA
//         h = g;
// #endif
        base::go();
    }

    void do_print_results(std::ostream& o)
    { untested();
        // auto &g=TWTHREAD<G>::_g;
		  // treedec::grtdprinter<G> P(o, _work);
		  // size_t numbags = boost::num_vertices(_work); // ask P?!
		  // P.head(numbags, _result);
		  // assert(_PP);
		  // _PP->get_tree_decomposition(P);
    }

    void run() { untested();
		_PP->do_it();

		  CFG::message(bLOG, "PP done, edges left: %d", _PP->num_edges());
        unsigned r = _PP->get_bagsize();
        // assert(boost::num_vertices(_work) || r==0);
		  incomplete();
        base::commit_result(r);
        base::unlock_results();
    }
	 ~PP_THREAD(){
		 delete _PP;
	 }
private:
#ifdef USE_GALA
   // work on svbs
   // INCOMPLETE, hardwire 16bit!
   GWORKFI _work;
	algo_type* _PP;
	//(_work);
   decomp_t<GWORKFI> _t;
#else
   incomplete
   decomp_t<G> _t;
   G _work;
#endif
   // why?
   //iorder_t _elimord;
};


#endif
