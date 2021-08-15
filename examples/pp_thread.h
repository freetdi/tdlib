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
    {
        _PP = new algo_type(_work);
        treedec::check(g);
        treedec::check(_work);
		trace2("PP_THREAD3", boost::num_vertices(_g), boost::num_edges(_g));
// #ifdef HAVE_GALA_GRAPH_H
//         h = g;
// #endif
        base::go();
    }

    void do_print_results(std::ostream& o) {
//        _PP->get_tree_decomposition(_printer); ?
        trace2("PP", base::_result, treedec::get_bagsize(_t));
        base::print_results_tree(o, _t, &_work);
    }

    void run() {
		_PP->do_it();

		  CFG::message(bLOG, "PP edges left %d\n", _PP->num_edges());
        unsigned r = _PP->get_bagsize();
        // assert(boost::num_vertices(_work) || r==0);
        _PP->get_tree_decomposition(_t);
        base::commit_result(r);
        base::unlock_results();
    }
	 ~PP_THREAD(){
		 delete _PP;
	 }
private:
#ifdef HAVE_GALA_GRAPH_H
   // work on svbs
   // INCOMPLETE, hardwire 16bit!
   GWORKFI _work;
	//(_work);
   decomp_t<GWORKFI> _t;
#else
   decomp_t<G> _t;
   G _work;
#endif
	algo_type* _PP;
   // why?
   //iorder_t _elimord;
};


#endif
