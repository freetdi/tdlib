#ifndef FITHREAD_H
#define FITHREAD_H

// #include <boost/graph/graph_utility.hpp>

#define GWORKFI G // ?

// TODO: merge into SEVERAL_FI
template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class FI_THREAD : public TWTHREAD<G, cfgt> {
public:
    typedef cfgt<G> CFG;
    typedef treedec::impl::fillIn<GWORKFI, treedec::algo::default_config> algo_type;
    typedef TWTHREAD<GWORKFI, cfgt> base;
	 using base::_g;
	 using base::_result;
    FI_THREAD( G const&g, const std::string& name )
        : base(g, name, 0), _work(g) // <= stored here
			 , _FI(NULL)
    {
        _FI = new algo_type(_work);
        treedec::check(g);
        treedec::check(_work);
		trace2("FI_THREAD3", boost::num_vertices(_g), boost::num_edges(_g));
// #ifdef HAVE_GALA_GRAPH_H
//         h = g;
// #endif
        base::go();
    }

    void do_print_results(std::ostream& o)
    {
        std::cerr<< "c size " << boost::num_vertices(_work) << "\n";
        // auto &g=TWTHREAD<G>::_g;
		  treedec::grtdprinter<G> P(o, _work);
		  size_t numbags = boost::num_vertices(_work); // ask P?!
		  P.head(numbags, _result);
		  assert(_FI);
		  _FI->get_tree_decomposition(P);
    }

    void run() {
#if 1
		_FI->do_it();

#else
        typedef typename boost::graph_traits<sg_dvv16>::vertex_descriptor vertex_descriptor;
        std::vector<vertex_descriptor> O;
        treedec::impl::endless_fillIn_ordering(_work, O);
#endif
        unsigned r = _FI->get_bagsize();
        assert(boost::num_vertices(_work) || r==0);
        base::commit_result(r);
        base::unlock_results();

//    NOT YET
//     G = H;
//     treedec::MSVS(G, T);
//     w = treedec::get_bagsize(T);
//     update_best_bagsize(w, "FI+MSVS");
//
    }
private:
#ifdef HAVE_GALA_GRAPH_H
   // work on svbs
   // INCOMPLETE, hardwire 16bit!
   GWORKFI _work;
	algo_type* _FI;
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

// FIXME: this is specific to FI
template<class G,
             template<class GG, template<class G_, class ...> class CFGT> class A,
             template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class SEVERAL_FI_THREAD : public TWTHREAD<G, cfgt> {
public:
    typedef cfgt<G> CFG;
    typedef treedec::impl::fillIn<GWORKFI, cfgt> algo_type;
    typedef TWTHREAD<GWORKFI, cfgt> base;
	 using base::_g;
	 using base::_result;
    SEVERAL_FI_THREAD( G const&g, const std::string& name )
        : base(g, name, 0), _work(g)
    { itested();
        base::go();
    }

    void do_print_results(std::ostream& o)
    {
        std::cerr<< "c size " << boost::num_vertices(_work) << "\n";
        // auto &g=TWTHREAD<G>::_g;
		  treedec::grtdprinter<G> P(o, _work);
		  size_t numbags = boost::num_vertices(_t);
		  P.head(numbags, _result);
		  boost::copy_graph(_t, P);
    }

    void run() {
		 A <GWORKFI, cfgt> a(_work);
		 a.do_it();
		 a.get_tree_decomposition(_t);


        unsigned r = get_bagsize(_t);
        assert(boost::num_vertices(_work) || r==0);
        base::commit_result(r);
        base::unlock_results();

//    NOT YET
//     G = H;
//     treedec::MSVS(G, T);
//     w = treedec::get_bagsize(T);
//     update_best_bagsize(w, "FI+MSVS");
//
    }
private:
   GWORKFI _work;
   decomp_t<GWORKFI> _t;
}; // SEVERAL_FI_THREAD

template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class PPFI_THREAD
    : public SEVERAL_FI_THREAD<G, treedec::pending::PP_FI, cfgt>
{
public:
	typedef SEVERAL_FI_THREAD<G, treedec::pending::PP_FI, cfgt> base;
	template<class A, class B>
	explicit PPFI_THREAD(A const& a, B const& b) : base(a,b){}
};

template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class PPFITM_THREAD
    : public SEVERAL_FI_THREAD<G, treedec::pending::PP_FI_TM, cfgt>
{
public:
	typedef SEVERAL_FI_THREAD<G, treedec::pending::PP_FI_TM, cfgt> base;
	template<class A, class B>
	explicit PPFITM_THREAD(A const& a, B const& b) : base(a,b){}
};

template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class PPMD_THREAD
    : public SEVERAL_FI_THREAD<G, treedec::comb::PP_MD, cfgt>
{
public:
	typedef SEVERAL_FI_THREAD<G, treedec::comb::PP_MD, cfgt> base;
	template<class A, class B>
	explicit PPMD_THREAD(A const& a, B const& b) : base(a,b){}
};

template<class G, template<class H, class ... >
                  class cfgt=treedec::algo::default_config>
class FITM_THREAD
    : public SEVERAL_FI_THREAD<G, treedec::comb::FI_TM, cfgt>
{
public:
	typedef SEVERAL_FI_THREAD<G, treedec::comb::FI_TM, cfgt> base;
	template<class A, class B>
	explicit FITM_THREAD(A const& a, B const& b) : base(a,b){}
};

#endif
