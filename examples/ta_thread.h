#include <treedec/algo.hpp>
#include <treedec/tuple_td.hpp>
#include <gala/boost_copy.h>
#include <treedec/exact_ta.hpp>

// TODO: is this needed?
namespace treedec {
    namespace detail {
        template <class G>
        using ta_ = treedec::exact_ta<G>;
    }

    // to combinations?
    template <typename G, typename T_t>
    void exact_decomposition__ta(G &g, T_t &T, unsigned lb_bs=2) { untested();
    auto alg = treedec::exact_ta<G>(g);
    return alg.do_it(T, lb_bs);
    }
}

#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
// #include <gala/examples/ssg16i.h>

#if 1 // slow preprocessing
typedef ssg_16i TA__GRAPH;
#else // broken preprocessing
typedef sg_dvv16 TR__GRAPH;
#endif

template <class G, template <class H, class...> class cfgt = treedec::algo::default_config>
class TA_THREAD : public TWTHREAD<G, cfgt> { //
public:
    typedef TWTHREAD<G, cfgt> base;
    typedef TA__GRAPH G_work;

#if 1 // tr_myset workaround
    typedef decomp_t<TA__GRAPH> T; // BUG.
#else
    typedef decomp_t<G_work> T;
#endif

private:
    TA_THREAD() { untested();
        incomplete();
    }
    TA_THREAD(const TA_THREAD &) { untested();
        incomplete();
    }
    TA_THREAD(const TA_THREAD &&) {
        incomplete();
    }

public:
    TA_THREAD(G &g, const std::string &name = "TA")
        : base(g, name, 0), _testg(g) {
        treedec::check(g);
        base::go();
    }

    void do_print_results(std::ostream &o) {
        trace2("TAR", base::_result, treedec::get_bagsize(_t));

        base::print_results_tree(o, _t);
    }

    void run() {
#if 1 // HACK
      // PP old(?) interface (used in exact.hpp) refuses gala
        G_work _work;
        boost::copy_graph(base::_g, _work);
//        treedec::preprocessing(base::_g, bags, low);
#else
        typedef G G_work;
#endif

        // auto alg = treedec::draft::exact_decomposition<G_work,
        //                                                treedec::algo::default_config, // FIXME
        //                                                treedec::detail::exact_ta_>(_work);
        // //     treedec::detail::exact_ta_>(base::_g);

        // // set_bagsize(lb,-1); do_it()?!
        // // get_treedec() ...
	//
	  auto alg=treedec::exact_ta<G_work,cfgt>(_work);
        unsigned lb_bs = 0;
	  //  a.store(_t);
	 alg.do_it(_t, lb_bs);
        unsigned r = treedec::get_bagsize(_t); // inefficient
        trace1("ta done", r);
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
