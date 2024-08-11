#include <treedec/algo.hpp>
#include <treedec/tuple_td.hpp>
#include <gala/boost_copy.h>

#include "../src/tr_iodec.h"
#include "../src/tr_bag.h"

// TODO: is this needed?
namespace treedec
{ //

    namespace detail{
        template <class G,
                  template <class G_, class...> class C>
        using tr_ = TR<G>;
    }

    // to combinations?
    template <typename G, typename T_t>
    void exact_decomposition_tr(G &g, T_t &T) { untested();
        TR<G> a(g);
        treedec::grtdprinter<G> P(std::cerr, g);
        // auto alg = draft::exact_decomposition<G, algo::default_config, detail::exact_ta_>(g);
         a.do_it(P);
    }
}

#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
// #include <gala/examples/ssg16i.h>

#if 1 // slow preprocessing
typedef ssg_16i TR__GRAPH;
#else // broken preprocessing
typedef sg_dvv16 TR__GRAPH;
#endif

template <class G, template <class H, class...> class cfgt = treedec::algo::default_config>
class TR_THREAD : public TWTHREAD<G, cfgt> { //
public:
    typedef TWTHREAD<G, cfgt> base;
    typedef TR__GRAPH G_work;

#if 1 // tr_myset workaround
    typedef decomp_t<TR__GRAPH> T; // BUG.
#else
    typedef decomp_t<G_work> T;
#endif

private:
    TR_THREAD() { untested();
        incomplete();
    }
    TR_THREAD(const TR_THREAD &) { untested();
        incomplete();
    }
    TR_THREAD(const TR_THREAD &&) {
        incomplete();
    }

public:
    TR_THREAD(G &g, const std::string &name = "TR")
        : base(g, name, 0), _testg(g) {
        treedec::check(g);
        base::go();
    }

    void do_print_results(std::ostream &o) {
        trace2("TRR", base::_result, treedec::get_bagsize(_t));

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
        TR<G> a(_work);
        a.do_it(_t);
        //  a.store(_t);

        unsigned r = treedec::get_bagsize(_t); // inefficient
        trace1("tr done", r);
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
