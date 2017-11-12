
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <treedec/treedec_traits.hpp>
#include <treedec/graph_traits.hpp>

#ifdef HAVE_GALA_GRAPH_H
#include <gala/graph.h>
#endif


#ifdef HAVE_GALA_GRAPH_H
#include <gala/boost.h>
#include <treedec/printer.hpp>
#include <boost/graph/copy.hpp>
#include <gala/boost_copy.h>
// undirected simple loopless graph
template<class G>
struct uvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=false;
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    // static constexpr bool force_symmetric=true; // meaninngless (undirected)
    // typedef tdDEGS<G> degs_type; // obsolete.
};
typedef gala::graph<std::vector, std::vector, uint16_t, uvv_config> sg_dvv16;
typedef gala::graph<std::vector, std::vector, uint32_t, uvv_config> sg_dvv32;

struct test{
    test(){
        sg_dvv16 a;
   treedec::draft::printer<sg_dvv16> x(std::cout, a);
   boost::add_vertex(x);
    }
} x;
#endif
#include <treedec/treedec.hpp>
#include <treedec/printer.hpp>
#include <boost/graph/copy.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> ALVVD;
template<class G>
using decomp_t = typename treedec::graph_traits<G>::treedec_type;

int main(){
	return 0;

	treedec::graph_traits<ALVVD>::treedec_type t;

	ALVVD g;
	treedec::grtdprinter<ALVVD> P(std::cerr, g);

	boost::copy_graph(t, P);
#ifdef HAVE_GALA_GRAPH_H
	 typedef gala::graph<std::vector, std::vector, unsigned int> GG;
	 gala::graph<std::vector, std::vector, unsigned int> gg;
	 decomp_t<GG> TT;
	boost::copy_graph(TT, P);
#endif

}
