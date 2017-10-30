#ifndef EXAMPLES_GALA_GRAPHS_H
#define EXAMPLES_GALA_GRAPHS_H

typedef simplegraph_vector_bs svbs;

template<class T, class...>
using tdDEGS = misc::DEGS<T>;


// possibly oriented directed simple loopless graph
template<class G>
struct odsvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=true;
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    static constexpr bool force_symmetric=false;
    static constexpr bool force_oriented=true; // not used, how?!
    typedef tdDEGS<G> degs_type; // obsolete.
};

// parsed raw data.
template<class G>
struct dpvv_config : uvv_config<G> {
    static constexpr bool force_simple=false;
    static constexpr bool force_symmetric=false;
    static constexpr bool is_directed=true;
};

typedef gala::graph<std::vector, std::vector, uint16_t, dpvv_config> sg_dpvv16;
typedef gala::graph<std::vector, std::vector, uint32_t, dpvv_config> sg_dpvv32;

typedef gala::graph<std::vector, std::vector, uint16_t, odsvv_config> sg_odsvv16;


//typedef gala::graph<std::set, std::vector, uint16_t> ssg16;
#endif
