//TODO: header
//TODO: this is algo_traits.hpp

#ifndef TREEDEC_CONFIG_TRAITS_HPP
#define TREEDEC_CONFIG_TRAITS_HPP

// get types from config objects

namespace treedec{

namespace config{

namespace get{
/*--------------------------------------------------------------------------*/
template<typename T>
struct tovoid { typedef void type; };
/*--------------------------------------------------------------------------*/
#define DC(a, b, c) a::b ## c
#define DECLARE_CONFIG_GETTER(what) \
template<class CFG, class fallback, class X=void> \
struct what{ \
	typedef fallback type; \
}; \
template<class CFG, class fallback> \
struct what<CFG, fallback, typename tovoid<typename DC(CFG, what, _type) >::type >{ \
	typedef typename DC(CFG, what, _type) type; \
}
/*--------------------------------------------------------------------------*/
DECLARE_CONFIG_GETTER(degree);
DECLARE_CONFIG_GETTER(degs);
DECLARE_CONFIG_GETTER(kernel);
DECLARE_CONFIG_GETTER(edge_overlay_graph); // yuck
/*--------------------------------------------------------------------------*/
#undef DC
#undef DECLARE_CONFIG_GETTER
}

}

}
#endif
