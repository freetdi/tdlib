
#warning "typocatching transitional header"
#include "separator_algorithm.hpp"

namespace treedec{
	// don't use this.
template <typename G_t, typename T_t>
void seperator_algorithm(G_t &G, T_t &T){
    return separator_algorithm(G, T);
}
}
