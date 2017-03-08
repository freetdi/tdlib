#ifndef GENERIC_ELIMINATION_SEARCH_OVERLAY
#define GENERIC_ELIMINATION_SEARCH_OVERLAY

#include <boost/graph/adjacency_list.hpp>

namespace treedec{

namespace gen_search{

template <typename UnderlyingG_t, typename OverlayG_t>
class overlay{
public:
    overlay(UnderlyingG_t &G_input)
      : G(G_input)
    {}

    UnderlyingG_t& underlying(){
        return &G;
    }

private:
    const UnderlyingG_t &G;
    OverlayG_t O;
};

} //namespace gen_search

} //namespace treedec

#endif //guard
