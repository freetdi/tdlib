//TODO: header

#ifndef TD_COPY_HPP
#define TD_COPY_HPP

#include <boost/graph/copy.hpp>

#ifdef USE_GALA
#include <gala/boost_copy.h>
#endif

namespace treedec{

using boost::num_vertices;

// copy s to t, but weed out multiedges
// rationale: when boost::copying a directed graph (0,1),(1,0) to an empty
// undirected multigraph, then we get two edges, while one might be more useful
// and more efficient sometimes.
template<class S, class T /* class M */>
void copy_trace(const S& s, T& t /*, M map=identity */)
{
    // TODO template metamagic...
    if(boost::num_vertices(t)){ untested();
        assert(false);
    }else if(!boost::is_multigraph<T>()){
        untested();
        boost::copy_graph(s, t);
    }else if(boost::is_directed(t)){
        boost::copy_graph(s, t);
    }else{
        t = MOVE(T(boost::num_vertices(s)));
        auto b=boost::edges(s);
        for(; b.first!=b.second; ++b.first){
            typename boost::graph_traits<S>::edge_descriptor f = *b.first;
            auto u=boost::source(f, s);
            auto v=boost::target(*b.first, s);
            if(u<v){
                boost::add_edge(u, v, t);
            }
        }
    }
}

} //namespace treedec

#endif

// vim:ts=8:sw=4:et:
