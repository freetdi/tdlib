// Felix Salfelder, 2017
//
// (c) 2014-2016 Goethe-Universität Frankfurt
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//
//
// almost-tree-decompositions
// a skeleton is a list of vertices with vertices attached to it.

#ifndef TD_SKELETON_HPP
#define TD_SKELETON_HPP

#include <boost/graph/graph_traits.hpp>

namespace treedec{ //

namespace draft{

template<class G, class N, class O>
class SKELETON{
public: // types
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
    SKELETON(G const& g, N const& n, O const& o)
      : _g(g), _n(n), _o(o)
    { untested();
    }
public:
    size_t size() const{ untested();
        return _o.size();
    }
//    bag_range operator[](vertices_size_type x){
//
//    }
private:
    G const& _g;
    N const& _n;
    O const& _o;
};

}


} // treedec

namespace boost{

template<class VD, class B>
std::pair<typename std::vector< std::pair<VD, B> >::const_iterator,
          typename std::vector< std::pair<VD, B> >::const_iterator >
vertices( std::vector< std::pair<VD, B> > const& skel )
{ untested();
    return std::make_pair(skel.begin(), skel.end());
}

template<class VD, class B>
VD get(vertex_owner_t, std::pair<VD, B> const& p,
         std::vector< std::pair<VD, B> > const& )
{
    return p.first;
}

template<class VD, class B>
B get(::treedec::bag_t, std::pair<VD, B> const& p,
         std::vector< std::pair<VD, B> > const& )
{
    return p.second;
}

template<class VD, class B>
B get(::treedec::bag_t, size_t u,
        std::vector< std::pair<VD, B> > const& X )
{
    return X[u].second;
}

} // boost

namespace treedec{

namespace detail{

// NOTE: G must be the fill in graph with respect to O (that is: N_G(v) = bag[v] in the algo above)
// if B is not provided
template <typename G_t, typename T_t, typename B_t, typename N_t>
class skeleton_helper{
public:
    skeleton_helper(G_t const &G, T_t &T, B_t const &B, N_t const &numbering, unsigned n)
      : _g(G), _t(T), _b(B), _numbering(numbering), _n(n)
    { untested();
    }

    template <typename X_t>
    void bag_to_treedec(std::set<X_t> const &b, T_t &T, unsigned idx){
        bag(idx, T) = MOVE(b);
    }

    template <typename X_t>
    void bag_to_treedec(std::vector<X_t> const &b, T_t &T, unsigned idx){
        for(auto bIt = b.begin(); bIt != b.end(); bIt++){
            insert(bag(idx, T), *bIt);
        }
    }

    void do_it(){
        if(_n == 0){
            return;
        }else{
        }

        //Bag for the u-th elimination vertex will be stored in T[u].
        for(unsigned u = 0; u < _n; u++){
            boost::add_vertex(_t);
        }

        //Since we made the neighbourhood N of the u-th vertex a clique,
        //the bag of the neighbour of this vertex with lowest elimination index
        //will have N as a subset.
        unsigned max = _n-1u;
        for(unsigned u = 0; u < max; u++){
            unsigned min_index = max; //note: if there's an empty bag, we can glue
                                      //it toghether with an arbitrary bag.

            auto b=boost::get(treedec::bag_t(), u, _b);
            for(auto bIt : b ){
                unsigned pos = get_pos(bIt, _g);
                unsigned index = _numbering.get_position( /*idmap?!*/ pos);
                if(index < min_index){
                    min_index = index;
                }else{
                }
            }
            //(min_index, u) will lead to a connected directed graph, if G_t is
            //directed.
            boost::add_edge(min_index, u, _t);
        }

        //Bag for the u-th elimination vertex will be stored in T[u].
        auto p=::boost::vertices(_b);
        size_t u=0;
        for(; p.first!=p.second; ++p.first){
            auto b=boost::get(treedec::bag_t(), *p.first, _b);
            auto v=boost::get(boost::vertex_owner, *p.first, _b);
            bag_to_treedec(b, _t, u);
            insert(bag(u, _t), v);
            ++u;
        }
    }

private:
    G_t const &_g;
    T_t &_t;
    B_t const &_b;
    N_t const &_numbering;
    unsigned _n;
}; // skeleton_helper

#if 0 // WIP
template <typename G_t, typename T_t, typename B_t, typename O_t>
void graph_and_numbering_to_treedec(G_t &const G, B_t const& B, N const& numbering)
{ untested();
    skeleton_helper<G_t, T_t, B_t, O_t> S(G, T, B, O, n_);
    S.do_it();
}
#endif

// create a treedecomposition of G from n_ bags B[0] ... B[n-1],
// but also put _o[i] into the ith bag.
// this is probably obsolete.
template <typename G_t, typename T_t, typename B_t, typename O_t>
void skeleton_to_treedec(G_t const &G, T_t &T, B_t const &B, O_t const &O, unsigned n_)
{ untested();

    // turn that into a numbering...
    typedef draft::NUMBERING_1<G_t> numbering_type;
    draft::NUMBERING_1<G_t> n(boost::num_vertices(G));
    for(unsigned i=0; i<n_; ++i){ itested();
        n.put(O[i]);
        n.increment();
    }
#ifndef NDEBUG
    auto b=boost::vertices(G);
    for(; b.first!=b.second; ++b.first){
        if(n.is_numbered(*b.first)){
            assert(n.get_position(*b.first)<n_);
            assert(O[n.get_position(*b.first)]==*b.first);
        }else{
        }
    }
#endif

    skeleton_helper<G_t, T_t, B_t, numbering_type> S(G, T, B, n, n_);
    S.do_it();
}

} //namespace detail

} // namespace treedec

#endif // guard

// vim:ts=8:sw=4:et
