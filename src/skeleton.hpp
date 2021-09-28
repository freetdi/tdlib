// Lukas Larisch, 2017
// Felix Salfelder, 2017
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
//
// almost-tree-decompositions
// a skeleton is a list of vertices with vertices attached to it.

#ifndef TREEDEC_SKELETON_HPP
#define TREEDEC_SKELETON_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "numbering.hpp"
#include "treedec_traits.hpp"

namespace treedec{

namespace draft{

template<class G, class N, class O>
class SKELETON{
public: // types
    typedef typename boost::graph_traits<G>::adjacency_iterator all_adj_it;
    typedef unsigned vertex_descriptor;
    typedef typename boost::graph_traits<G>::vertex_descriptor owner_type;
    typedef typename boost::graph_traits<G>::vertices_size_type vertices_size_type;
    struct sgm{
        sgm(N const& n, owner_type v)
          : _n(n), _x(v)
        {
        }
        bool operator()(owner_type v) const{
            // idea: add to bag, if not eliminated earlier.
#ifdef DEBUG
            if( !_n.is_before(_x, v)){
                trace1("skip", v);
            }else{
                trace1("noskip", v);
            }
#endif

            return _n.is_before(_x, v);
        }
        N const& _n;
        owner_type _x;
    };
    typedef boost::filter_iterator<sgm, all_adj_it> bag_iterator;
    typedef boost::iterator_range<bag_iterator> bag_type;
public: // construct
    SKELETON(G const& g, N const& n, O const& o)
      : _g(g), _n(n), _o(o)
    {
    }
public:
    size_t size() const{
        trace2("skeleton::size", _n.total(), _o.size());
        assert(_n.total() == _o.size());
        return _n.total();
        return _o.size();
    }
//    bag_range operator[](vertices_size_type x){
//
//    }
    boost::iterator_range<bag_iterator> bag(vertex_descriptor i) const{
        assert(i<_o.size());
        auto p=boost::adjacent_vertices(_o[i], _g);

        typedef boost::filter_iterator<sgm, all_adj_it> FilterIter;
        sgm fP(_n, _o[i]);

        FilterIter fb(fP, p.first, p.second);
        FilterIter fe(fP, p.second, p.second);

        return  make_iterator_range(std::make_pair(fb, fe));
    }
    owner_type owner(unsigned u) const {
        return _o[u];
    }
private:
    G const& _g;
    N const& _n;
    O const& _o;
}; // SKELETON

}


} // treedec

namespace boost{

template<class VD, class B>
size_t num_vertices( std::vector< std::pair<VD, B> > const& skel )
{
    trace1("boost wrap", skel.size());
    return skel.size();
}

template<class VD, class B>
std::pair<typename std::vector< std::pair<VD, B> >::const_iterator,
          typename std::vector< std::pair<VD, B> >::const_iterator >
vertices( std::vector< std::pair<VD, B> > const& skel )
{
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

template<class A, class B, class C>
typename treedec::draft::SKELETON<A, B, C>::bag_type
get(::treedec::bag_t, size_t u, treedec::draft::SKELETON<A, B, C> const& X )
{
    return X.bag(u);
}

template<class A, class B, class C>
typename treedec::draft::SKELETON<A, B, C>::owner_type
get(vertex_owner_t, size_t u, treedec::draft::SKELETON<A, B, C> const& X )
{
    return X.owner(u);
}

template<class A, class B, class C>
typename treedec::draft::SKELETON<A, B, C>::vertices_size_type
num_vertices( treedec::draft::SKELETON<A, B, C> const& skel )
{
    return skel.size();
}

template<class A, class B, class C>
std::pair<boost::counting_iterator<unsigned>,
          boost::counting_iterator<unsigned> >
vertices( treedec::draft::SKELETON<A, B, C> const& skel )
{
    typedef counting_iterator<unsigned> ci;
    return std::make_pair(ci(0), ci(skel.size()));
}

} // boost

namespace treedec{

namespace detail{

// G is the fill in graph with respect to O (that is: N_G(v) = bag[v] in the algo above)
// if B is not provided
template <typename G_t, typename T_t, typename B_t, typename N_t>
class skeleton_helper{
public:
    skeleton_helper(G_t const& G, T_t &T, B_t const &B, N_t const &numbering)
      : _g(G), _t(T), _b(B), _numbering(numbering)
    {
    }

private: // impl
    template <typename X_t, class B>
    void bag_to_treedec(std::set<X_t> const &b, B &T){
        T = MOVE(b);
    }

    template <typename X_t, class B>
    void bag_to_treedec(std::vector<X_t> const &b, B &T){
        for(auto bIt = b.begin(); bIt != b.end(); bIt++){
            insert(T, *bIt);
        }
    }
    template <typename X_t, class B>
    void bag_to_treedec(boost::iterator_range<X_t> const &b, B &T){
        for(auto x : b){
            push(T, x);
        }
    }
public:

    void do_it(){
        if(boost::num_vertices(_b) == 0){ untested();
            return;
        }else if(boost::num_vertices(_b) == boost::num_vertices(_t)){
        }else if(!boost::num_vertices(_t)){
            //Bag for the u-th elimination vertex will be stored in T[u].
            for( size_t x=0 ; x<_numbering.total(); ++x){
                boost::add_vertex(_t);
            }
        }else{
            trace2("skeleton::do_it mismatch", boost::num_vertices(_b), boost::num_vertices(_t));
            incomplete();
        }

        //Since we made the neighbourhood N of the u-th vertex a clique,
        //the bag of the neighbour of this vertex with lowest elimination index
        //will have N as a subset.
        unsigned max = _numbering.total()-1u; // boost::num_vertices(_b)-1u;

#if 0 // swapping ordering (HACK). td format... (incomplete tdprinter)
        for(unsigned u = 0; u < max; u++){
            unsigned min_index = max; //note: if there's an empty bag, we can glue
                                      //it together with an arbitrary bag.

            auto b=boost::get(treedec::bag_t(), u, _b);
            for(auto bIt : b ){
                auto index=_numbering.get_position(bIt);

                // BUG: get_position applies map
                // auto pos=boost::get(boost::vertex_index, _g, bIt);
                // unsigned index = _numbering.get_position( /*idmap?!*/ pos);
                if(index < min_index){
                    min_index = index;
                }
            }
            //(min_index, u) will lead to a connected directed graph, if G_t is
            //directed.
            boost::add_edge(min_index, u, _t);
        }
#endif

        //Bag for the u-th elimination vertex will be stored in T[u].
        auto p=::boost::vertices(_b);
        size_t u=0; // does not work.
//
        for(; p.first!=p.second; ++p.first){
            auto v=boost::get(boost::vertex_owner, *p.first, _b);

            if(_numbering.is_numbered(v)){
                auto b=boost::get(treedec::bag_t(), *p.first, _b);
                auto& target_bag=boost::get(treedec::bag_t(), _t, u);
                bag_to_treedec(b, target_bag);
                push(target_bag, v); // insert?
                ++u;
             }else{ untested();
             }
        }
        assert(u==_numbering.total());

        for(unsigned u = 0; u < max; u++){
            unsigned min_index = max; //note: if there's an empty bag, we can glue
                                      //it toghether with an arbitrary bag.

            auto b=boost::get(treedec::bag_t(), u, _b);
            for(auto bIt : b ){
                auto index=_numbering.get_position(bIt);

                if(index < min_index){
                    min_index = index;
                }else{
                }
            }
            //(min_index, u) will lead to a connected directed graph, if G_t is
            //directed.
            boost::add_edge(min_index, u, _t);
        }
    }

private:
    G_t const &_g;
    T_t &_t;
    B_t const &_b;
    N_t const &_numbering;
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
{

    // turn that into a numbering...
    typedef draft::NUMBERING_1<G_t> numbering_type;
    draft::NUMBERING_1<G_t> n(boost::num_vertices(G));
    for(unsigned i=0; i<n_; ++i){
        n.put(O[i]);
        n.increment();
    }
#ifndef NDEBUG
    auto b=boost::vertices(G);
    for(; b.first!=b.second; ++b.first){
        if(n.is_numbered(*b.first)){
            assert(n.get_position(*b.first)<n_);
            assert(O[n.get_position(*b.first)]==*b.first);
        }
    }
#endif

    assert(n_ == boost::num_vertices(B));
    skeleton_helper<G_t, T_t, B_t, numbering_type> S(G, T, B, n);
    S.do_it();
}

} //namespace detail

} // namespace treedec

#endif // guard

// vim:ts=8:sw=4:et
