// Felix Salfelder 2016
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
// exact decomposition, common denominator
//

// does not work (yet?) because the size is not known in advance.
// ... immutable_clone needs multipass iterator (?)
// #define NEWRANGE
// ( TODO: move to pp_base, this way it does not work )

#ifndef TREEDEC_EXACT_BASE_HPP
#define TREEDEC_EXACT_BASE_HPP

#include "lower_bounds.hpp"
#include "preprocessing.hpp"
#include "overlay.hpp"
#include "treedec_misc.hpp"
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/bandwidth.hpp>
#include "treedec.hpp"
#include <boost/graph/copy.hpp>

namespace treedec{

namespace draft{

template<typename G_t,
	template<class G_, class ...> class config,
	template<class X, template<class G__, class ...> class Y> class kernel>
class exact_decomposition /* algo? */ {
public:
  typedef kernel<G_t, config> kern_t;
public:
  template<class GraphType, class... rest>
  struct cfg16 : config<GraphType>{
      static constexpr unsigned max_vertex_index=15;
      // sparsity?
      // bandwidth?
  };
  template<class GraphType, class... rest>
  struct cfg32 : config<GraphType>{
      static constexpr unsigned max_vertex_index=31;
  };
  template<class GraphType, class... rest>
  struct cfg64 : config<GraphType>{
      static constexpr unsigned max_vertex_index=63;
  };
  template<class GraphType, class... rest>
  struct cfg128 : config<GraphType>{
      static constexpr unsigned max_vertex_index=127;
  };
  template<class GraphType, class... rest>
  struct cfg192 : config<GraphType>{
      static constexpr unsigned max_vertex_index=191;
  };
  template<class GraphType, class... rest>
  struct cfg256 : config<GraphType>{
      static constexpr unsigned max_vertex_index=255;
  };
  template<class GraphType, class... rest>
  struct cfg512 : config<GraphType>{
      static constexpr unsigned max_vertex_index=511;
  };
  template<class GraphType, class... rest>
  struct cfg1024 : config<GraphType>{
      static constexpr unsigned max_vertex_index=1023;
  };
public:
    exact_decomposition(G_t &G)
        : _g(G), _cleanup_g(false)
    {
    }
    exact_decomposition(G_t const &G)
        : _g( *new G_t(G)), _cleanup_g(true)
    { untested();
    }
    ~exact_decomposition() {
        if(_cleanup_g){ untested();
            delete &_g;
        }else{
        }
    }
public:
    void do_it(unsigned lb_bs=0){
        try_it(_t, lb_bs);
    }
    template<class T>
    void get_tree_decomposition(T& t) const{
        boost::copy_graph(_t, t);
    }
    template<class T>
    void try_it(T&, unsigned lb_bs);
    template<class G, class T>
    void run_kernel(G const&, T&, unsigned& lb_bs);
private:
    template<class T_t>
    void do_components(T_t&, unsigned lb_bs);
private:
    G_t& _g;
    typename graph_traits<G_t>::treedec_type _t;
    bool _cleanup_g;
};


template<typename G_t,
	template<class G_, class ...> class config,
	template<class X, template<class G__, class... > class Y> class kernel>
template<class G_t_, class T_t>
inline void exact_decomposition<G_t, config, kernel>::run_kernel(
        G_t_ const& G, T_t& T, unsigned& lb_bs)
{

    auto numv=boost::num_vertices(G);
    // std::cout << "c kernel " << numv << "\n";
    if(numv<=32){
        kernel<G_t_, cfg32> kern(G);
        kern.do_it(T, lb_bs);
    }else if(numv<=64){ untested();
        kernel<G_t_, cfg64> kern(G);
        kern.do_it(T, lb_bs);
    }else if(numv<=128){ untested();
        kernel<G_t_, cfg128> kern(G);
        kern.do_it(T, lb_bs);
    }else if(numv<=192){ untested();
        kernel<G_t_, cfg192> kern(G);
        kern.do_it(T, lb_bs);
    }else if(numv<=256){ untested();
        kernel<G_t_, cfg256> kern(G);
        kern.do_it(T, lb_bs);
    }else if(numv<=512){ untested();
        kernel<G_t_, cfg512> kern(G);
        kern.do_it(T, lb_bs);
    }else{ incomplete();
        kernel<G_t_, cfg1024> kern(G);
        kern.do_it(T, lb_bs);
    }
    assert(boost::num_vertices(T) == boost::num_edges(T)+1);
}

template<typename G_t,
	template<class G_, class ...> class config,
	template<class X, template<class G__, class ...> class Y> class kernel>
template<class T_t>
inline void exact_decomposition<G_t, config, kernel>::do_components(
        T_t& t, unsigned lb_bs)
{

    // Compute a tree decomposition for each connected component of G and glue
    // the decompositions together.
    // FIXME: move to preprocessor.
#ifndef NEWRANGE
    typedef std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components_t;
    components_t components;
    treedec::get_components(_g, components);
#else
    auto R=boost::vertices(_g);
    std::vector<EXCUT_BOOL> visited(n);
    auto cmps_range = make_components_range(
                R.first, R.second,
                _g, &visited, NULL, EXCUT_BOOL());
#endif

    // root
    boost::add_vertex(t);
    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
#ifndef NEWRANGE
    typename components_t::iterator i = components.begin();
    for(; i!=components.end(); ++i)
#else
    for(; cmps_range.first != cmps_range.second; ++cmps_range.first)
#endif
    {
        //Ignore isolated vertices (already included in 'bags').
#ifndef NEWRANGE
        trace2("found component ", i->size(), components.size());
        if(i->size() == 1){
            continue;
            auto nv=boost::add_vertex(t);
            auto& B=boost::get(bag_t(), t, nv);
            treedec::push(B, *(*i).begin());
            trace2("isolated node ", nv,  *(*i).begin());
            if(nv!=0){ untested();
                // uuh hack
                boost::add_edge(nv, nv-1, t);
            }else{ untested();
            }
        }else{
        }
#endif

        typedef typename graph_traits<G_t>::immutable_type immutable_type;

#ifndef NEWRANGE
        unsigned compsize = i->size();
        auto comp_range = *i;
#else
        unsigned compsize = n; // don't know... incomplete!!
        auto comp_range = *(cmps_range.first);
#endif
        immutable_type H(compsize);
        const immutable_type& G_=draft::immutable_clone(_g, H,
                std::begin(comp_range), std::end(comp_range), compsize,
                &vdMap);

#ifdef NEWRANGE
        if(boost::num_vertices(H)==1){
            incomplete();
        }
#endif
        trace3("excut comp", lb_bs, boost::num_vertices(G_), boost::num_edges(G_));

        assert_connected(G_);

//        std::cout << "raw.\n";
//        print_matrix(G_, boost::identity_property_map(), std::cout);

#ifdef USE_RCMK
        incomplete();
	typedef typename boost::graph_traits<immutable_type>::vertex_descriptor Vertex;
	typedef std::vector<Vertex> M;
	auto nv=boost::num_vertices(G_);

	auto ne2=boost::num_edges(G_);
	auto ne3=boost::num_edges(_g);

	M m(nv);
	M minv(nv);
        unsigned bw=rcmk_(G_, m, minv);

        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> balu;
        balu OG;

        auto G_cmk=treedec::make_mapped_graph(G_, m);
	auto ne1=boost::num_edges(G_cmk);
        assert(ne1==ne2);
#endif

        // not yet
        // rcmk_(G_);
        T_t T_;

#ifdef USE_RCMK
        incomplete();
        run_kernel(G_cmk, T_, lb_bs);

        // permute_vector(vdMap, m); // hmm
        concat_maps(minv, vdMap);
        draft::append_decomposition(t, std::move(T_), G_, minv);
#else
        run_kernel(G_, T_, lb_bs);
        assert(boost::num_vertices(T_) == boost::num_edges(T_)+1);
        assert(boost::num_vertices(t) == boost::num_edges(t)+1);

#ifndef NDEBUG
        unsigned tn=boost::num_vertices(t);
        unsigned tn_=boost::num_vertices(T_);
#endif
        draft::append_decomposition(t, std::move(T_), G_, vdMap);
        assert(boost::num_vertices(t) == boost::num_edges(t)+1);
        assert( tn+tn_ == boost::num_vertices(t));
#endif
    }
    assert(boost::num_vertices(t) == boost::num_edges(t)+1);
} // do_it

template<typename G_t,
	template<class G_, class ...> class config,
	template<class X, template<class Y, class ...> class c> class kernel>
template<class T_t>
void exact_decomposition<G_t, config, kernel>::try_it(T_t& T, unsigned lb_bs)
{
    typedef config<G_t> CFG;
    int lb_tw=lb_bs-1;

    auto n=boost::num_vertices(_g);
    auto e=boost::num_edges(_g);
    trace3("exact_decomposition", lb_tw, n, e);
    if(n==0){
        boost::add_vertex(T);
        return;
    }else{
    }

    //Preprocessing.
    // this is really inefficient...
    int low = -1;

    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags;

    // if config.preprocessing?
    // fixme: instanciate.
    treedec::preprocessing(_g, bags, low);

    if(boost::num_edges(_g) == 0){
        treedec::glue_bags(bags, T);
        return;
    }else{
    }

    CFG::message(0, "PP said tw %d\n", low);

    //Lower bound on the treewidth of the reduced instance of G.
    G_t H(_g);
    incomplete(); //deltac does not seem to work
    int tw_lb_deltaC = 0; // treedec::lb::deltaC_least_c(H);

    CFG::message(0, "deltaC said tw %d\n", tw_lb_deltaC);

    trace3("excut comb", lb_tw, low, tw_lb_deltaC);
    if(low > lb_tw){
        lb_tw = low;
    }else{ untested();
    }

    if (tw_lb_deltaC > lb_tw){ untested();
        lb_tw = tw_lb_deltaC;
    }else{
    }

    trace3("excut comb", lb_tw, boost::num_vertices(_g), boost::num_edges(_g));

    do_components(T, lb_tw+1);
    trace1("did components", bags.size());
    assert(boost::num_vertices(T) == boost::num_edges(T)+1);

    treedec::glue_bags(bags, T);
    trace2("done", boost::num_vertices(T), boost::num_edges(T));
    assert(boost::num_vertices(T) == boost::num_edges(T)+1);
} // try_it

}// draft

}// treedec

#endif
// vim:ts=8:sw=4:et
