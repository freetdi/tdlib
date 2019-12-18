// Felix Salfelder 2016-2017
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
// tree decomposition operations

#ifndef TREEDEC_MISC_H
#define TREEDEC_MISC_H

#include "treedec.hpp"

namespace treedec{
namespace draft{

    // do I need a class here?!
template<class G, class M>
class applyGmap{//
public:
    applyGmap(G const& g, M const& m) : _g(g), _m(m)
    {
    }

    typedef unsigned result_type;

    unsigned operator()(unsigned x) const
    {
        auto p=boost::get(boost::vertex_index, _g, x);
        if(p == x){
            // simple numbering perhaps, needs more thought
        }else{ untested();
        }
        return _m[p];
    }
private:
    G const& _g;
    M const& _m;
};

// append one decomposition to another.
// (take care of node mapping, mainly)
// BUG: use boost::graph_copy with a sensible vertex_copier!!!1
template <typename T_t, class S_t, class G_t, class M_t>
void append_decomposition(T_t &tgt, S_t const&& src, G_t const& /*GR*/, M_t const& map)
{

	// find a leaf in tgt (hack)
	auto back=boost::vertices(tgt).first + (boost::num_vertices(tgt) - 1);
	if(boost::out_degree(*back, tgt)){
		back = boost::vertices(tgt).first;
	}else{
	}


#ifdef TREEDEC_DEBUG
	std::cout << " appending\n";
	boost::print_graph(src);
	treedec::draft::dump_tree_decomposition(src);
	std::cout << " ---- \n";
#endif
	typedef typename boost::graph_traits<T_t>::vertex_descriptor treenode_descriptor;
   if(boost::is_directed(src)){ itested();
		if(boost::is_directed(tgt)){
		}else{
		}
	}else{ itested();
		if(boost::is_directed(tgt)){ untested();
		}else{
		}
	}
#ifdef DEBUG
	for(auto& i: map){ untested();
		std::cerr << "map " << i << "\n";
	}
#endif
//        std::cout<<"tgt\n";
  //      boost::print_graph(tgt);
    //    std::cout<<"appending\n";
      //  boost::print_graph(src);
    assert(boost::num_vertices(tgt) == boost::num_edges(tgt)+1);
    assert(boost::num_vertices(src) == boost::num_edges(src)+1);
    auto op=applyGmap<S_t, M_t>(src, map);
//    assert(undirected_tree)...
    unsigned offset=boost::num_vertices(tgt);
    unsigned sn=boost::num_vertices(src);
    auto SR=boost::vertices(src);
    auto TR=boost::vertices(tgt);
    trace4("append", offset, SR.second-SR.first, TR.second-TR.first, map.size());

    if(SR.first==SR.second){ untested();
        // no bags to fetch. done
    }else{
		 treenode_descriptor new_tv;
		  auto const& M=boost::get(bag_t(), src);
		  unsigned n=0;
		  std::vector<treenode_descriptor> tmap(sn);
        for(; SR.first!=SR.second; ++SR.first){
            new_tv = boost::add_vertex(tgt);
				tmap[n] = new_tv;
				++n;
		  }
		  n = 0;
		  SR = boost::vertices(src);
        for(; SR.first!=SR.second; ++SR.first) {
            trace1("appendfound", *SR.first);
            new_tv=tmap[n];
            auto& B=boost::get(bag_t(), tgt, new_tv);
            auto const& SB=boost::get(M, *SR.first);

            assert(boost::degree(*SR.first, src)>0);
//            assert(boost::out_degree(*SR.first, src)>0); //?
            auto srcnp=boost::adjacent_vertices(*SR.first, src);
				assert(srcnp.first!=srcnp.second || boost::is_directed(src));
				for(;srcnp.first!=srcnp.second; ++srcnp.first){
					auto target_in_copy = *srcnp.first;
					assert(target_in_copy!=*SR.first);
					if(boost::is_directed(src)){
						// assert(target_in_copy+offset <= new_tv);
						boost::add_edge(new_tv, target_in_copy+offset, tgt);
					}else if(target_in_copy<*SR.first){
						trace2("copy edge", *SR.first, target_in_copy);
						assert(!boost::edge(new_tv, target_in_copy+offset, tgt).second);
						// this one seems correct?
						boost::add_edge(target_in_copy+offset, new_tv, tgt);
						//boost::add_edge(new_tv, target_in_copy+offset, tgt);
					}else{
						/// ???
						trace1("no copy edge", *SR.first);
					}
				}

#if 0 // vector version (how to select?)
            B = std::move(bag(*SR.first, src));
            // now apply map
            //
            for(auto & x: B){ untested();
                x = map[get_pos(x,src)];
            }
#else

            using draft::applyGmap;
//            using boost::bind1st;
            using boost::make_transform_iterator;
            push(B, make_transform_iterator(SB.begin(), op ),
                    make_transform_iterator(SB.end(),   op ));
#endif
				++n;
        }

        if(offset){
            // connect new stuff to existing.
				trace2("connecting existing ", new_tv, *boost::vertices(tgt).first);
            assert(!boost::edge(new_tv, *boost::vertices(tgt).first, tgt).second);
            //boost::add_edge(new_tv, *boost::vertices(tgt).first, tgt);
				//
				// hack/heuristic: find a root in src
				// should use rooted graphs instead
				// (currently, the root is either front or back in all outputs...)
				if(boost::out_degree(tmap.back(), tgt)){
					boost::add_edge(*back, tmap.back(), tgt);
				}else{
					boost::add_edge(*back, tmap.front(), tgt);
				}
        }else{ untested();
        }
    }
	 trace2("done", boost::num_vertices(src), boost::num_edges(src));
	 trace2("done", boost::num_vertices(tgt), boost::num_edges(tgt));
    assert(boost::num_vertices(tgt) == boost::num_edges(tgt)+1);
//        boost::print_graph(tgt);
} // append_decomposition

#ifdef NDEBUG
#define assert_permutation(P)
#else
#define assert_permutation(P) \
{ \
	std::vector<bool> check(P.size()); \
		for(unsigned i : P){ \
		    check[i]=true; \
		} \
		for(auto i : check){ \
		    assert(i); \
		} \
}
#endif

template<class V, class P>
void permute_vector(V& vec, P &perm)
{ untested();
	assert_permutation(perm);
	unsigned seek=0;
	unsigned n=vec.size();

	for(; seek<n;){ untested();

		if(perm[seek]<=seek){ untested();
			// already done.
		}else{ untested();
			typename V::value_type tmp = vec[seek];
			unsigned i = perm[seek];
			while(i!=seek){ untested();
				std::swap(tmp, vec[i]);
				perm[i]=seek;
				i=perm[i];
			}
			vec[i]=tmp;
		}
		++seek;
	}
}

template<class V, class P>
void concat_maps(V& first, P const&second){ untested();
	unsigned n=first.size();
	for(unsigned i=0; i<n; ++i){ untested();
		first[i]=second[first[i]];
	}
}

} // draft


} // treedec

#endif
