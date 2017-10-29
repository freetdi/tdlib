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
    { itested();
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
template <typename T_t, class S_t, class G_t, class M_t>
void append_decomposition(T_t &tgt, S_t const&& src, G_t const& /*GR*/, M_t const& map)
{
#ifdef DEBUG
	for(auto& i: map){
		std::cerr << "map " << i << "\n";
	}
#endif
    assert(boost::num_vertices(tgt) == boost::num_edges(tgt)+1);
    assert(boost::num_vertices(src) == boost::num_edges(src)+1);
    auto op=applyGmap<S_t, M_t>(src, map);
//    assert(undirected_tree)...
    unsigned offset=boost::num_vertices(tgt);
    auto SR=boost::vertices(src);
    auto TR=boost::vertices(tgt);
    trace4("append", offset, SR.second-SR.first, TR.second-TR.first, map.size());

    if(SR.first==SR.second){ untested();
        // no bags to fetch. done
    }else{
        unsigned new_tv;
        auto next=SR.first;
		  auto const& M=boost::get(bag_t(), src);
        for(; SR.first!=SR.second; SR.first=next){
            trace1("appendfound", *SR.first);
            ++next;
            new_tv=boost::add_vertex(tgt);
            auto& B=boost::get(bag_t(), tgt, new_tv);
            auto const& SB=boost::get(M, *SR.first);

            assert(boost::degree(*SR.first, src)>0);
            assert(boost::out_degree(*SR.first, src)>0);
            auto srcnp=boost::adjacent_vertices(*SR.first, src);
				assert(srcnp.first!=srcnp.second);
				for(;srcnp.first!=srcnp.second; ++srcnp.first){
            auto target_in_copy = *srcnp.first;
            assert(target_in_copy!=*SR.first);
            if(target_in_copy<*SR.first){
					trace2("copy edge", *SR.first, target_in_copy);
					assert(!boost::edge(new_tv, target_in_copy+offset, tgt).second);
					boost::add_edge(new_tv, target_in_copy+offset, tgt);
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
        }

        if(offset){
            // connect new stuff to existing.
				trace2("connecting existing ", new_tv, *boost::vertices(tgt).first);
            assert(!boost::edge(new_tv, *boost::vertices(tgt).first, tgt).second);
            boost::add_edge(new_tv, *boost::vertices(tgt).first, tgt);
        }else{ untested();
        }
    }
	 trace2("done", boost::num_vertices(src), boost::num_edges(src));
	 trace2("done", boost::num_vertices(tgt), boost::num_edges(tgt));
    assert(boost::num_vertices(tgt) == boost::num_edges(tgt)+1);
}

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
{
	assert_permutation(perm);
	unsigned seek=0;
	unsigned n=vec.size();

	for(; seek<n;){ untested();

		if(perm[seek]<=seek){ untested();
			// already done.
		}else{
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
void concat_maps(V& first, P const&second){
	unsigned n=first.size();
	for(unsigned i=0; i<n; ++i){
		first[i]=second[first[i]];
	}
}

} // draft

template <typename T_t>
void thicken(T_t &T){
    unsigned int maxsize = (unsigned int) treedec::get_width(T)+1;
    bool modified = true;

    //Fill bags such that they all have size 'maxsize'.
    while(modified){
        modified = false;

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if((int)bag(*tIt, T).size() == maxsize){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                    typename treedec_traits<T_t>::bag_type::iterator bIt = bag(*tIt, T).begin();
                    while(bag(*nIt, T).size() < maxsize){
                        bag(*nIt, T).insert(*(bIt++));
                        modified = true;
                    }
                }
            }
        }
    }

    modified = true;

    //Remove duplicated bags.
    while(modified){
        modified = false;

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                if(bag(*tIt, T) == bag(*nIt, T)){
                    typename boost::graph_traits<T_t>::adjacency_iterator nIt2, nEnd2;
                    for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(*tIt, T); nIt2 != nEnd2; nIt2++){
                        if(*nIt2 != *nIt){
                            boost::add_edge(*nIt2, *nIt, T);
                        }
                    }
                    boost::clear_vertex(*tIt, T);
                    boost::remove_vertex(*tIt, T);
                    modified = true;

                    goto NEXT_ITER1;
                }
            }
        }
        NEXT_ITER1: ;
    }

    //Adjacent bags B1, B2 must fulfill |(B1 ^ B2)| = maxwidth-1.
    while(modified){
        modified = false;

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(bag(*tIt, T).begin(), bag(*tIt, T).end(),
                                      bag(*nIt, T).begin(), bag(*nIt, T).end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() != maxsize-1){
                    typename boost::graph_traits<T_t>::vertex_descriptor new_vertex = boost::add_vertex(T);
                    bag(new_vertex, T) = intersection;

                    typename treedec_traits<T_t>::bag_type::iterator bIt = bag(*tIt, T).begin();
                    while(bag(new_vertex, T).size() < maxsize){
                        bag(new_vertex, T).insert(*(bIt++));
                    }

                    boost::remove_edge(*tIt, *nIt, T);
                    boost::add_edge(*tIt, new_vertex, T);
                    boost::add_edge(new_vertex, *nIt, T);

                    modified = true;
                    goto NEXT_ITER2;
                }
            }
        }
        NEXT_ITER2: ;
    }
}

// legacy. does not make anything.
template <typename T_t>
void make_thick(T_t &T){ untested();
	return thicken(T);
}

} // treedec

#endif
