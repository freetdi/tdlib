// Felix Salfelder, 2021, 2022
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
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
#ifndef TREEDEC_ELIM_UTIL_H
#define TREEDEC_ELIM_UTIL_H

#include "treedec_traits.hpp"

namespace treedec {
namespace draft {

template<class G, class O>
class io_smaller_than{
    typedef typename boost::property_map< G, boost::vertex_index_t >::const_type::value_type vertex_index_type;
public:
    explicit io_smaller_than(size_t k, O const& o, G const& g) : _k(k), _o(o), _g(g) {
    }
    template<class E>
    bool operator()(E const& e) const{
        auto t = boost::target(e, _g);
        auto idm = boost::get(boost::vertex_index, _g);
        if(idm[_o[t]] < _k){
            return true;
        }else{
            return false;
        }
    }
private:
    vertex_index_type _k;
    O const& _o;
    G const& _g;
};

template<class I, class N, class G, class O, class M>
void cleanup_bmdo(I j, N const& numbering, G& g, O const&, M const& my_numbering_order)
{
    io_smaller_than<G, O> P(numbering[j], numbering, g);
    vertex_descriptor_G b(j);
    boost::remove_out_edge_if(b, P, g);

    std::sort(g->vertices()[j].begin(), g->vertices()[j].end(), my_numbering_order);
}

template<class N>
class mapped_order{
public:
    mapped_order(N const& n):_n(n){}

    template<class V, class W>
    bool operator()(V const& a, W const& b) const{
        return _n[a] < _n[b];
    }
private:
    N const& _n;
};

template<class T>
void set(std::vector<T>& v, size_t idx, T& val)
{
	v[idx] = val;
}
template<class N>
void set(N& v, size_t idx, typename N::value_type val)
{
	incomplete();
}

template <typename G, typename O_t, class T, class N>
void inplace_bmdo_tree(G &g, O_t const& O, T& t, size_t bagsize, N const& io_)
{ untested();
    auto const* io = &io_;
//    typedef typename boost::property_map<G, boost::vertex_index_t>::type::value_type vertex_index_type;
    size_t num_vert = boost::num_vertices(g);

    if(num_vert == 0){ untested();
        boost::add_vertex(t);
    }else{ untested();

        assert(num_vert == O.size());

#if 0 // later.
        N iOlocal;
        if(io){ untested();
            trace2("DBG", io->size(), num_vert);
            assert(io->size()==num_vert);
            for(vertex_index_type i = 0; i < num_vert; i++){ untested();
                // bug? does O map to vertex_descriptors?
                assert(vertex_index_type((*io)[O[i]]) == i);
            }
        }else{ untested();
            iOlocal.resize(num_vert);
            io=&iOlocal;
            for(unsigned i = 0; i < num_vert; i++){ untested();
					set(iOlocal, O[i], i);
//                iOlocal[O[i]] = i;
            }
        }
        //    O_t const& iO=*io;
#else
		  assert(io);
#endif
        N const& numbering=*io;

        auto invalid = num_vert;
        std::vector<unsigned> edges(num_vert-1u, invalid);
        assert(edges.size()==num_vert-1);

        mapped_order<N> my_numbering_order(numbering);

        for(unsigned j = 0; j < num_vert-bagsize; j++){ untested();
            cleanup_bmdo(O[j], numbering, g, O, my_numbering_order);
        }

        std::vector<vertex_descriptor_G> buf;
        auto nodes_left = O.size();

        for(auto oi_ = O.begin(); ; ++oi_){ untested();
            auto oi = *oi_;
            --nodes_left;
            auto i = numbering[oi];
            auto R = boost::adjacent_vertices(oi, g);
            // auto D = boost::out_degree(oi, g);

            for(;R.first!=R.second;++R.first) { untested();
                auto j = *R.first;
                unsigned iO_n_node = numbering[j];
                if(iO_n_node < edges[i]){ untested();
                    edges[i] = iO_n_node;
                }else{ untested();
                }
            }

            if(bagsize == nodes_left + 1){ untested();
                // the rest is one big bag, no matter what.
                auto& last_adj = g->vertices()[oi];
                last_adj.resize(bagsize-1);
                size_t f = 0;
                while(++oi_ != O.end()){ untested();
                    last_adj[f++] = *oi_;
                }
                assert(last_adj.size() + 1 == bagsize);
                break;
            }else{ untested();
            }

            auto const& NN = g->vertices()[oi];
            // NN is now "bag O[i]"
            // R = boost::adjacent_vertices(oi, g);

            // rewire bag at O[i]. will become bag i.
            for( auto j_ = NN.begin(); j_ != NN.end(); ++j_ ) { untested();
                auto k_ = j_;
                ++k_;
                auto j = *j_;

                buf.resize(0);

                auto Aj = boost::adjacent_vertices(j, g);

                // could try canonical order... but then NN is wrong.
                //
                if(numbering[j] + bagsize < num_vert){ untested();
                    std::set_union(k_, NN.end(), Aj.first, Aj.second, std::back_inserter(buf), my_numbering_order);
                    std::swap(g->vertices()[j], buf);
                }else{ untested();
                }

                for( ; k_!=NN.end() ; ++k_ ){ untested();
                    auto k = *k_;
                    assert(numbering[j] < numbering[k]);

                    if((unsigned)numbering[k] < edges[numbering[j]]){ untested();
                        edges[numbering[j]] = numbering[k];
                    }else{ untested();
                    }
                }
            }
        }

        trace2("loop done", nodes_left, num_vert);
        assert(! boost::num_vertices(t));

        for(unsigned i = 0; i < num_vert-nodes_left; ++i){ untested();
            boost::add_vertex(t);
            assert(i+1 == boost::num_vertices(t));
            auto& b = boost::get(treedec::bag_t(), t, i);

            auto& NN = g->vertices()[O[i]];
            assign(b, std::move(NN));
            push(b, O[i]);
        }

        assert(boost::num_vertices(t) == num_vert-nodes_left);

        // invert edge direction?
        for(unsigned i = 0; i < num_vert-nodes_left-1u; ++i){ untested();
            assert(edges[i]>i || edges[i]==invalid);
            if(edges[i] >= num_vert-nodes_left){ untested();
                edges[i] = num_vert-nodes_left-1;
                boost::add_edge(i, edges[i], t);
            }else if(edges[i]!=invalid){ untested();
                // normal edge, as computed above.
                boost::add_edge(i, edges[i], t);
            }else if(i+1!=num_vert){ untested();
                // edge to next component
                boost::add_edge(i, i+1, t);
            }else{ untested();
                incomplete(); // ?
                trace2("edging exit", i, edges[i]);
                // exiting last connected component.
                // ... dont connect
            }
        }

        assert(boost::num_vertices(t) == num_vert-nodes_left);
//        assert(boost::num_edges(t) +1 == boost::num_vertices(t));

        for(unsigned i = 0; i < num_vert-nodes_left; i++){ untested();
            auto& b=boost::get(treedec::bag_t(), t, i);
//            assert(b.size());
            treedec::sort(b); // bug in check_tree_decomp, need sorted bags...
        }
    }

} // inplace_bmdo

}
}

#endif
