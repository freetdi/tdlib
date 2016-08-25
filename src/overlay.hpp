// Felix Salfelder, 2015 - 2016
//
// (c) 2016 Goethe-Universität Frankfurt
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
// graph overlays through views or (partial) copies.
//

#ifndef OVERLAY_H
#define OVERLAY_H

namespace treedec {
namespace draft {

// immutable overlay
// make a vertex subset of a graph look like a graph
// this graph is immutable, to allow for efficient storage
// a callback can be used to add even more edges.
template<class G_t, class I_t, class S_t, class IG_t, class M_t, class CB_t>
inline IG_t const& immutable_clone(
     G_t const &G,
     IG_t& ig,
     I_t bbegin,
     I_t bend,
     S_t bag_nv,
    //   URGHS. no default types without c++11.
     M_t* vdMap, /*=NULL*/
     CB_t* cb
     )
{
//    typedef typename graph_traits<G_t>::immutable_type immutable_type;
    typedef typename boost::graph_traits<IG_t>::vertex_descriptor vertex_descriptor_ig;

    BOOST_AUTO(nv, boost::num_vertices(G));
    ig = MOVE(IG_t(bag_nv));

    assert(bag_nv == boost::num_vertices(ig));

    // map ig vertices (positions) to bag elements (= vertices in G)
    M_t local_vd_map;
    // std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> local_vd_map;
    if(vdMap){
        // use that...
    }else{ untested();
        vdMap = &local_vd_map;
    }
    vdMap->resize(bag_nv);
    // map vertex positions in G to vertices in ig
    std::vector<vertex_descriptor_ig> reverse_map(nv);

    BOOST_AUTO(bi, bbegin);
    BOOST_AUTO(be, bend);
    unsigned i=0;
    for(; bi!=be; ++bi){
        // FIXME: pos, vertex_index?
        assert(i < vdMap->size());
        (*vdMap)[i] = *bi;
        reverse_map[get_pos(*bi, G)] = i;
        ++i;
    }
    assert(i==bag_nv);


    bi = bbegin;
    unsigned s=-1;
    unsigned t=-1;
    unsigned vertices_count;
    for(; bi!=be; ++bi){
        ++vertices_count;
        
        if(!cb){
            BOOST_AUTO(s, get_pos(*bi, G));
            auto A=boost::adjacent_vertices(*bi,G);
            for(;A.first!=A.second;++A.first){
                BOOST_AUTO(t, get_pos(*A.first, G));
                boost::add_edge(reverse_map[s], reverse_map[t], ig);
            }
        }else{
            BOOST_AUTO(vi, bi);
            ++vi; // skip self loop

            for(; vi!=be; ++vi){
                bool edg=false;
                if(boost::edge(*bi, *vi, G).second){
                    edg = true;
                }else if(!cb){
                }else if((*cb)(*bi, *vi)){
                    edg = true;
                }else{
                    // no edge.
                }

                if(edg){
                    BOOST_AUTO(s, get_pos(*bi, G));
                    BOOST_AUTO(t, get_pos(*vi, G));
                    boost::add_edge(reverse_map[s], reverse_map[t], ig);
                }else if(s==-1u){
                    assert(get_pos(*bi, G)!=-1u);
                    s = get_pos(*bi, G);
                    t = get_pos(*vi, G);
                }else{
                }
            }
        }
    }
    // HACK. not here.
    if(cb && s!=-1u){
        /// let MSVS know about a particular new edge
        cb->a = reverse_map[s];
        cb->b = reverse_map[t];
    }else{
        // assert(is_clique(ig));
    }

    return ig;
}

// FIXME: must be more implicit...
namespace dummy_hack{
template<class VD_t>
class cb{ //
public:
    cb(){unreachable();}
    bool operator()(VD_t, VD_t){unreachable(); return false;}
public: // HACK
    unsigned a, b;
};
}

template<class G_t, class I_t, class S_t, class IG_t, class M_t>
inline IG_t const& immutable_clone(
     G_t const &G,
     IG_t& ig,
     I_t bbegin,
     I_t bend,
     S_t bag_nv,
     M_t* vdMap /*=NULL*/)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
    dummy_hack::cb<vd>* c=NULL;
    return immutable_clone(G, ig, bbegin, bend, bag_nv, vdMap, c);
}

} // draft

} // treedec
#endif // guard
// vim:ts=8:sw=4:et
