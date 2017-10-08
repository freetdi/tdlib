// Lukas Larisch, 2014 - 2016
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

#ifndef TD_VALIDATION
#define TD_VALIDATION

namespace treedec{

namespace validation{

template <typename G_t>
bool is_valid_clique(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    for(typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt1 = C.begin(); sIt1 != C.end(); sIt1++){
        typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != C.end(); sIt2++){
            if(!boost::edge(*sIt1, *sIt2, G).second){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_independent_set(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() && C.find(boost::target(*eIt, G)) != C.end()){
            return false;
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_vertex_cover(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &C){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(C.find(boost::source(*eIt, G)) != C.end() || C.find(boost::target(*eIt, G)) != C.end()){
        }
        else{
            return false;
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_dominating_set(G_t &G, typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type &set){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(set.find(*vIt) == set.end()){
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            bool hit = false;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                if(set.find(*nIt) != set.end()){
                    hit = true;
                    break;
                }
            }
            if(!hit){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_coloring(G_t &G, std::vector<typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type> &vec)
{
    std::vector<BOOL> visited(boost::num_vertices(G), false);
    for(unsigned i = 0; i < vec.size(); i++){
        for(typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt = vec[i].begin(); sIt != vec[i].end(); sIt++){
            unsigned pos = get_pos(*sIt, G);
            visited[pos] = true;
        }
    }
    for(unsigned i = 0; i < visited.size(); i++){
        if(!visited[i]){
            return false;
        }
    }

    std::vector<unsigned> col(boost::num_vertices(G));
    for(unsigned i = 0; i < vec.size(); i++){
        for(typename treedec::treedec_traits<typename treedec::treedec_chooser<G_t>::type>::bag_type::iterator sIt = vec[i].begin(); sIt != vec[i].end(); sIt++){
            col[get_pos(*sIt, G)] = i;
        }
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(col[get_pos(boost::source(*eIt, G), G)] == col[get_pos(boost::target(*eIt, G), G)]){
            return false;
        }
    }
    return true;
}

} //namespace validation

} //namespace treedec

#endif //ifdef TD_VALIDATION
