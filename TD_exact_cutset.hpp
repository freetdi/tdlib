// Lukas Larisch, 2014 - 2015
//
// (c) 2014-2015 Goethe-Universität Frankfurt
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
// Offers functionality to compute a tree decomposition of exact width. 
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, tree_dec_node> tree_dec_t;
//
// Vertices of the input graph have to provide the attribute 'id', e.g.:
//
// struct Vertex
// {
//  unsigned int id;
// };
// typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, Vertex> TD_graph_t;
//
//
//
// These functions are most likely to be interesting for outside use:
//
// void exact_cutset(G_t &G, T_t &T, int lb)
// void exact_cutset(G_t &G, T_t &T)
//

#define TEST

namespace treedec{

template <typename G_t>
bool explore_cutsets(G_t &G, std::set<typename boost::graph_traits<G_t>::vertex_descriptor> cut, std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component, std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &results, unsigned int k){

    if(cut.size() > k)
        return false;

    else if(cut.size() + component.size() <= k+1){
        component.insert(cut.begin(), cut.end());
        results.push_back(component);
        results.push_back(cut);

        return true;
    }

    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> candidates;
    
    //it is not sufficient to just consider the neighbourhood of cut in component
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = component.begin(); sIt != component.end(); sIt++)
            candidates.push_back(*sIt);

    for(unsigned int i = 0; i < candidates.size(); i++){
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> cut_ext = cut;
        cut_ext.insert(candidates[i]);

        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component_red = component;
        component_red.erase(candidates[i]);

        std::vector<bool> visited(boost::num_vertices(G), true);

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = component_red.begin(); sIt != component_red.end(); sIt++)
            visited[G[*sIt].id] = false;

        typename std::vector<typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > new_components;
        get_components_provided_map(G, new_components, visited);

        bool all_successful = true;

        for(unsigned int t = 0; t < new_components.size(); t++){
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> cut_red;

            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = cut_ext.begin(); sIt != cut_ext.end(); sIt++){
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, G); nIt != nEnd; nIt++){
                    if(new_components[t].find(*nIt) != new_components[t].end())
                        cut_red.insert(*sIt);
                }
            }

         unsigned int idx = results.size();

            if(!explore_cutsets(G, cut_red, new_components[t], results, k)){
                all_successful = false;
                results.erase(results.begin()+idx, results.end());
                break;
            }
            results.push_back(cut_red);
            results.push_back(cut_ext);
            results.push_back(cut_ext);
            results.push_back(cut);
        }

        if(all_successful)
            return true;

    }

    return false;
}

template <typename T_t>
void glue_bags(T_t &T, std::set<unsigned int> bag1, std::set<unsigned int> &bag2){
    if(bag1 == bag2)
        return;

    typename boost::graph_traits<T_t>::vertex_iterator vIt1, vIt2, vEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor b1,b2;
    bool bag1_found = false; 
    bool bag2_found = false;
    for(boost::tie(vIt1, vEnd) = boost::vertices(T); vIt1 != vEnd; vIt1++){
        if(T[*vIt1].bag == bag1){
            b1 = *vIt1;
            bag1_found = true;
            break;
        }
    }
    for(boost::tie(vIt2, vEnd) = boost::vertices(T); vIt2 != vEnd; vIt2++){
        if(T[*vIt2].bag == bag2){
            b2 = *vIt2;
            bag2_found = true;
            break;
        }
    }

    if(!bag1_found){
        b1 = boost::add_vertex(T);
        T[b1].bag = bag1;
    }

    if(!bag2_found){
        b2 = boost::add_vertex(T);
        T[b2].bag = bag2;
    }

    if(!(boost::edge(b1, b2, T).second || boost::edge(b2, b1, T).second))
        boost::add_edge(b1, b2, T);
}


template <typename G_t, typename T_t>
void exact_cutset(G_t &G, T_t &T, int lb){
    if(boost::num_vertices(G) == 0)
        return;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);

    if(boost::num_vertices(G) == 1){
        std::set<unsigned int> bag;
        bag.insert(G[*vIt].id);
        typename boost::graph_traits<T_t>::vertex_descriptor t = boost::add_vertex(T);
        T[t].bag = bag;
        return;
    }
        

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> cut;
    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component;
    

    cut.insert(*vIt);
    
    vIt++;
    for(; vIt != vEnd; vIt++)
        component.insert(*vIt);

    unsigned int k = (lb >= 0) ? (unsigned int)lb : 0;

    typename std::vector<typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > results;
    
    while(!explore_cutsets(G, cut, component, results, k++)){
        results.clear();
    }

    for(unsigned int i = 0; i < results.size()-1; i++){
        std::set<unsigned int> bag1;
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = results[i].begin(); sIt != results[i].end(); sIt++)
            bag1.insert(G[*sIt].id);

        std::set<unsigned int> bag2;
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = results[i+1].begin(); sIt != results[i+1].end(); sIt++)
            bag2.insert(G[*sIt].id);           

        glue_bags(T, bag1, bag2);
        i++;
    }
}

template <typename G_t, typename T_t>
void exact_cutset(G_t &G, T_t &T){
    int lb = -1;
    exact_cutset(G, T, lb);
}

}
