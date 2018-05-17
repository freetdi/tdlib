

namespace obsolete{

/* The top-down computation on tree decompositions is equal for some problems.
 *
 *   - 'cur' denotes the current node in T.
 *   - 'results' stores the computed table for all nodes of T.
 *   - 'val' denotes the value of a set that has to be choosen from tables of
 *       nodes in the subtree of 'T' with root 'cur'.
 *   - 'S' denotes the so far choosen set.
 *   - 'S_comp' denotes the complement of 'S' with respect to V(G).
 *   - according to (nice treedecomposition-) node types, the choice in the
 *       table of the child of 'cur' must be restricted:
 *         take_flag = 0 -> no restriction
 *         take_flag = 1 -> 'subset' must be choosen
 *         take_flag = 2 -> 'subset' must be a subset
 */

template <typename T_t>
void top_down_computation(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    typename std::vector<typename std::map<typename treedec_traits<T_t>::bag_type, int> > &results,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &S,
                    typename treedec_traits<T_t>::bag_type &S_comp,
                    typename treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                results[cur].begin(); it != results[cur].end(); it++)
        {
            if(it->second == (int)val){
                S.insert(it->first.begin(), it->first.end());
                return;
            }
        }
    }
    else if(node_type == treedec::nice::INTRODUCE || node_type == treedec::nice::FORGET){
        if(take_flag == 1){
            val = results[cur][subset];
        }

        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                   results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val 
               && std::includes(it->first.begin(), it->first.end(),
                                subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
                                        it->first.begin(), it->first.end(),
                                        std::inserter(S_comp, S_comp.begin()));
                    S.insert(it->first.begin(), it->first.end());
                    subset = it->first;
                    break;
                }
            }
        }

        if(node_type == treedec::nice::INTRODUCE){
            if(S.find(treedec::nice::get_introduced_vertex(cur, T)) != S.end()){
                val = val - 1;
            }
        }

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        if(node_type == treedec::nice::FORGET){
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex(cur, T));
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                          results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val
              && std::includes(it->first.begin(), it->first.end(),
                               subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    typename treedec_traits<T_t>::bag_type must_take = it->first;
                    S.insert(it->first.begin(), it->first.end());

                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
                                        it->first.begin(), it->first.end(),
                                        std::inserter(S_comp, S_comp.begin()));

                    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur, T); nIt != nEnd; nIt++){
                        top_down_computation(T, *nIt, results, results[*nIt][must_take], S, S_comp, must_take, 1);
                    }
                    return;
                }
            }
        }
    }
}

} //namespace obsolete
