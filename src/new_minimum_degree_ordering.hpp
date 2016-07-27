#include <boost/graph/minimum_degree_ordering.hpp>

namespace treedec{

namespace detail{

    template<class Graph, class DegreeMap, 
             class InversePermutationMap, 
             class PermutationMap,
             class SuperNodeMap, 
             class VertexIndexMap>
    class mmd_impl
    {
      // Typedefs
      typedef boost::graph_traits<Graph> Traits;
      typedef typename Traits::vertices_size_type size_type;
      typedef typename boost::detail::integer_traits<size_type>::difference_type 
        diff_t;
      typedef typename Traits::vertex_descriptor vertex_t;
      typedef typename Traits::adjacency_iterator adj_iter;
      typedef boost::iterator_property_map<vertex_t*, 
        boost::identity_property_map, vertex_t, vertex_t&> IndexVertexMap;
      typedef boost::detail::Stacks<diff_t> Workspace;
      typedef boost::bucket_sorter<size_type, vertex_t, DegreeMap, VertexIndexMap> 
        DegreeLists;
      typedef boost::detail::Numbering<InversePermutationMap, diff_t, vertex_t,VertexIndexMap>
        NumberingD;
      typedef boost::detail::degreelists_marker<diff_t, vertex_t, VertexIndexMap> 
        DegreeListsMarker;
      typedef boost::detail::Marker<diff_t, vertex_t, VertexIndexMap> MarkerP;

      // Data Members

      // input parameters
      Graph& G;
      DegreeMap degree;
      InversePermutationMap inverse_perm;
      PermutationMap perm;
      SuperNodeMap supernode_size;
      VertexIndexMap vertex_index_map;

      // internal data-structures
      std::vector<vertex_t> index_vertex_vec;
      size_type n;
      IndexVertexMap index_vertex_map;
      DegreeLists degreelists;
      NumberingD numbering;
      DegreeListsMarker degree_lists_marker;
      MarkerP marker;
      Workspace work_space;
      unsigned ub;
      unsigned count1, count2, count3, count4;
      std::vector<bool> eliminated;
    public:
      unsigned width;

      mmd_impl(Graph& g, size_type n_, DegreeMap degree, 
               InversePermutationMap inverse_perm, 
               PermutationMap perm,
               SuperNodeMap supernode_size, 
               VertexIndexMap id,
               unsigned ub_) 
        : G(g), degree(degree), 
        inverse_perm(inverse_perm), 
        perm(perm), 
        supernode_size(supernode_size), 
        vertex_index_map(id),
        index_vertex_vec(n_), 
        n(n_),
        degreelists(n_ + 1, n_, degree, id),
        numbering(inverse_perm, n_, vertex_index_map),
        degree_lists_marker(n_, vertex_index_map), 
        marker(n_, vertex_index_map),
        work_space(n_),
        ub(ub_),
        width(0)
      {
        typename boost::graph_traits<Graph>::vertex_iterator v, vend;
        size_type vid = 0;
        for (boost::tie(v, vend) = vertices(G); v != vend; ++v, ++vid)
          index_vertex_vec[vid] = *v;
        index_vertex_map = IndexVertexMap(&index_vertex_vec[0]);

        // Initialize degreelists.  Degreelists organizes the nodes
        // according to their degree.
        for (boost::tie(v, vend) = vertices(G); v != vend; ++v) {
          put(degree, *v, out_degree(*v, G));
          degreelists.push(*v);
        }
      }

      void do_mmd()
      {
        // Eliminate the isolated nodes -- these are simply the nodes
        // with no neighbors, which are accessible as a list (really, a
        // stack) at location 0.  Since these don't affect any other
        // nodes, we can eliminate them without doing degree updates.
        typename DegreeLists::stack list_isolated = degreelists[0];
        while (!list_isolated.empty()) {
          vertex_t node = list_isolated.top();
          marker.mark_done(node);
          numbering(node); count3++;
          numbering.increment();
          list_isolated.pop();
        }
        size_type min_degree = 1;
        typename DegreeLists::stack list_min_degree = degreelists[min_degree];

        while (list_min_degree.empty()) {
          ++min_degree;
          list_min_degree = degreelists[min_degree];
        }

        // check if the whole eliminating process is done
        while (!numbering.all_done()) {

          typename Workspace::stack llist = work_space.make_stack();

            // Find the next non-empty degree
            for (list_min_degree = degreelists[min_degree]; 
                 list_min_degree.empty(); 
                 ++min_degree, list_min_degree = degreelists[min_degree])
              ;

            const vertex_t node = list_min_degree.top();
            const size_type node_id = get(vertex_index_map, node);

            unsigned actual_deg = min_degree+supernode_size[node]-1;

            if (actual_deg > ub)
              throw exception_unsuccessful();

            width = (actual_deg > width)? actual_deg : width;


            list_min_degree.pop();
            numbering(node); count1++;

            // check if node is the last one
            if (numbering.all_done(supernode_size[node])) {
              numbering.increment(supernode_size[node]);
              break;
            }
            marker.increment_tag();
            marker.mark_tagged(node);

/*
            std::cout << "eliminate: " << G[node].id << "(deg:" << min_degree << ")" << std::endl;
            std::cout << "supernodeinfo: " << supernode_size[node] << std::endl;
*/

            this->eliminate(node);

            numbering.increment(supernode_size[node]);
            llist.push(node_id);

          if (numbering.all_done()) 
            break;

          this->update( llist, min_degree);
        } 
      } // do_mmd()

      void eliminate(vertex_t node)
      {
        typename Workspace::stack element_neighbor = work_space.make_stack();

        // Create two function objects for edge removal
        typedef typename Workspace::stack WorkStack;
        boost::detail::predicateRemoveEdge1<Graph, MarkerP, NumberingD, 
                             WorkStack, VertexIndexMap>
          p(G, marker, numbering, element_neighbor, vertex_index_map);

        boost::detail::predicate_remove_tagged_edges<Graph, MarkerP> p2(G, marker);

        // Reconstruct the adjacent node list, push element neighbor in a List.
        boost::remove_out_edge_if(node, p, G);
        //during removal element neighbors are collected.

        while (!element_neighbor.empty()) {
          // element absorb
          size_type e_id = element_neighbor.top();
          vertex_t element = get(index_vertex_map, e_id);
          adj_iter i, i_end;
          for (boost::tie(i, i_end) = adjacent_vertices(element, G); i != i_end; ++i){
            vertex_t i_node = *i;
            if (!marker.is_tagged(i_node) && !numbering.is_numbered(i_node)) {
              marker.mark_tagged(i_node);
              add_edge(node, i_node, G);
            }
          }
          element_neighbor.pop();
        }
        adj_iter v, ve;
        for (boost::tie(v, ve) = boost::adjacent_vertices(node, G); v != ve; ++v) {
          vertex_t v_node = *v;
          if (!degree_lists_marker.need_update(v_node) 
              && !degree_lists_marker.outmatched_or_done(v_node)) {
            degreelists.remove(v_node);
          }
          //update out edges of v_node
          boost::remove_out_edge_if(v_node, p2, G);

          if ( boost::out_degree(v_node, G) == 0 ) { // indistinguishable nodes
            supernode_size[node] += supernode_size[v_node];
            supernode_size[v_node] = 0;
            numbering.indistinguishable(v_node, node);
            marker.mark_done(v_node);
            degree_lists_marker.mark(v_node);
          } else {                           // not indistinguishable nodes
            add_edge(v_node, node, G);
            degree_lists_marker.mark_need_update(v_node);
          }
        }
      } // eliminate()


      template <class Stack>
      void update(Stack llist, size_type& min_degree)
      {
        size_type min_degree0 = min_degree + 1;

        while (! llist.empty()) {
          size_type deg, deg0 = 0;

          marker.set_multiple_tag(min_degree0);
          typename Workspace::stack q2list = work_space.make_stack();
          typename Workspace::stack qxlist = work_space.make_stack();

          vertex_t current = boost::get(index_vertex_map, llist.top());
          adj_iter i, ie;
          for (boost::tie(i,ie) = boost::adjacent_vertices(current, G); i != ie; ++i) {
            vertex_t i_node = *i;
            const size_type i_id =boost:: get(vertex_index_map, i_node);
            if (supernode_size[i_node] != 0) {
              deg0 += supernode_size[i_node];
              marker.mark_multiple_tagged(i_node);
              if (degree_lists_marker.need_update(i_node)) {
                if (out_degree(i_node, G) == 2)
                  q2list.push(i_id);
                else
                  qxlist.push(i_id);
              }
            }
          }

          while (!q2list.empty()) {
            const size_type u_id = q2list.top();
            vertex_t u_node = boost::get(index_vertex_map, u_id);
            // if u_id is outmatched by others, no need to update degree
            if (degree_lists_marker.outmatched_or_done(u_node)) {
              q2list.pop();
              continue;
            }
            marker.increment_tag();
            deg = deg0;

            adj_iter nu = boost::adjacent_vertices(u_node, G).first;
            vertex_t neighbor = *nu;
            if (neighbor == u_node) {
              ++nu;
              neighbor = *nu;
            }
            if (numbering.is_numbered(neighbor)) {
              adj_iter i, ie;
              for (boost::tie(i,ie) = boost::adjacent_vertices(neighbor, G);
                   i != ie; ++i) {
                const vertex_t i_node = *i;
                if (i_node == u_node || supernode_size[i_node] == 0)
                  continue;
                if (marker.is_tagged(i_node)) {
                  if (degree_lists_marker.need_update(i_node)) {
                    if ( boost::out_degree(i_node, G) == 2 ) { // is indistinguishable
                      supernode_size[u_node] += supernode_size[i_node];
                      supernode_size[i_node] = 0;
                      numbering.indistinguishable(i_node, u_node);
                      marker.mark_done(i_node);
                      degree_lists_marker.mark(i_node);
                    } else                        // is outmatched
                      degree_lists_marker.mark(i_node);
                  }
                } else {
                  marker.mark_tagged(i_node);
                  deg += supernode_size[i_node];
                }
              }
            } else
              deg += supernode_size[neighbor];

            deg -= supernode_size[u_node];
            degree[u_node] = deg; //update degree
            degreelists[deg].push(u_node);
            //u_id has been pushed back into degreelists
            degree_lists_marker.unmark(u_node);
            if (min_degree > deg) 
              min_degree = deg;
            q2list.pop();
          } // while (!q2list.empty())

          while (!qxlist.empty()) {
            const size_type u_id = qxlist.top();
            const vertex_t u_node = get(index_vertex_map, u_id);

            // if u_id is outmatched by others, no need to update degree
            if (degree_lists_marker.outmatched_or_done(u_node)) {
              qxlist.pop();
              continue;
            }
            marker.increment_tag();
            deg = deg0;
            adj_iter i, ie;
            for (boost::tie(i, ie) = boost::adjacent_vertices(u_node, G); i != ie; ++i) {
              vertex_t i_node = *i;
              if (marker.is_tagged(i_node)) 
                continue;
              marker.mark_tagged(i_node);

              if (numbering.is_numbered(i_node)) {
                adj_iter j, je;
                for (boost::tie(j, je) = boost::adjacent_vertices(i_node, G); j != je; ++j) {
                  const vertex_t j_node = *j;
                  if (marker.is_not_tagged(j_node)) {
                    marker.mark_tagged(j_node);
                    deg += supernode_size[j_node];
                  }
                }
              } else
                deg += supernode_size[i_node];
            } // for adjacent vertices of u_node
            deg -= supernode_size[u_node];
            degree[u_node] = deg;
            degreelists[deg].push(u_node);
            // u_id has been pushed back into degreelists
            degree_lists_marker.unmark(u_node); 
            if (min_degree > deg)
              min_degree = deg;
            qxlist.pop();
          } // while (!qxlist.empty()) {

          marker.set_tag_as_multiple_tag();
          llist.pop();
        } // while (! llist.empty())

      } // update()

      void build_permutation(InversePermutationMap next,
                             PermutationMap prev) 
      {
        // collect the permutation info
        size_type i;
        for (i = 0; i < n; ++i) {
          diff_t size = supernode_size[get(index_vertex_map, i)];
          if ( size <= 0 ) {
            prev[i] = next[i];
            supernode_size[get(index_vertex_map, i)]
              = next[i] + 1;  // record the supernode info
          } else
            prev[i] = - next[i];
        }
        for (i = 1; i < n + 1; ++i) {
          if ( prev[i-1] > 0 )
            continue;
          diff_t parent = i;
          while ( prev[parent - 1] < 0 ) {
            parent = - prev[parent - 1];
          }

          diff_t root = parent;
          diff_t num = prev[root - 1] + 1;
          next[i-1] = - num;
          prev[root-1] = num;

          parent = i;
          diff_t next_node = - prev[parent - 1];
          while (next_node > 0) {
            prev[parent-1] = - root;
            parent = next_node;
            next_node = - prev[parent - 1];
          }
        }
        for (i = 0; i < n; i++) {
          diff_t num = - next[i] - 1;
          next[i] = num;
          prev[num] = i;
        }
      } // build_permutation()

    };

  } //namespace detail


  // MMD algorithm
  // 
  //The implementation presently includes the enhancements for mass
  //elimination, incomplete degree update, multiple elimination, and
  //external degree.
  //
  //Important Note: This implementation requires the BGL graph to be
  //directed.  Therefore, nonzero entry (i, j) in a symmetrical matrix
  //A coresponds to two directed edges (i->j and j->i).
  //
  //see Alan George and Joseph W. H. Liu, The Evolution of the Minimum
  //Degree Ordering Algorithm, SIAM Review, 31, 1989, Page 1-19
  template<class Graph, class DegreeMap, 
           class InversePermutationMap, 
           class PermutationMap,
           class SuperNodeMap, class VertexIndexMap>
  unsigned minimum_degree_ordering
    (Graph& G, 
     DegreeMap degree, 
     InversePermutationMap inverse_perm, 
     PermutationMap perm, 
     SuperNodeMap supernode_size, 
     VertexIndexMap vertex_index_map,
     unsigned ub = UINT_MAX)
  {
    detail::mmd_impl<Graph,DegreeMap,InversePermutationMap,
      PermutationMap, SuperNodeMap, VertexIndexMap> 
      impl(G, num_vertices(G), degree, inverse_perm, 
           perm, supernode_size, vertex_index_map, ub);
    impl.do_mmd();
    impl.build_permutation(inverse_perm, perm); //not necessary for an ub.

    return impl.width;
  }

} //namespace treedec
