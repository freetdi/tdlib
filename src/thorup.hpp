// Philipp Klaus Krause, philipp@informatik.uni-frankfurt.de, pkk@spth.de, 2010 - 2011
//
// some changes:
// Felix Salfelder 2017
//
// (c) 2010-2011 Goethe-Universität Frankfurt
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
// Some routines for tree-decompositions.
//
// A tree decomposition is a graph that has a set of vertex indices as bundled
// property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, tree_dec_node> tree_dec_t;
//
// The following are the routines that are most likely to be interesting for outside use:
//
// void nicify(T_t &T)
// Transforms a tree decomposition T into a nice tree decomposition
//
// void thorup_tree_decomposition(T_t &tree_decomposition, const G_t &cfg)
// Creates a tree decomposition T from a graph cfg using Thorup's heuristic.
//
// void tree_decomposition_from_elimination_ordering(
//           T_t &T, std::list<unsigned int>& l, const G_t &G)
// Creates a tree decomposition T of a graph G from an elimination ordering l.
//
// void thorup_elimination_ordering(l_t &l, const J_t &J)
// Creates an elimination ordering l of a graph J using Thorup's heuristic.
//
// Changelog:
//
// 2011 - initial
// 2017 - merged into tdlib
//        - namespace
//        - bag access
//        - list->vector
//
#ifndef TDLIB_PKK_THORUP_H
#define TDLIB_PKK_THORUP_H

#include <map>
#include <vector>
#include <set>
#include <stack>
#include <list>

#include <boost/tuple/tuple_io.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <tdlib/treedec.hpp> // get_bagsize

// namespace boost {
// // some hacks. perhaps wrong place
// using HACK_Bag=std::vector<unsigned>;
// using HACK_Treedec=adjacency_list<setS, vecS, undirectedS, HACK_Bag,
// 		  no_property, no_property, listS>;
// 
// HACK_Bag& get(treedec::bag_t, HACK_Treedec& t, unsigned u)
// { untested();
// 	return t[u];
// }
// 
// }

namespace treedec{

namespace pkk_thorup{

struct forget_properties
{
  template<class T1, class T2>
  void operator()(const T1&, const T2&) const
  {
  }
};

// Thorup algorithm D.
// The use of the multimap makes the complexity of this O(|I|log|I|), which could be reduced to O(|I|).
template <class l_t>
void thorup_D(l_t &l, const std::multimap<unsigned, unsigned> &MJ,
		const std::multimap<unsigned, unsigned> &MS,
		const unsigned n)
{
  std::map<unsigned int, unsigned int> m;

  l.clear();

  unsigned int i = 0;
  for (unsigned int j = n; j > 0;)
    {
      j--;
      if (m.find(j) == m.end())
        m[j] = i++;
        
      std::multimap<unsigned int, unsigned int>::const_iterator k, k_end;

      for (boost::tie(k, k_end) = MS.equal_range(j); k != k_end; ++k)
        if (m.find(k->second) == m.end())
          m[k->second] = i++;

      for (boost::tie(k, k_end) = MJ.equal_range(j); k != k_end; ++k)
        if (m.find(k->second) == m.end())
          m[k->second] = i++;
    }

  std::vector<unsigned int> v(n);

  std::map<unsigned int, unsigned int>::iterator mi;

  for (mi = m.begin(); mi != m.end(); ++mi){
    v[mi->second] = mi->first;
  }

  for (i = 0; i < n; i++){
    l.push_back(v[i]);
  }
}

// Thorup algorithm E.
// The use of the multimap makes the complexity of this O(|I|log|I|), which
// could be reduced to O(|I|).
template <class I_t>
void thorup_E(std::multimap<unsigned, unsigned> &M, I_t const& I)
{
  typedef typename boost::graph_traits<I_t>::adjacency_iterator adjacency_iter_t;
  typedef typename boost::property_map<I_t, boost::vertex_index_t>::type index_map;
  index_map index = boost::get(boost::vertex_index, I);

  std::stack<std::pair<int, unsigned int> > s;

  M.clear();

  s.push(std::pair<int, unsigned int>(-1, boost::num_vertices(I)));

  for (unsigned int i = 0; i < boost::num_vertices(I); i++)
    {
      unsigned int j = i;
      adjacency_iter_t j_curr, j_end;

      for (boost::tie(j_curr, j_end)=boost::adjacent_vertices(i, I);
			  	j_curr != j_end; ++j_curr){
        if (index[*j_curr] > j){
          j = index[*j_curr];
		  }else{
		  }
		}

      if (j == i){
        continue;
		}else{
		}

      while (s.top().second <= i) {
			M.insert(std::pair<unsigned int, unsigned int>(s.top().second, s.top().first));
			s.pop();
		}

      unsigned int i2 = i;
		while (j >= s.top().second && s.top().second > i2) {
			i2 = s.top().first;
			s.pop();
		}

      s.push(std::pair<int, unsigned int>(i2, j));
    }
    
    // Thorup forgot this in his paper. Without it, some maximal chains are omitted.
    while(s.size() > 1) {
        M.insert(std::pair<unsigned int, unsigned int>(s.top().second, s.top().first));
        s.pop();
    }
}

// Heuristically give an elimination ordering for a directed graph.
// For a description of this, including algorithms D and E, see
// Mikkel Thorup, "All Structured Programs have Small Tree-Width and Good
// Register Allocation", Appendix A.  The use of the multimap makes the
// complexity of this O(|I|log|I|), could be reduced to O(|I|).
//
// it seems that G is expected to be oriented. otherwise result is worse
template <class l_t, class G_t>
void thorup_elimination_ordering(l_t &l, const G_t &G)
{
  // Remove edges to immediately following instruction. By "each statement can
  // have at most one jump" in the last paragraph of Appendix A it is clear
  // that Thorup does not consider the implicit next-instruction-edges as
  // jumps.
  boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> J;
  boost::copy_graph(G, J, boost::vertex_copy(forget_properties()).edge_copy(forget_properties()));
  auto N=boost::num_vertices(J);
  for (unsigned int i = 1; i < N; i++){
    remove_edge(i-1, i, J);
  }

  // Todo: Implement a graph adaptor for boost that allows to treat directed
  // graphs as undirected graphs.
  //
  // TODO: is this copy necessary, if G is symmetric?
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> S;
  boost::copy_graph(J, S);

  std::multimap<unsigned int, unsigned int> MJ, MS;

  thorup_E(MJ, J);

  thorup_E(MS, S);

  thorup_D(l, MJ, MS, num_vertices(J));
}

// Finds a (the newest) bag that contains all vertices in X in the tree decomposition T.
template <class T_t>
typename boost::graph_traits<T_t>::vertex_iterator
find_bag(const std::set<unsigned int> &X, const T_t &T)
{
  typedef typename boost::graph_traits<T_t>::vertex_iterator T_vertex_iter_t;
  typedef typename std::set<unsigned int>::const_iterator vertex_index_iter_t;

  T_vertex_iter_t t, t_end, t_found;
  vertex_index_iter_t v;

  for (boost::tie(t, t_end) = vertices(T), t_found = t_end; t != t_end; ++t) {
	  for (v = X.begin(); v != X.end(); ++v){
		  auto const& b=boost::get(treedec::bag_t(), T, *t);
		  if (b.find(*v) == b.end()){ itested();
			  break;
		  }else{ itested();
		  }
	  }

	  if (v == X.end()){
		  t_found = t;
	  }
  }

  if (t_found == t_end) // Todo: Better error handling (throw exception?)
    {
      std::cerr << "find_bag() failed.\n";
      std::cerr.flush();
    }

  return(t_found);
}

// Add edges to make the vertices in X a clique in G.
template <class G_t>
void make_clique(const std::set<unsigned int> &X , G_t &G)
{
  std::set<unsigned int>::const_iterator n1, n2;
  for (n1 = X.begin(); n1 != X.end(); n1++){
    for (n2 = n1, ++n2; n2 != X.end(); ++n2){
      if(!edge(*n1, *n2, G).second){
			add_edge(*n1, *n2, G);
		}else{
		}
	 }
  }
}

template <class T_t, class v_t, class G_t>
void add_vertices_to_tree_decomposition(T_t &T, const v_t v, const v_t v_end, G_t &G, std::vector<bool> &active)
{
  // Base case: Empty graph. Create an empty bag.
  if (v == v_end)
    {
      boost::add_vertex(T);
      return;
    }

  // Todo: A more elegant solution, e.g. using subgraphs or filtered graphs.

  typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iter_t;
  typedef typename boost::property_map<G_t, boost::vertex_index_t>::type index_map;
  index_map index = boost::get(boost::vertex_index, G);

  // Get the neigbours
  adjacency_iter_t n, n_end;
  std::set<unsigned int> neighbours;
  for (boost::tie(n, n_end) = boost::adjacent_vertices(*v, G); n != n_end; ++n)
    if (active[index[*n]])
      neighbours.insert(index[*n]);

  // Recurse
  active[*v] = false;
  make_clique(neighbours, G);
  v_t v_next = v;
  add_vertices_to_tree_decomposition(T, ++v_next, v_end, G, active);

  // Add new bag
  typename boost::graph_traits<T_t>::vertex_iterator t;
  typename boost::graph_traits<T_t>::vertex_descriptor s;
  t = find_bag(neighbours, T);
  s = boost::add_vertex(T);
  boost::add_edge(*t, s, T);
  // T[s].bag = neighbours;
  // T[s].bag.insert(*v);
  auto& b=boost::get(treedec::bag_t(), T, s);
  b = neighbours;
  treedec::insert(b, *v);

}

// Create a tree decomposition from an elimination ordering.
// FIXME. not here?
template <class T_t, class G_t>
void tree_decomposition_from_elimination_ordering(T_t &T, const std::vector<unsigned int>& l, const G_t &G)
{
  auto v = l.rbegin();
  auto v_end = l.rend();

  // Todo: Implement a graph adaptor for boost that allows to treat directed
  // graphs as undirected graphs.
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G_sym;
  boost::copy_graph(G, G_sym, boost::vertex_copy(forget_properties()).edge_copy(forget_properties()));

  std::vector<bool> active(boost::num_vertices(G), true);

  add_vertices_to_tree_decomposition(T, v, v_end, G_sym, active);
}

template <class T_t, class G_t>
void thorup_tree_decomposition(T_t &tree_decomposition, const G_t &cfg)
{
	std::list<unsigned int> elimination_ordering;

	thorup_elimination_ordering(elimination_ordering, cfg);

	tree_decomposition_from_elimination_ordering(tree_decomposition, elimination_ordering, cfg);
}

// Ensure that all joins are at proper join nodes: Each node that has two children has the same bag as its children.
// Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_joins(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t)
{
  typedef typename boost::graph_traits<T_t>::adjacency_iterator adjacency_iter_t;

  adjacency_iter_t c, c_end;
  typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

  boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

  switch (out_degree(t, T))
    {
    case 0:
      return;
    case 1:
      nicify_joins(T, *c);
      return;
    case 2:
      break;
    default:
      c0 = *c++;
      c1 = *c;
      typename boost::graph_traits<T_t>::vertex_descriptor d;
      d = boost::add_vertex(T);
      add_edge(d, c0, T);
      add_edge(d, c1, T);
      boost::remove_edge(t, c0, T);
      boost::remove_edge(t, c1, T);
      T[d].bag = T[t].bag;
      boost::add_edge(t, d, T);
      nicify_joins(T, t);
      return;
    }

  c0 = *c++;
  c1 = *c;
  nicify_joins(T, c0);
  if (T[t].bag != T[c0].bag)
    {
      typename boost::graph_traits<T_t>::vertex_descriptor d;
      d = boost::add_vertex(T);
      boost::add_edge(d, c0, T);
      boost::remove_edge(t, c0, T);
      T[d].bag = T[t].bag;
      boost::add_edge(t, d, T);
    }
  nicify_joins(T, c1);
  if (T[t].bag != T[c1].bag)
    {
      typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
      boost::add_edge(d, c1, T);
      boost::remove_edge(t, c1, T);
      T[d].bag = T[t].bag;
      boost::add_edge(t, d, T);
    }
}

// Ensure that all nodes' bags are either a subset or a superset of their successors'.
// Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_diffs(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t)
{
  typedef typename boost::graph_traits<T_t>::adjacency_iterator adjacency_iter_t;

  adjacency_iter_t c, c_end;
  typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

  boost::tie(c, c_end) = adjacent_vertices(t, T);

  switch (boost::out_degree(t, T))
    {
    case 0:
      if (T[t].bag.size())
        boost::add_edge(t, boost::add_vertex(T), T);
      return;
    case 1:
      break;
    case 2:
      c0 = *c++;
      c1 = *c;
      nicify_diffs(T, c0);
      nicify_diffs(T, c1);
      return;
    default:
      std::cerr << "nicify_diffs error.\n";
      return;
    }

  c0 = *c;
  nicify_diffs(T, c0);

  if (std::includes(T[t].bag.begin(), T[t].bag.end(), T[c0].bag.begin(), T[c0].bag.end()) || std::includes(T[c0].bag.begin(), T[c0].bag.end(), T[t].bag.begin(), T[t].bag.end()))
    return;

  typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);

  boost::add_edge(d, c0, T);
  boost::remove_edge(t, c0, T);
  std::set_intersection(T[t].bag.begin(), T[t].bag.end(), T[c0].bag.begin(), T[c0].bag.end(), std::inserter(T[d].bag, T[d].bag.begin()));
  boost::add_edge(t, d, T);
}

// // Ensure that all nodes' bags' sizes differ by at most one to their successors'.
template <class T_t>
void nicify_diffs_more(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t)
{
  typedef typename boost::graph_traits<T_t>::adjacency_iterator adjacency_iter_t;

  adjacency_iter_t c, c_end;
  typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

  boost::tie(c, c_end) = adjacent_vertices(t, T);

  switch (boost::out_degree(t, T))
    {
    case 0:
      if (T[t].bag.size() > 1)
        {
          typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
          T[d].bag = T[t].bag;
          T[d].bag.erase(T[d].bag.begin());
          T[d].weight = 0;
          boost::add_edge(t, d, T);
          nicify_diffs_more(T, t);
        }
      else
        T[t].weight = 0;
      return;
    case 1:
      break;
    case 2:
      c0 = *c++;
      c1 = *c;
      nicify_diffs_more(T, c0);
      nicify_diffs_more(T, c1);
      T[t].weight = std::min(T[c0].weight, T[c1].weight) + 1;
      return;
    default:
      std::cerr << "nicify_diffs_more error.\n";
      return;
    }

  c0 = *c;

  size_t c0_size, t_size;
  t_size = T[t].bag.size();
  c0_size = T[c0].bag.size();

  if (t_size <= c0_size + 1 && t_size + 1 >= c0_size)
    {
      nicify_diffs_more(T, c0);
      T[t].weight = T[c0].weight;
      return;
    }

  typename boost::graph_traits<T_t>::vertex_descriptor d = add_vertex(T);
  boost::add_edge(d, c0, T);
  boost::remove_edge(t, c0, T);
  T[d].bag = T[t_size > c0_size ? t : c0].bag;
  std::set<unsigned int>::iterator i;

  // YIKES.
  for (i = T[d].bag.begin(); T[t_size < c0_size ? t : c0].bag.find(*i) != T[t_size < c0_size ? t : c0].bag.end(); ++i);

  T[d].bag.erase(i);
  boost::add_edge(t, d, T);

  nicify_diffs_more(T, t);
}

// Find a root of an acyclic graph T
// Complexity: Linear in the number of vertices of T.
template <class T_t>
typename boost::graph_traits<T_t>::vertex_descriptor find_root(T_t &T)
{
  typename boost::graph_traits<T_t>::vertex_descriptor t;
  typename boost::graph_traits<T_t>::in_edge_iterator e, e_end;

  t = *(boost::vertices(T).first);

  for (boost::tie(e, e_end) = boost::in_edges(t, T);
		  e!=e_end; boost::tie(e, e_end) = boost::in_edges(t, T)){
    t = boost::source(*e, T);
  }

  return(t);
}

// Transform a tree decomposition into a nice tree decomposition.
template <class T_t>
void nicify(T_t &T) {
  typename boost::graph_traits<T_t>::vertex_descriptor t;

  t = find_root(T);

  // Ensure we have an empty bag at the root.
  if(T[t].bag.size()) {
    typename boost::graph_traits<T_t>::vertex_descriptor d = t;
    t = add_vertex(T);
    boost::add_edge(t, d, T);
  }else{
  }

  nicify_joins(T, t);
  nicify_diffs(T, t);
  nicify_diffs_more(T, t);
}

} // pkk_thorup

template<typename G, template<class G_, class ...> class config=algo::default_config>
class thorup : public ::treedec::algo::draft::algo1{
public: // types
	// FIXME: use treedec_t<G>
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
			  boost::property<treedec::bag_t, std::set<unsigned> > > decomp_type;
public:
	thorup(G const& g)
		: algo1("thorup"), _g(g){
	}

	void do_it(){
		pkk_thorup::thorup_elimination_ordering(_o, _g);

		// check: do we have better algorithms?
		pkk_thorup::tree_decomposition_from_elimination_ordering(
				_t, _o, _g);

#ifdef DEBUG
		G h;
		boost::copy_graph(_g, h);
		auto w=get_width_of_elimination_ordering(h, _o);
		std::cerr << w << " " << treedec::get_width(_t) << "\n";
		assert(w==treedec::get_width(_t));
#endif
	}

	template<class T>
	void get_tree_decomposition(T& t) const{
		boost::copy_graph(_t, t);
	}
	decomp_type const& get_tree_decomposition() const{
		return _t;
	}
	size_t get_bagsize() const{
		return treedec::get_bagsize(_t);
	}
	std::list<unsigned> const&  get_elimord() const{
		return _o;
	}
private:
	G const& _g;
	std::vector<unsigned> _o;
	decomp_type _t;
};

} // treedec

#endif
