// Holger Dell & co
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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <stack>
#include <string>
#include <unordered_set>
#include <vector>
#include "trace.hpp"

/*
 * An unsigned variant of std::stoi that also checks if the input consists
 * entirely of base-10 digits
 */
unsigned pure_stou(const std::string& s) {
  if(s.empty()) { untested();
    throw std::invalid_argument("Empty string");
  }
  unsigned long result = 0;
  unsigned long power = 1;
  for (int i = s.length()-1; i >= 0; i--) {
    char c = s[i];
    if (c < '0' || c > '9') { untested();
      throw std::invalid_argument("Non-numeric entry '" + s + "'");
    }
    result += power*(c-'0');
    power *= 10;
  }
  if (result > std::numeric_limits<unsigned>::max()) { untested();
    throw std::out_of_range("stou");
  }
  return result;
}

/*
 * The most basic graph structure you could imagine, in adjacency list
 * representation.  This is not supposed to be any kind of general-purpose graph
 * class, but only implements what we need for our computations. E.g., removing
 * vertices is not supported.
 */
struct graph {
  std::vector<std::vector<unsigned> > adj_list;
  unsigned num_vertices = 0;
  unsigned num_edges = 0;

  /*
   * Adds a vertex to the vertex set and returns its index.
   * Vertices in this class are indexed starting with 0 (!).
   * This is important because the vertices in the file format are not.
   */
  unsigned add_vertex() {
    adj_list.push_back(std::vector<unsigned>());
    return num_vertices++;
  }

  /*
   * The neighbors of a given vertex, where the index is, as usual, starting
   * with 0. Does *not* include the vertex itself
   */
  std::vector<unsigned>& neighbors(unsigned vertex_index) {
    return adj_list[vertex_index];
  }

  /*
   * Adds an undirected edge between two vertices, identified by their index; no
   * checks are performed and bad indices can cause segfaults.
   */
  void add_edge(unsigned u, unsigned v) {
    adj_list.at(u).push_back(v);
    adj_list.at(v).push_back(u);
    num_edges++;
  }

  /*
   * Removes an undirected edge
   */
  void remove_edge(unsigned u, unsigned v) {
    bool end1 = false;
    bool end2 = false;
    auto v_it = std::find(adj_list.at(u).begin(), adj_list.at(u).end(), v);
    if (v_it != adj_list.at(u).end()) {
      adj_list.at(u).erase(v_it);
      end1 = true;
    }

    auto u_it = std::find(adj_list.at(v).begin(), adj_list.at(v).end(), u);
    if (u_it != adj_list.at(v).end()) {
      adj_list.at(v).erase(u_it);
      end2 = true;
    }

    if (end1 && end2) {
      num_edges--;
    }
  }
};

/*
 * Type for a tree decomposition; i.e. a set of bags (vertices) together with
 * adjacency lists on this set.  As the graph class, this is highly minimalistic
 * and implements all operations in such a way that our algorithms may work on
 * it.
 */
struct tree_decomposition {
  typedef std::set<unsigned> vertex_t;
  std::vector<vertex_t> bags;
  std::vector<std::vector<unsigned>> adj_list;

  /*
   * The number of *relevant* bags currently in the tree; bags that were removed
   * using remove_vertex and continue to exist empty and isolated are not
   * counted.
   */
  unsigned num_vertices = 0;
  unsigned num_edges = 0;

  /*
   * Adds a given bag to the vertex set of the decomposition and returns its
   * index.  Vertices in this class are indexed starting with 0 (!).  This is
   * important because the vertices in the file format are not.
   */
  unsigned add_bag(vertex_t& bag) {
    bags.push_back(bag);
    adj_list.push_back(std::vector<unsigned>());
    return num_vertices++;
  }

  /*
   * See the graph class
   */
  std::vector<unsigned>& neighbors(unsigned vertex_index) {
    return adj_list.at(vertex_index);
  }

  /*
   * See the graph class
   */
  void add_edge(unsigned u, unsigned v) {
    adj_list.at(u).push_back(v);
    adj_list.at(v).push_back(u);
    num_edges++;
  }

  /*
   * See the graph class
   */
  void remove_edge(unsigned u, unsigned v) {
    auto v_it = std::find(adj_list.at(u).begin(),adj_list.at(u).end(),v);
    auto u_it = std::find(adj_list.at(v).begin(),adj_list.at(v).end(),u);
    if (v_it != adj_list.at(u).end() && u_it != adj_list.at(v).end()) {
      adj_list.at(u).erase(v_it);
      adj_list.at(v).erase(u_it);
      num_edges--;
    }
  }

  /*
   * Removes a vertex, in the following sense: The bag corresponding to the
   * index is emptied and all adjacencies are removed, i.e., the bag will
   * contain 0 vertices and have no incident edges. Nevertheless, the number of
   * vertices is reduced (hence num_vertices--).
   */
  void remove_vertex(unsigned u) {
    bags.at(u).clear();
    std::vector<unsigned> remove;
    for (auto it = adj_list.at(u).begin(); it != adj_list.at(u).end(); it++) {
      remove.push_back(*it);
    }
    for (auto it = remove.begin(); it != remove.end(); it++) {
      remove_edge(u,*it);
    }
    num_vertices--;
  }

  /*
   * Get the u-th bag
   */
  vertex_t& get_bag(unsigned u) {
    return bags.at(u);
  }

  /*
   * Checks if the given decomposition constitutes a tree using DFS
   */
  bool is_tree() {
    if ((num_vertices > 0) && num_vertices - 1 != num_edges) { untested();
      return false;
    } else if (num_vertices == 0) { untested();
      return (num_edges == 0);
    }

    std::vector<int> seen(num_vertices, 0);
    unsigned seen_size = 0;

    bool cycle = !tree_dfs(seen,0,-1,seen_size);
    if (cycle || seen_size != num_vertices) { untested();
      return false;
    }
    return true;
  }

  /*
   * Helper method for is_tree(); not to be called from outside of the class
   */
  bool tree_dfs(std::vector<int>& seen, unsigned root, unsigned parent, unsigned& num_seen) {
    if (seen[root] != 0) { untested();
      return false;
    }

    seen[root] = 1;
    num_seen++;

    for (auto it = adj_list[root].begin(); it != adj_list[root].end(); it++) {
      if (*it != parent) {
        if (!tree_dfs(seen, *it, root,num_seen)) { untested();
          return false;
        }
      }
    }
    return true;
  }
};

/*
 * The different states the syntax checker may find itself in while checking the
 * syntax of a *.gr- or *.td-file.
 */
enum State {
  COMMENT_SECTION,
  S_LINE,
  BAGS,
  EDGES,
  P_LINE
};

/*
 * Messages associated with the respective exceptions to be thrown while
 * checking the files
 */
enum td_error{
	E_OK=0,
	E_INV_FMT,
	E_INV_SOLN,
	E_INV_SOLN_BAGSIZE,
	E_INV_EDGE,
	E_INV_BAG,
	E_INV_BAG_INDEX,
	E_INC_SOLN,
	E_NO_BAG_INDEX,
	E_BAG_MISSING,
	E_FILE_ERROR,
	E_EMPTY_LINE,
	E_INV_PROB,
	E_NUM_ERRORS
};

const char* errstr[E_NUM_ERRORS] = { //
	"OK",
	"Invalid format",
	"Invalid s-line",
	"Invalid s-line: Reported bagsize and actual bagsize differ",
	"Invalid edge",
	"Invalid bag",
	"Invalid bag index",
	"Inconsistent values in s-line",
	"No bag index given",
	"Bag missing",
	"Could not open file",
	"No empty lines allowed",
	"Invalid p-line",
};

class exception_invalid_td : public std::exception{
public:
	~exception_invalid_td() throw() {}
   exception_invalid_td(td_error w)
      : std::exception(), _w(w) { untested(); }

   const char* what() const throw(){ itested();
      return errstr[_w];
   }

private:
   td_error _w;
};

/*
 * The state the syntax checker is currently in; this variable is used for both
 * checking the tree *and* the graph
 */
State current_state = COMMENT_SECTION;

/*
 * A collection of global variables that are relentlessly manipulated and read
 * from different positions of the program.  Initially, they are set while
 * reading the input files.
 */

/* The number of vertices of the graph underlying the decomposition, as stated
 * in the *.td-file
 */
unsigned n_graph;

/* The number of bags as stated in the *.td-file */
unsigned n_decomp;

/* The width of the decomposition as stated in the *.td-file */
unsigned width;

/* The number of vertices of the graph as stated in the *.gr-file */
unsigned n_vertices;

/* The number of edges of the graph as stated in the *.gr-file */
unsigned n_edges;

/* The maximal width of some bag, to compare with the one stated in the file */
unsigned real_width;

/* A vector to record which of the bags were seen so far, while going through
 * the *.td-file, as to ensure uniqueness of each bag.
 */
std::vector<int> bags_seen;

/* Temporary storage for all the bags before they are inserted into the actual
 * decomposition; we might as well directly manipulate the tree TODO
 */
std::vector<std::set<unsigned> > bags;

/*
 * Given the tokens from one line (split on whitespace), this reads the
 * s-line from these tokens and initializes the corresponding globals
 */
void read_solution(const std::vector<std::string>& tokens)
{
  if (current_state != COMMENT_SECTION) { untested();
    throw exception_invalid_td(E_INV_FMT);
  }
  current_state = S_LINE;
  if(tokens.size() != 5 || tokens[1] != "td") { untested();
    throw exception_invalid_td(E_INV_SOLN);
  }

  n_decomp = pure_stou(tokens[2]);
  width = pure_stou(tokens[3]);
  n_graph = pure_stou(tokens[4]);
  if (width > n_graph) { untested();
    throw exception_invalid_td(E_INC_SOLN);
  }
}

/*
 * Given the tokens from one line (split on whitespace), this reads the bag
 * represented by these tokens and manipulates the global bags accordingly
 */
void read_bag(const std::vector<std::string>& tokens)
{
  if (current_state == S_LINE) {
    current_state = BAGS;
    bags.resize(n_decomp);
    bags_seen.resize(n_decomp,0);
  }

  if (current_state != BAGS) { untested();
    throw exception_invalid_td(E_INV_FMT);
  }

  if(tokens.size() < 2) { untested();
    throw exception_invalid_td(E_NO_BAG_INDEX);
  }

  unsigned bag_num = pure_stou(tokens[1]);
  if (bag_num < 1 || bag_num > n_decomp || bags_seen[bag_num-1] != 0) { untested();
    throw exception_invalid_td(E_INV_BAG_INDEX);
  }
  bags_seen[bag_num-1] = 1;

  for(unsigned i = 2; i < tokens.size(); i++) {
    if (tokens[i] == "") break;
    unsigned id = pure_stou(tokens[i]);
    if (id < 1 || id > n_graph) { untested();
      throw exception_invalid_td(E_INV_BAG);
    }

    bags[bag_num-1].insert(id);
  }

  if(bags[bag_num-1].size() > real_width) {
    real_width = bags[bag_num-1].size();
  }
}

/*
 * Given the tokens from one line (split on whitespace) and a tree
 * decomposition, this reads the edge represented by this line (in the
 * decomposition) and adds the respective edge to the tree decomposition
 */
void read_decomp_edge(const std::vector<std::string>& tokens, tree_decomposition &T)
{
  if (current_state == BAGS) {
    for (auto it = bags_seen.begin(); it != bags_seen.end(); it++) {
      if (*it == 0) { untested();
        throw exception_invalid_td(E_BAG_MISSING);
      }
    }

    for (auto it = bags.begin(); it != bags.end(); it++) {
      T.add_bag(*it);
    }
    current_state = EDGES;
  }

  if (current_state != EDGES) { untested();
    throw exception_invalid_td(E_INV_FMT);
  }

  unsigned s = pure_stou(tokens[0]);
  unsigned d = pure_stou(tokens[1]);
  if(s < 1 || d < 1 || s > n_decomp || d > n_decomp) { untested();
    throw exception_invalid_td(E_INV_EDGE);
  }
  T.add_edge(s-1, d-1);
}

/*
 * Given the tokens from one line (split on whitespace), this reads the
 * p-line from these tokens and initializes the corresponding globals
 */
void read_problem(const std::vector<std::string>& tokens, graph& g) {
  if (current_state != COMMENT_SECTION) { untested();
    throw exception_invalid_td(E_INV_FMT);
  }
  current_state = P_LINE;

  if(tokens.size() != 4 || tokens[1] != "tw") { untested();
    throw exception_invalid_td(E_INV_PROB);
  }

  n_vertices = pure_stou(tokens[2]);
  n_edges = pure_stou(tokens[3]);

  while (g.add_vertex()+1 < n_vertices) {};
}

/*
 * Given the tokens from one line (split on whitespace) and a tree
 * decomposition, this reads the edge (in the graph) represented by this line
 * and adds the respective edge to the graph
 */
void read_graph_edge(const std::vector<std::string>& tokens, graph& g)
{
  if (current_state == P_LINE) {
    current_state = EDGES;
  }

  if (current_state != EDGES) { untested();
    throw exception_invalid_td(E_INV_FMT);
  }

  unsigned s = pure_stou(tokens[0]);
  unsigned d = pure_stou(tokens[1]);
  if(s < 1 || d < 1 || s > n_vertices || d > n_vertices) { untested();
    throw exception_invalid_td(E_INV_EDGE);
  }
  g.add_edge(s-1, d-1);
}

/*
 * Given a stream to the input file in the *.gr-format, this reads from the file
 * the graph represented by this file.  If the file is not conforming to the
 * format, it throws a corresponding exception_invalid_td with one of the error
 * messages defined above.
 */
void read_graph(std::ifstream& fin, graph& g) {
  current_state = COMMENT_SECTION;
  n_edges = -1;
  n_vertices = -1;

  if(!fin.is_open()){ untested();
    throw std::invalid_argument("cannot open file");
  }else{
  }

  std::string line;
  std::string delimiter = " ";

  while(std::getline(fin, line)) {
    if(line == "" || line == "\n") { untested();
      throw std::invalid_argument("empty line");
    }

    std::vector<std::string> tokens;
    tokens.reserve(line.length());
    size_t oldpos = 0;
    size_t newpos = 0;

    while(newpos != std::string::npos) {
      newpos = line.find(delimiter, oldpos);
      tokens.push_back(line.substr(oldpos, newpos-oldpos));
      oldpos = newpos + delimiter.size();
    }
    if (tokens[0] == "c") {
      continue;
    } else if (tokens[0] == "p") {
      read_problem(tokens,g);
    } else if (tokens.size() == 2){
      read_graph_edge(tokens, g);
    } else { untested();
      throw std::invalid_argument(std::string(errstr[E_INV_EDGE]) + " (an edge has exactly two endpoints)");
    }
  }

  if (g.num_edges != n_edges) { untested();
    throw std::invalid_argument(std::string(errstr[E_INV_PROB]) + " (incorrect number of edges)");
  }
}

/*
 * Given a stream to the input file in the *.td-format, this reads from the file
 * the decomposition represented by this file.  If the file is not conforming to
 * the format, it throws a corresponding std::invalid_argument with one of the
 * error messages defined above.
 */
void read_tree_decomposition(std::ifstream& fin, tree_decomposition& T)
{
  current_state = COMMENT_SECTION;
  n_graph = -1;
  n_decomp = 0;
  width = -2;
  bags_seen.clear();
  bags.clear();

  if(!fin.is_open()){ untested();
    throw std::invalid_argument("stream error");
  }

  std::string line;
  std::string delimiter = " ";

  while(std::getline(fin, line)) {
    if(line == "" || line == "\n") { untested();
      throw std::invalid_argument("empty line");
    }

    std::vector<std::string> tokens;
    tokens.reserve(line.length());
    size_t oldpos = 0;
    size_t newpos = 0;

    while(newpos != std::string::npos) {
      newpos = line.find(delimiter, oldpos);
      tokens.push_back(line.substr(oldpos, newpos-oldpos));
      oldpos = newpos + delimiter.size();
    }

    if (tokens[0] == "c") {
      continue;
    } else if (tokens[0] == "s") {
      read_solution(tokens);
    } else if (tokens[0] == "b") {
      read_bag(tokens);
    } else {
      read_decomp_edge(tokens, T);
    }
  }

  if (current_state == BAGS) { untested();
    for (auto it = bags.begin(); it != bags.end(); it++) { untested();
      T.add_bag(*it);
    }
  }

  if (width != real_width) { untested();
    throw exception_invalid_td(E_INV_SOLN_BAGSIZE);
  }
}

/*
 * Given a graph and a decomposition, this checks whether or not the set of
 * vertices in the graph equals the union of all bags in the decomposition.
 */
bool check_vertex_coverage(tree_decomposition& T)
{
  if (!(n_graph == n_vertices)) { untested();
    std::cerr << "Error: .gr and .td disagree on how many vertices the graph has" << std::endl;
    return false;
  } else if (n_vertices == 0) return true;

  std::vector<unsigned> occurrence_nums(n_graph, 0);
  for (unsigned i = 0; i < n_decomp; i++) {
    for (auto it = T.get_bag(i).begin(); it != T.get_bag(i).end(); it++) {
      occurrence_nums[*it - 1]++;
    }
  }

  for (unsigned i = 0; i < occurrence_nums.size(); i++) {
    if (occurrence_nums[i] == 0) { untested();
      std::cerr << "Error: vertex " << (i+1) << " appears in no bag" << std::endl;
      return false;
    }
  }
  return true;
}

/*
 * Given a graph and a decomposition, this checks whether or not each edge is
 * contained in at least one bag.  This has the side effect of removing from the
 * graph all those edges of the graph that are in fact contained in some bag.
 * The emptiness of the edge set of the resulting pruned graph is used to decide
 * whether the decomposition has this property.
 */
bool check_edge_coverage(tree_decomposition& T, graph& g)
{
  std::vector<std::pair<unsigned,unsigned> > to_remove;
  /*
   * We go through all bags, and for each vertex in each bag, we remove all its
   * incident edges.  If all the edges are indeed covered by at least one bag,
   * this will leave the graph with an empty edge-set, and it won't if they
   * aren't.
   */
  for(unsigned i = 0; i < T.bags.size() && g.num_edges > 0; i++){
    std::set<unsigned>& it_bag = T.get_bag(i);
    for (std::set<unsigned>::iterator head = it_bag.begin(); head != it_bag.end(); head++) {
      for(auto tail = g.neighbors(*head-1).begin(); tail != g.neighbors(*head-1).end(); tail++) {
        if(it_bag.find(*tail+1) != it_bag.end()) {
          to_remove.push_back(std::pair<unsigned, unsigned>(*head,*tail));
        }
      }
      for (std::vector<std::pair<unsigned, unsigned> >::iterator rem_it = to_remove.begin(); rem_it != to_remove.end(); rem_it++) {
        g.remove_edge(rem_it->first-1,rem_it->second);
      }
      to_remove.clear();
    }
  }

  if (g.num_edges > 0)
  { untested();
    for (unsigned u = 0; u < g.num_vertices; u++)
    { untested();
      if (! g.adj_list.at(u).empty() )
      { untested();
        unsigned v=g.adj_list.at(u).front();
        std::cerr << "Error: edge {"<< (u+1) << ", " << (v+1) << "} appears in no bag" << std::endl;
        break;
      }
    }
  }
  return (g.num_edges == 0);
}

/*
 * Given a graph and a decomposition, this checks whether or not the set of bags
 * in which a given vertex from the graph appears forms a subtree of the tree
 * decomposition.  It does so by successively removing leaves from the
 * decomposition until the tree is empty, and hence the decomposition will
 * consist only of isolated empty bags after calling this function.
 */
bool check_connectedness(tree_decomposition& T)
{
  if (T.bags.size() == 0) return true;
  /*
   * At each leaf, we first check whether it contains some forgotten vertex (in
   * which case we return false), and if it doesn't, we compute whether there
   * are any vertices that appear in the leaf, but not in its parent, and if so,
   * we add those vertices to the set of forgotten vertices (Now that I come to
   * think of it, 'forbidden' might be a better term TODO.) We then remove the
   * leaf and continue with the next leaf.  We return true if we never encounter
   * a forgotten vertex before the tree is entirely deleted.
   *
   * It is quite easy to see that a tree decomposition satisfies this property
   * (given that it satisfies the others) if and only if it satisfies the
   * condition of the bags containing any given vertex forming a subtree of the
   * decomposition.
   */

  // forgotten[i] != 0 if and only if the vertex with index i (starting from 0)
  // is currently forbidden
  std::vector<bool> forgotten(n_graph,0);
  std::stack<unsigned> parent_stack;
  unsigned NIL = T.bags.size()+1;
  unsigned next = NIL;
  unsigned head = NIL;
  std::stack<unsigned> stack;
  stack.push(0);
  parent_stack.push(NIL);
  while (!stack.empty()) {
    head = next;
    next = stack.top();

    // Once all subtrees at next are processed, do the operation at next (i.e.,
    // this is post-order traversal)
    if (head == next) {
      stack.pop();
      if (head==parent_stack.top()) parent_stack.pop();
      unsigned parent = parent_stack.top();
      for (auto it = T.get_bag(head).begin(); it != T.get_bag(head).end(); ++it) {
        if (forgotten[*it-1] == true) return false;
        if (parent != NIL && T.get_bag(parent).find(*it) == T.get_bag(parent).end()) {
          forgotten[*it-1] = true;
        }
      }
      T.remove_vertex(head);
    } else {
      parent_stack.push(next);
      for (unsigned i = 0; i < T.neighbors(next).size(); i++) {
        if (T.neighbors(next).at(i) != head) stack.push(T.neighbors(next).at(i));
      }
    }
  }
  return true;
}

/*
 * Using all of the above functions, this simply checks whether a given
 * purported tree decomposition for some graph actually is one.
 */
bool is_valid_decomposition(tree_decomposition& T, graph& g)
{
  if (!T.is_tree()) { untested();
    std::cerr << "Error: not a tree" << std::endl;
    return false;
  } else if (!check_vertex_coverage(T)) { untested();
    return false;
  } else if (!check_edge_coverage(T,g)) { untested();
    return false;
  } else if (!check_connectedness(T)) { untested();
    std::cerr << "Error: some vertex induces disconnected components in the tree" << std::endl;
    return false;
  }
  else {
    return true;
  }
}

int main(int argc, char** argv) {
  bool is_valid = true;
  bool empty_td_file = false;

  if (argc < 2 || argc > 3) { untested();
    std::cerr << "Usage: " << argv[0] << " input.gr [input.td]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Validates syntactical and semantical correctness of .gr and .td files. If only a .gr file is given, it validates the syntactical correctness of that file." << std::endl;
    exit(1);
  }
  tree_decomposition T;
  graph g;

  std::ifstream fin(argv[1]);
  try {
    read_graph(fin,g);
  } catch (const std::invalid_argument& e) { untested();
    std::cerr << "Invalid format in " << argv[1] << ": " << e.what() << std::endl;
    is_valid = false;
  }
  fin.close();

  if(argc==3 && is_valid) {
    fin.open(argv[2]);
    if (fin.peek() == std::ifstream::traits_type::eof()) { untested();
      is_valid = false;
      empty_td_file = true;
    } else {
      try {
        read_tree_decomposition(fin, T);
      } catch (const std::invalid_argument& e) { untested();
        std::cerr << "Invalid format in " << argv[2] << ": " << e.what() << std::endl;
        is_valid = false;
      }
    }
    fin.close();
    if (is_valid) {
      is_valid = is_valid_decomposition(T,g);
    }
  }

  if (is_valid) {
    std::cerr << "valid" << std::endl;
    return 0;
  } else if (empty_td_file) { untested();
    std::cerr << "invalid: empty .td file" << std::endl;
    return 2;
  } else { untested();
    std::cerr << "invalid" << std::endl;
    return 1;
  }
}
