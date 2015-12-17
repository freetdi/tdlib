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
// Offers functionality to compute a minimal seperator of two vertex sets
//
// These functions are most likely to be interesting for outside use:
//
// - void seperate_vertices(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &Y, std::set<unsigned int> &S)
//
//
// computes a seperator S up to size k, aborts and returns false if S would be greater than k:
//
// - bool seperate_vertices(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &Y, std::set<unsigned int> &S, unsigned int k)
//

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include "TD_simple_graph_algos.hpp"

#ifndef TD_STRUCT_VERTEX_VI_PD
#define TD_STRUCT_VERTEX_VI_PD

struct Vertex_VI_PD{
    unsigned int id;
    //bool visited;
    //int predecessor;
};

#endif

#ifndef TD_TYPEDEF_DIGRAPH
#define TD_TYPEDEF_DIGRAPH

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex_VI_PD> digraph_t;
typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, Vertex_VI_PD> H_t;

#endif

#ifndef TD_NETWORKFLOW
#define TD_NETWORKFLOW

//Makes G a digraph, resulting in the graph H. Adds a source to H, connected with vertices in X, that are not on paths in P
template <typename G_t>
digraph_t::vertex_descriptor make_digraph(G_t &G, std::vector<bool> &disabled, digraph_t &H, std::vector<digraph_t::vertex_descriptor> &idxMap, unsigned int max){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;

    unsigned int i = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(!disabled[G[*vIt].id]){
            idxMap[G[*vIt].id] = boost::add_vertex(H);
            H[i++].id = G[*vIt].id;
        }
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        if(!disabled[G[boost::source(*eIt, G)].id] && !disabled[G[boost::target(*eIt, G)].id]){
            boost::add_edge(idxMap[G[boost::source(*eIt, G)].id], idxMap[G[boost::target(*eIt, G)].id], H);
            boost::add_edge(idxMap[G[boost::target(*eIt, G)].id], idxMap[G[boost::source(*eIt, G)].id], H);
        }
    }


    digraph_t::vertex_descriptor source = boost::add_vertex(H);
    H[source].id = max+1;
    idxMap[max+1] = source;

    return source;
}

static void modify_digraph(digraph_t &H, std::set<unsigned int> &X, std::vector<std::vector<std::vector<unsigned int> > > &P, std::vector<digraph_t::vertex_descriptor> &idxMap, digraph_t::vertex_descriptor source){
    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++)
        boost::add_edge(source, idxMap[*sIt], H);

    //remove the edges on paths
    for(unsigned int i = 0; i < P.size(); i++){
        for(unsigned int j = 0; j < P[i].size(); j++)
            boost::remove_edge(idxMap[P[i][j][0]], idxMap[P[i][j][1]], H);
    }

    //remove the edges from the source to vertices in X, that are on a path
    for(unsigned int i = 0; i < P.size(); i++)
        boost::remove_edge(source, idxMap[P[i][0][0]], H);
}

//follows the edge set of H beginning on the start_vertex (a vertex in Y), aborting, if a vertex in X is reached, saving the path in new_path
static void follow_path(H_t &H, boost::graph_traits<H_t>::vertex_descriptor v, std::vector<std::vector<unsigned int> > &new_path, std::set<unsigned int> &X){
    boost::graph_traits<H_t>::adjacency_iterator nIt, nEnd;
    while(true){
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, H);
        if(nIt == nEnd)
            break;
        std::vector<unsigned int> edge;
        edge.push_back(H[*nIt].id);
        edge.push_back(H[v].id);
        new_path.push_back(edge);
        boost::remove_edge(v, *nIt, H);
        v = *nIt;
        if(X.find(H[*nIt].id) != X.end())
            break;
    }
    std::reverse(new_path.begin(), new_path.end());
}

//computes the symmetric difference between edges on the walk and on paths in P
//inserts the edges in the symmetric difference into an empty graph and follows the paths starting in Y
//function is called, if there exists a P-alternating walk
static void modify_paths(std::set<unsigned int> &X, std::set<unsigned int> Y,
                  std::vector<std::vector<unsigned int> > &walk, std::vector<std::vector<std::vector<unsigned int> > > &P, unsigned int source){
    //unfold the edges on P
    std::vector<std::vector<unsigned int> > edgesP;
    for(unsigned int i = 0; i < P.size(); i++){
        for(unsigned int j = 0; j < P[i].size(); j++)
            edgesP.push_back(P[i][j]);
    }

    //the new edge set is the symmetric difference of edgesP and walk
    std::vector<std::vector<unsigned int> > symDiff;
    if(edgesP.size() == 0)
        symDiff = walk;
    else{
        for(unsigned int i = 0; i < edgesP.size(); i++){
            bool contained = false;
            for(unsigned int j = 0; j < walk.size(); j++){
                if(edgesP[i][0] == walk[j][1] && edgesP[i][1] == walk[j][0]){
                    contained = true;
                    break;
                }
            }
            if(!contained)
                symDiff.push_back(edgesP[i]);
        }
        for(unsigned int i = 0; i < walk.size(); i++){
            bool contained = false;
            for(unsigned int j = 0; j < edgesP.size(); j++){
                if(edgesP[j][0] == walk[i][1] && edgesP[j][1] == walk[i][0]){
                    contained = true;
                    break;
                }
            }
            if(!contained)
                symDiff.push_back(walk[i]);
        }
    }

    unsigned int P_size = P.size();
    P.clear();
    P.resize(P_size+1);

    //there exist |P| + 1 disjoint paths, coded in the edges "symDiff".
    //add the edge-set "symDiff" to an empty graph and follow
    //the edges starting on a vertex contained in Y until a vertex
    //in X is reached. This will be a new path. Construct all such paths.
    H_t H;

    boost::graph_traits<H_t>::vertex_descriptor hdesc;
    std::set<unsigned int> added_vertices;

    unsigned int max = 0;

    for(unsigned int i = 0; i < symDiff.size(); i++){
        max = (symDiff[i][0] > max)? symDiff[i][0] : max;
        max = (symDiff[i][1] > max)? symDiff[i][1] : max;
    }

    std::vector<boost::graph_traits<H_t>::vertex_descriptor> vertexMap(max+1);

    for(unsigned int i = 0; i < symDiff.size(); i++){
        if(added_vertices.find(symDiff[i][0]) == added_vertices.end()){
            hdesc = boost::add_vertex(H);
            H[hdesc].id = symDiff[i][0];
            added_vertices.insert(symDiff[i][0]);
            vertexMap[symDiff[i][0]] = hdesc;
        }
        if(added_vertices.find(symDiff[i][1]) == added_vertices.end()){
            hdesc = boost::add_vertex(H);
            H[hdesc].id = symDiff[i][1];
            added_vertices.insert(symDiff[i][1]);
            vertexMap[symDiff[i][1]] = hdesc;
        }
    }
    for(unsigned int i = 0; i < symDiff.size(); i++)
        boost::add_edge(vertexMap[symDiff[i][1]], vertexMap[symDiff[i][0]], H);

    //follow the paths starting in Y until a sink is reached (follow_path)
    boost::graph_traits<H_t>::vertex_iterator hIt, hEnd;

    for(unsigned int i = 0; i < P_size+1; i++){
        std::vector<std::vector<unsigned int> > new_path;
        for(boost::tie(hIt, hEnd) = boost::vertices(H); hIt != hEnd; hIt++){
            std::set<unsigned int>::iterator yIt = Y.find(H[*hIt].id);
            if(yIt != Y.end()){
                Y.erase(yIt);
                follow_path(H, *hIt, new_path, X);
                break;
            }
        }

        P[i] = new_path;
    }

    walk.clear();
}

//takes the last vertex on each path, that could be reached by a P-alternating walk, that does not end in a vertex in Y
//if no such one exists, the first vertex on a path is taken
static void calculate_seperator(std::vector<std::vector<std::vector<unsigned int> > > &P, std::set<unsigned int> &cS, std::set<unsigned int> &S){
    for(unsigned int i = 0; i < P.size(); i++){
        for(unsigned int j = P[i].size(); j > 0; j--){
            if(cS.find(P[i][j-1][1]) != cS.end()){
                S.insert(P[i][j-1][1]);
                break;
            }
            else if(j == 1)
                S.insert(P[i][j-1][0]);
        }
    }
}

//collects the neighbourhood of v and modifies the neighbourhood if v is on a path in P
static bool t_search_disjoint_ways(digraph_t &diG, digraph_t::vertex_descriptor v, std::set<unsigned int> &Y, std::vector<boost::tuple<bool, int> > &visited, 
   std::vector<std::vector<unsigned int> > &walk, std::vector<digraph_t::vertex_descriptor> &idxMap, bool edge_used, int source, std::set<unsigned int> &dangerous,
    std::set<unsigned int> &cS
){
    visited[diG[v].id].get<0>() = true;
    bool on_a_path = visited[diG[v].id].get<1>() != -1;

    //the walk has reached a valid vertex in Y
    if(Y.find(diG[v].id) != Y.end() && !on_a_path)
        return true;

    std::set<unsigned int> neighbours;

    //v is on a path, save this vertex for computing a seperator
    if(on_a_path)
        cS.insert(diG[v].id);

    //case, that v is on a path in P and the last visited vertex was not the successor of v on this path in P
    //this vertex could be possibly reached by the successor of v on the path at a later time, without the condition
    //to use the edge on the path in the opposite direction. For this reason, v is stored in dangerous
    if(on_a_path && !edge_used){
        if(visited[diG[v].id].get<1>() != source)
            neighbours.insert(visited[diG[v].id].get<1>());

        dangerous.insert(diG[v].id);
        visited[diG[v].id].get<0>() = false;
    }
    else{
        digraph_t::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, diG); nIt != nEnd; nIt++)
            neighbours.insert(diG[*nIt].id);
    }

    for(std::set<unsigned int>::iterator sIt = neighbours.begin(); sIt != neighbours.end(); sIt++){
        if(!visited[*sIt].get<0>()){
            std::set<unsigned int>::iterator it = dangerous.find(*sIt);
            //do not visit a vertex, if it is dangerous and we do not come from its successor on the path
            if(it != dangerous.end() && (int)*it != visited[diG[v].id].get<1>())
                continue;
            //visit a vertex, if it is dangerous and we come from its successor on the path, the vertex is no more dangerous
            //and will be marked as visited
            if((int)*it == visited[diG[v].id].get<1>())
                dangerous.erase(*it);
            //we would use an edge of a path if *sIt is the successor of v on some path
            bool edge_used_ = visited[diG[v].id].get<1>() == (int)*sIt;

            //recursivly builds the walk, if a valid vertex in Y could be reached
            if(t_search_disjoint_ways(diG, idxMap[*sIt], Y, visited, walk, idxMap, edge_used_, source, dangerous, cS)){
                std::vector<unsigned int> edge;
                edge.push_back(diG[v].id);
                edge.push_back(*sIt);
                walk.push_back(edge);
                return true;
            }
        }
    }
    return false;
}

template <typename G_t>
bool _disjoint_ways(G_t &G, std::vector<bool> &disabled, std::set<unsigned int> &X, std::set<unsigned int> &Y, std::set<unsigned int> &S, unsigned int k){ 
    //evaluate maximum id
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(!disabled[G[*vIt].id])
            max = (G[*vIt].id > max)? G[*vIt].id : max;
    }

    std::vector<boost::tuple<bool, int> > visited(max+2);

    //make G a digraph, add a new source, connect the source with vertices in X
    digraph_t diG;
    std::vector<digraph_t::vertex_descriptor> idxMap(max+2);
    digraph_t::vertex_descriptor source = make_digraph(G, disabled, diG, idxMap, max);

    std::vector<std::vector<unsigned int> > walk;
    std::vector<std::vector<std::vector<unsigned int> > > P;

    std::set<unsigned int> cS;

    //main loop of algorithm
    //since at most one path can begin in a distinct vertex in X, |X| + 1 iterations are sufficient (one for the unavailing try) 
    for(unsigned int iter = 0; iter < X.size()+1; iter++){
        if(S.size()+iter > k)
            return false;

        cS.clear();

        digraph_t diG_;
        boost::copy_graph(diG, diG_);
        modify_digraph(diG_, X, P, idxMap, source);

        for(unsigned int i = 0; i < visited.size(); i++){
            visited[i].get<0>() = false;
            visited[i].get<1>() = -1;
        }

        for(unsigned int i = 0; i < P.size(); i++){
            for(unsigned int j = 0; j < P[i].size(); j++)
                visited[P[i][j][1]].get<1>() = P[i][j][0];
        }
        for(unsigned int i = 0; i < P.size(); i++)
            visited[P[i][0][0]].get<1>() = max+1;

        //start extended DFS at source
        std::set<unsigned int> dangerous;
        if(!t_search_disjoint_ways(diG_, source, Y, visited, walk, idxMap, false, (int)diG[source].id, dangerous, cS))
            break;

        modify_paths(X, Y, walk, P, max+1);
    }
    calculate_seperator(P, cS, S);

    return true;
}

template <typename G_t>
bool seperate_vertices(G_t &G, std::vector<bool> &disabled, std::set<unsigned int> &X, std::set<unsigned int> &Y, std::set<unsigned int> &S, unsigned int k){ 
    //common neighbours must be contained in a seperator
    std::set_intersection(X.begin(), X.end(), Y.begin(), Y.end(), std::inserter(S, S.begin()));

    for(std::set<unsigned int>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
        X.erase(*sIt);

    for(std::set<unsigned int>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
        Y.erase(*sIt);

    if(S.size() > k)
        return false;

    if(X.size() == 0 || Y.size() == 0)
        return true;

    for(std::set<unsigned int>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
        disabled[*sIt] = true;

    return _disjoint_ways(G, disabled, X, Y, S, k);
}

template <typename G_t>
void seperate_vertices(G_t &G, std::vector<bool> &disabled, std::set<unsigned int> &X, std::set<unsigned int> &Y, std::set<unsigned int> &S){  
    seperate_vertices(G, disabled, X, Y, S, UINT_MAX);
}

#endif
