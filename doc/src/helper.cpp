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
// Offers some helper functions e.g. reading/writing graphs/decompositions from/to the dot format
//

#ifndef TD_HELPER
#define TD_HELPER

#include <sys/time.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>
#include <boost/algorithm/string/replace.hpp>

double time_diff(struct timeval x, struct timeval y){
    double x_ms, y_ms, diff;
    
    x_ms = (double)x.tv_sec*1000000 + (double)x.tv_usec;
    y_ms = (double)y.tv_sec*1000000 + (double)y.tv_usec;
    
    diff = (double)y_ms - (double)x_ms;
    
    return diff/1000;
}

template <typename G_t>
void read_dot_graph_directly(const std::string &src, G_t &G){
    G.clear();

    boost::dynamic_properties dp(boost::ignore_other_properties); 
    std::ifstream dot_graph(src.c_str());
    read_graphviz(dot_graph, G, dp);
}

template <typename G_t>
void read_dot_graph(const std::string &src, G_t &G){
    G.clear();
    //try to read as a directed graph
    try{
        typename boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> dotG;

        boost::dynamic_properties dp(boost::ignore_other_properties); 

        std::ifstream dot_graph(src.c_str());
        read_graphviz(dot_graph, dotG, dp);
    
        boost::copy_graph(dotG, G);
    }
    //try to read as a undirected graph
    catch(...){
        typename boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> dotG;

        boost::dynamic_properties dp(boost::ignore_other_properties); 

        std::ifstream dot_graph(src.c_str());
        read_graphviz(dot_graph, dotG, dp);
    
        boost::copy_graph(dotG, G);
    }
}


template <typename G_t>
void write_dot_graph(const std::string &dest, G_t &G){
    std::ofstream dot_graph(dest.c_str());
    
    boost::dynamic_properties dp(boost::ignore_other_properties);

    boost::write_graphviz_dp(dot_graph, G, dp);
}

template <typename T_t>
void write_dot_td(const std::string &dest, T_t &T){
    std::ofstream dot_graph(dest.c_str());
    
    std::string *name = new std::string[boost::num_vertices(T)];
    for (unsigned int i = 0; i < boost::num_vertices(T); i++){
        std::ostringstream os;
        std::set<unsigned int>::const_iterator v1;
        for (v1 = T[i].bag.begin(); v1 != T[i].bag.end(); ++v1)
            os << *v1 << " ";

        name[i] = os.str();
    }

    write_graphviz(dot_graph, T, boost::make_label_writer(name));
    delete[] name;
}

#endif
