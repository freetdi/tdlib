#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include "Thorup.hpp"

#include <tdlib/combinations.hpp>
#include <tdlib/misc.hpp>

#ifdef RESULTS_HERE
#define RESULTS_DIR
#else
#define RESULTS_DIR "./Results/"
#endif

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, treedec::bag_t> TD_tree_dec_t;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Thorup_graph_t;

template <typename Thorup_G_t, typename G_t>
void thorup_graph(Thorup_G_t &G2, G_t &G1){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<typename boost::graph_traits<Thorup_G_t>::vertex_descriptor> G2vertices;
    unsigned int currentID = 0;
    while (currentID != boost::num_vertices(G1)){
        for(boost::tie(vIt, vEnd) = boost::vertices(G1); vIt != vEnd; vIt++){
            if((unsigned int) *vIt == currentID){
                G2vertices.push_back(boost::add_vertex(G2));
                currentID++;
                break;
            }
        }
    }
    for(boost::tie(vIt, vEnd) = boost::vertices(G1); vIt != vEnd; vIt++){
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G1); nIt != nEnd; nIt++){
            if(!boost::edge(G2vertices[*vIt], G2vertices[*nIt], G2).second)
                boost::add_edge(G2vertices[*vIt], G2vertices[*nIt], G2);
        }
    }
}


double totaltime = 0;
int max_width = -2;

unsigned graphs_count = 0;
double avg_width = 0.0;

void test_single(std::string filename){
    struct timeval t1, t2;
    double time1;

    TD_graph_t G;
    TD_graph_t H;

    read_dot_graph(filename, G);
    read_dot_graph(filename, H);

    TD_tree_dec_t T;

    std::cout << filename << " [" << boost::num_vertices(G) << ";" << boost::num_edges(G) << "]" << std::endl;

    Thorup_graph_t TG;
    thorup_graph(TG, G);

    gettimeofday(&t1, NULL);
    thorup_tree_decomposition(T, TG);
    gettimeofday(&t2, NULL);
    time1 = time_diff(t1, t2);

    graphs_count++;

    int width = treedec::get_width(T);

    avg_width += width;

    max_width = (width > max_width)? width : max_width;

    if(!treedec::is_valid_treedecomposition(H, T)){
        std::cout << "invalid tree decomposition: " << filename << std::endl;
        write_dot_td("error_td.dot", T);
        exit(false);
    }

    totaltime += time1;

}

int main(int argc, char * const * argv){
    std::cout.setf( std::ios_base::unitbuf );

    if(argc<3){
        std::cerr << "usage\n"
                  << "    " << argv[0] << " [package] [*.dot]" << std::endl;
        exit(1);
    }

    std::string package(argv[1]);

    for(int i = 2; i < argc; i++){
        std::ifstream fin(argv[i]);
        if(!fin.is_open()){
            std::cerr << "error: could not open file " << argv[i] << std::endl;
        }
        fin.close();

        test_single(std::string(argv[i]));
    }

    avg_width /= graphs_count;


    std::ofstream texresults(RESULTS_DIR "13_thorup.tex", std::ios::app);

    tex_tabular_entry(texresults, package, graphs_count, avg_width, max_width, totaltime);

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

