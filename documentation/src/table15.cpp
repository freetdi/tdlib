#include <fstream>
#include <string>
#include <set>
#include <boost/algorithm/string/replace.hpp>

#include "tex.cpp"
#include "helper.cpp"

#include <tdlib/postprocessing.hpp>
#include <tdlib/misc.hpp>

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, treedec::bag_t> TD_tree_dec_t;


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

    gettimeofday(&t1, NULL);
    treedec::trivial_decomposition(G, T);
    treedec::MSVS(G, T);
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


    std::ofstream texresults("./Results/15_msvs_trivial.tex", std::ios::app);

    tex_tabular_entry(texresults, package, graphs_count, avg_width, max_width, totaltime);

    texresults.close();

    std::cout << "totaltime: " << totaltime << std::endl;
}

